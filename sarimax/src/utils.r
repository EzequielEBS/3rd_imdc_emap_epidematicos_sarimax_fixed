library(forecast)
library(dplyr)
library(tidyverse)
library(tidyr)
library(purrr)
library(lubridate)
library(scoringutils)
library(parallel)
library(pbapply)

source("sarimax/src/forecast_covariates.r")

#' Convert a YYYYWW epiweek integer to the Date of its opening Sunday (MMWR)
#'
#' Anchors on Jan 4, which is always in MMWR epidemiological week 1.
#'
#' @param yw  Integer (or vector) epiweek in YYYYWW format, e.g. 202541.
#'
#' @return Date (or vector of Dates): the Sunday opening that epiweek.
epiweek_to_date <- function(yw) {
  yr   <- yw %/% 100L
  wk   <- yw  %% 100L
  jan4 <- as.Date(paste0(yr, "-01-04"))
  dow_jan4  <- as.integer(format(jan4, "%w"))  # %w: 0 = Sunday
  sunday_w1 <- jan4 - dow_jan4                 # Sunday on or before Jan 4
  sunday_w1 + (wk - 1L) * 7L
}

#' Check whether a given year has 53 epiweeks (MMWR/Brazilian calendar)
#'
#' Derived from `epiweek_to_date()` for consistency: a year has 53 weeks iff
#' the Sunday that would open week 53 falls before week 1 of the next year.
#'
#' @param year  Integer year (e.g. 2025).
#'
#' @return Logical: TRUE if `year` has 53 epiweeks, FALSE if it has 52.
has_53_weeks <- function(year) {
  epiweek_to_date(year * 100L + 53L) < epiweek_to_date((year + 1L) * 100L + 1L)
}

#' Enumerate all epiweeks (inclusive) between two YYYYWW integers
#'
#' Correctly rolls over year boundaries and accounts for 52- vs 53-week years
#' via `has_53_weeks()`.
#'
#' @param start  Integer epiweek in YYYYWW format marking the start (inclusive).
#' @param end    Integer epiweek in YYYYWW format marking the end (inclusive).
#'
#' @return Integer vector of YYYYWW epiweeks from `start` to `end`.
enumerate_epiweeks <- function(start, end) {
  start_year <- start %/% 100L
  start_week <- start  %% 100L
  end_year   <- end   %/% 100L
  end_week   <- end    %% 100L

  epiweeks <- integer(0)
  yr <- start_year
  wk <- start_week

  repeat {
    epiweeks <- c(epiweeks, yr * 100L + wk)
    if (yr == end_year && wk == end_week) break
    n_weeks <- if (has_53_weeks(yr)) 53L else 52L
    if (wk < n_weeks) {
      wk <- wk + 1L
    } else {
      yr <- yr + 1L
      wk <- 1L
    }
  }
  epiweeks
}

#' Load one unit's (state or municipality) processed data file, normalized to
#' a single canonical column schema
#'
#' State-level aggregate files (`processed_data/<disease>/<disease>_<UF>_agg.csv.gz`)
#' and municipality-level files (`processed_data/<disease>/sel_cities/<disease>_<geocode>.csv.gz`)
#' come from different upstream pipelines and don't share column names: the
#' response is `cases` at state level but `casos` (Portuguese) at municipality
#' level, and municipality weather variables (`temp_min`, `precip_med`, ...)
#' lack the `_mean` suffix state-level variables carry for the same
#' already-aggregated quantity (`temp_min_mean`, `precip_med_mean`, ...).
#' Every other function in this pipeline (`get_candidates()`, `fit_sarimax()`,
#' `pca_all()`, ...) is written against the state-level names, so rather than
#' teaching each of them about two schemas, this function normalizes once at
#' load time: after this, `cases` and the `_mean`-suffixed weather names are
#' always present regardless of level, and every downstream function stays
#' completely level-agnostic.
#'
#' Note: municipality files carry no lagged (`_lag<k>`) or rolling-window
#' (`_mean_<k>mo`) covariates at all -- that feature engineering has only
#' been done at state level upstream -- so `get_candidates()` on
#' municipality data will only ever return contemporaneous candidates (plus
#' any of enso/iod/pdo that survive filtering), never lagged/rolling ones.
#' That's a real difference in available signal, not a bug this function
#' works around.
#'
#' @param disease  "dengue" or "chikungunya".
#' @param level    "state" or "city".
#' @param id       Two-letter UF code (level = "state") or IBGE municipality
#'                 geocode (level = "city").
#' @param processed_data_dir  Base directory containing the disease
#'                 subfolders (default "processed_data").
#'
#' @return A tibble with canonical column names (`cases`, `epiweek`,
#'   `train_1..4`/`target_1..4`, `<var>_mean` weather columns, etc.).
load_unit_data <- function(disease, level = c("state", "city"), id,
                            processed_data_dir = "processed_data") {
  level <- match.arg(level)

  file_name <- if (level == "state") {
    file.path(processed_data_dir, disease, paste0(disease, "_", id, "_agg.csv.gz"))
  } else {
    file.path(processed_data_dir, disease, "sel_cities", paste0(disease, "_", id, ".csv.gz"))
  }
  if (!file.exists(file_name)) {
    stop("load_unit_data(): file not found: ", file_name)
  }
  data <- read_csv(file_name, show_col_types = FALSE)

  # -- Response column: municipality files use the Portuguese "casos" -------
  if (!"cases" %in% names(data) && "casos" %in% names(data)) {
    data <- dplyr::rename(data, cases = casos)
  }

  # -- Weather variables: municipality files omit the "_mean" suffix that ---
  #    state-level files use for the same (already-aggregated) quantity.
  weather_bases <- c("temp_min", "temp_med", "temp_max",
                      "precip_min", "precip_med", "precip_max",
                      "pressure_min", "pressure_med", "pressure_max",
                      "rel_humid_min", "rel_humid_med", "rel_humid_max",
                      "thermal_range", "rainy_days")
  for (base in weather_bases) {
    mean_name <- paste0(base, "_mean")
    if (base %in% names(data) && !(mean_name %in% names(data))) {
      names(data)[names(data) == base] <- mean_name
    }
  }

  # -- Population: cosmetic alignment only; not a default candidate column --
  if (!"pop" %in% names(data) && "population" %in% names(data)) {
    data <- dplyr::rename(data, pop = population)
  }

  stopifnot(
    "epiweek" %in% names(data),
    "cases"   %in% names(data),
    all(paste0("train_",  1:4) %in% names(data)),
    all(paste0("target_", 1:4) %in% names(data))
  )

  data
}

#' Extract right-hand-side term names from a SARIMAX covariate formula spec
#'
#' Accepts the shapes used across this repo: NULL/empty (no covariates), a
#' `formula` object, or a character string/vector of term names (optionally
#' "+"-joined, optionally with a leading "~").
#'
#' @param formula  NULL, a formula, or a character string/vector.
#'
#' @return Character vector of term names (possibly empty).
parse_rhs_terms <- function(formula) {
  if (is.null(formula)) return(character(0))
  if (inherits(formula, "formula")) {
    return(attr(stats::terms(formula), "term.labels"))
  }
  if (is.character(formula)) {
    if (length(formula) == 0) return(character(0))
    rhs <- paste(formula, collapse = "+")
    if (nchar(trimws(rhs)) == 0) return(character(0))
    fo <- stats::as.formula(paste("~", rhs))
    return(attr(stats::terms(fo), "term.labels"))
  }
  stop("parse_rhs_terms(): formula must be NULL, a formula, or a character string/vector")
}

#' Fit a SARIMAX model over an explicit epiweek range and forecast forward
#'
#' Fits `Arima()` on `data[epiweek >= train_start & epiweek <= train_end, ]`
#' and forecasts every epiweek in `forecast_start..forecast_end`. Every
#' regressor the formula uses is itself forecast forward with its own
#' seasonal time-series model (`forecast_covariates()`,
#' sarimax/src/forecast_covariates.r), fit only on data through `train_end`
#' -- never real future values, and never another year's realized value --
#' so the covariate side of the forecast is exactly as ex-ante as the
#' case-count model itself (see root README.md, Section 6, "Data Usage
#' Restriction"). This single function now covers every use in this repo:
#' the four retrospective validation splits (train_1..train_4) and the real
#' submission / forecast-phase window alike -- they differ only in which
#' epiweek range is passed in, not in how the fit or the covariate
#' forecasting works.
#'
#' @param data                A data frame containing `epiweek`, `cases`,
#'                            and every column referenced in `formula`.
#' @param formula             NULL, a formula, or a character string/vector
#'                            of covariate names (RHS only; `cases` is the
#'                            implicit response) -- see `parse_rhs_terms()`.
#' @param train_start         Integer YYYYWW: first training epiweek
#'                            (inclusive).
#' @param train_end           Integer YYYYWW: last training epiweek
#'                            (inclusive) -- also the cutoff used to fit
#'                            each regressor's own forecasting model.
#' @param forecast_start      Integer YYYYWW: first forecast epiweek
#'                            (inclusive).
#' @param forecast_end        Integer YYYYWW: last forecast epiweek
#'                            (inclusive).
#' @param method              Estimation method passed to `Arima()`
#'                            (default "CSS-ML").
#' @param lambda              Box-Cox lambda; NULL (default) uses
#'                            log1p/expm1 instead.
#' @param optim.control       List passed to `Arima()`'s optimizer (default
#'                            list(maxit = 500)).
#' @param optim.method        Optimization method passed to `Arima()`
#'                            (default "BFGS").
#' @param levels              Numeric vector of prediction-interval coverage
#'                            levels (default c(50, 80, 90, 95)).
#' @param order               ARIMA (p,d,q) order (default c(1,1,1)).
#' @param seasonal            Seasonal ARIMA list, e.g.
#'                            list(order = c(1,0,1), period = 52).
#' @param bootstrap           Logical: use bootstrapped simulation paths for
#'                            prediction intervals instead of Gaussian
#'                            normal-theory intervals (default TRUE).
#' @param npaths              Number of bootstrap simulation paths (default
#'                            1000).
#' @param covariate_forecasts Optional numeric matrix (forecast horizon rows
#'                            x regressor columns, named) of already-computed
#'                            covariate forecasts to use instead of calling
#'                            `forecast_covariates()` fresh. Meant for hot
#'                            loops that call this function many times for
#'                            the same fold with different formulas/orders
#'                            (see `run_grid_search()`'s per-fold precompute
#'                            step) -- forecasting each candidate regressor
#'                            once per fold, rather than once per
#'                            (formula, order, fold) combination, is what
#'                            keeps that grid search computationally
#'                            tractable. Must have >= the columns named in
#'                            the formula and exactly `forecast_end -
#'                            forecast_start` (inclusive) rows. NULL
#'                            (default) computes forecasts on the fly --
#'                            the right choice for a one-off call, as in
#'                            sarimax/src/fit.r.
#' @param table_path          Optional: passed straight through to
#'                            `forecast_covariates()`'s own `table_path`
#'                            when `covariate_forecasts` isn't already
#'                            supplied -- looks each needed root series up
#'                            in the persistent covariate-forecast table
#'                            (see `covariate_forecast_table_path()`)
#'                            first, only fitting live (and saving the
#'                            result) on a miss. Ignored when
#'                            `covariate_forecasts` is supplied, since
#'                            there's nothing left to forecast. Default
#'                            `NULL` preserves the original always-fit-live
#'                            behavior exactly.
#' @param id                  Required when `table_path` is supplied: this
#'                            unit's id (UF code or municipality geocode).
#'
#' @return A tibble with columns `date`, `pred`, and `lower_*`/`upper_*` for
#'         each level in `levels`, with `warnings` and `fit` attributes
#'         attached.
fit_sarimax <- function(data,
                         formula,
                         train_start,
                         train_end,
                         forecast_start,
                         forecast_end,
                         method        = "CSS-ML",
                         lambda        = NULL,
                         optim.control = list(maxit = 500),
                         optim.method  = "BFGS",
                         levels        = c(50, 80, 90, 95),
                         order         = c(1, 1, 1),
                         seasonal      = list(order = c(1, 0, 1), period = 52),
                         bootstrap     = TRUE,
                         npaths        = 1000,
                         covariate_forecasts = NULL,
                         table_path    = NULL,
                         id            = NULL) {

  stopifnot(is.data.frame(data))
  stopifnot("epiweek" %in% names(data))
  stopifnot("cases"   %in% names(data))

  # ── 1. Training rows ────────────────────────────────────────────────────
  train_rows <- data[data$epiweek >= train_start & data$epiweek <= train_end, ]
  if (nrow(train_rows) == 0) {
    stop("No training rows found for epiweeks ", train_start, "-", train_end)
  }

  # ── 2. Enumerate forecast epiweeks & build forecast dates ───────────────
  fc_epiweeks    <- enumerate_epiweeks(forecast_start, forecast_end)
  n_fc           <- length(fc_epiweeks)
  forecast_dates <- epiweek_to_date(fc_epiweeks)  # already Sundays

  # ── 3. Log-transform response (log1p to handle zero counts) ─────────────
  y <- if (is.null(lambda)) log1p(train_rows$cases) else train_rows$cases

  # ── 4. Build & standardize regressor matrices ───────────────────────────
  rhs_terms <- parse_rhs_terms(formula)

  if (length(rhs_terms) == 0) {
    xreg_train  <- NULL
    xreg_future <- NULL
  } else {
    raw_train <- as.matrix(train_rows[, rhs_terms, drop = FALSE])

    # Compute mean/sd from training data only (no data leakage)
    col_means <- colMeans(raw_train, na.rm = TRUE)
    col_sds   <- apply(raw_train, 2, sd, na.rm = TRUE)
    col_sds[col_sds == 0] <- 1

    xreg_train <- scale(raw_train, center = col_means, scale = col_sds)

    # Forecast each regressor forward with its own seasonal model, fit only
    # on train_start..train_end -- unless a precomputed matrix was supplied
    # (see the @param doc above), in which case just look the columns up.
    raw_future <- if (!is.null(covariate_forecasts)) {
      missing_terms <- setdiff(rhs_terms, colnames(covariate_forecasts))
      if (length(missing_terms) > 0) {
        stop("fit_sarimax(): covariate_forecasts is missing column(s): ",
             paste(missing_terms, collapse = ", "))
      }
      if (nrow(covariate_forecasts) != n_fc) {
        stop("fit_sarimax(): covariate_forecasts has ", nrow(covariate_forecasts),
             " row(s) but the forecast horizon needs ", n_fc)
      }
      covariate_forecasts[, rhs_terms, drop = FALSE]
    } else {
      forecast_covariates(
        data        = data,
        rhs_terms   = rhs_terms,
        train_start = train_start,
        train_end   = train_end,
        h           = n_fc,
        table_path  = table_path,
        id          = id
      )
    }
    xreg_future <- scale(raw_future, center = col_means, scale = col_sds)
  }

  # ── 5. Fit SARIMAX ───────────────────────────────────────────────────────
  y_ts         <- ts(y, frequency = seasonal$period)
  fit_warnings <- character(0)

  # NOTE: Arima()/forecast() decide whether xreg was supplied by inspecting
  # the *unevaluated call expression*, not the argument's value -- passing
  # `xreg = xreg_train` when xreg_train is NULL still counts as "xreg used"
  # and later breaks forecast() with "object 'xreg_train' not found" (it
  # tries to re-evaluate that symbol from a different call frame). Building
  # the argument list with do.call() and omitting `xreg` entirely when there
  # are no regressors avoids this.
  arima_args <- list(
    y_ts,
    order         = order,
    seasonal      = seasonal,
    method        = method,
    lambda        = lambda,
    optim.control = optim.control,
    optim.method  = optim.method
  )
  if (!is.null(xreg_train)) arima_args$xreg <- xreg_train

  fit <- withCallingHandlers(
    do.call(Arima, arima_args),
    error = function(e) stop("Arima() failed: ", conditionMessage(e)),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # ── 6. Detect NaN standard errors ───────────────────────────────────────
  se_warnings <- character(0)
  ses <- withCallingHandlers(
    sqrt(diag(fit$var.coef)),
    warning = function(w) {
      se_warnings <<- c(se_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  nan_se_params <- names(which(is.nan(ses)))
  if (length(nan_se_params) > 0) {
    fit_warnings <- c(
      fit_warnings,
      se_warnings,
      paste0(
        "NaN standard errors for: ",
        paste(nan_se_params, collapse = ", "),
        ". Prediction intervals may be unreliable. ",
        "Consider reducing model order or using auto.arima()."
      )
    )
  }

  # ── 7. Forecast & back-transform ─────────────────────────────────────────
  levels      <- sort(unique(levels))
  fc_warnings <- character(0)

  fc_args <- list(
    fit,
    h         = n_fc,
    level     = levels,
    bootstrap = bootstrap,
    npaths    = npaths
  )
  if (!is.null(xreg_future)) fc_args$xreg <- xreg_future

  fc <- withCallingHandlers(
    do.call(forecast::forecast, fc_args),
    warning = function(w) {
      fc_warnings <<- c(fc_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  bt <- if (is.null(lambda)) expm1 else identity

  # ── 8. Assemble output tibble ─────────────────────────────────────────────
  out <- tibble::tibble(
    date = forecast_dates,
    pred = bt(as.numeric(fc$mean))
  )

  for (lv in levels) {
    lv_char <- paste0(lv, "%")
    out[[paste0("lower_", lv)]] <- bt(as.numeric(fc$lower[, lv_char]))
    out[[paste0("upper_", lv)]] <- bt(as.numeric(fc$upper[, lv_char]))
  }

  out <- out |>
    dplyr::mutate(dplyr::across(
      c(pred, dplyr::starts_with("lower_"), dplyr::starts_with("upper_")),
      \(x) pmax(x, 0)
    ))

  # ── 9. Attach metadata ────────────────────────────────────────────────────
  all_warnings <- c(fit_warnings, fc_warnings)
  attr(out, "warnings") <- if (length(all_warnings) > 0) all_warnings else NULL
  attr(out, "fit")      <- fit

  out
}


#' Build all covariate combinations where no pair exceeds a correlation threshold
#'
#' @param data        A data frame containing the candidate covariate columns.
#' @param covariates  Character vector of candidate covariate names.
#' @param threshold   Absolute Pearson correlation threshold (default 0.7).
#'                    Any pair with |r| >= threshold is considered collinear.
#' @param min_size    Minimum number of covariates per combination (default 1).
#' @param max_size    Maximum number of covariates per combination (default Inf).
#'
#' @return A list of character vectors, each a valid (uncorrelated) covariate set.
build_covariate_combinations <- function(data,
                                         covariates,
                                         threshold = 0.7,
                                         min_size  = 1,
                                         max_size  = Inf) {

  stopifnot(is.data.frame(data))
  stopifnot(all(covariates %in% names(data)))
  stopifnot(threshold > 0 & threshold <= 1)

  # ── 1. Correlation matrix ─────────────────────────────────────────────────
  cor_mat  <- cor(data[, covariates, drop = FALSE], use = "pairwise.complete.obs")
  collinear <- abs(cor_mat) >= threshold
  diag(collinear) <- FALSE

  # ── 2. Pre-compute total combinations for progress tracking ───────────────
  max_size  <- min(max_size, length(covariates))
  sizes     <- seq(min_size, max_size)
  combos_per_size <- sapply(sizes, \(k) choose(length(covariates), k))
  total     <- sum(combos_per_size)

  # ── 3. Progress bar setup ─────────────────────────────────────────────────
  start_time  <- proc.time()["elapsed"]
  checked     <- 0L
  valid       <- vector("list", 0)

  fmt_time <- function(secs) {
    if (secs < 60)  return(sprintf("%.0fs",       secs))
    if (secs < 3600) return(sprintf("%.0fm %.0fs", secs %/% 60, secs %% 60))
    sprintf("%.0fh %.0fm", secs %/% 3600, (secs %% 3600) %/% 60)
  }

  draw_progress <- function(checked, total, n_valid, start_time) {
    elapsed  <- proc.time()["elapsed"] - start_time
    pct      <- checked / total
    rate     <- if (elapsed > 0) checked / elapsed else 0
    eta      <- if (rate > 0) (total - checked) / rate else NA

    bar_width <- 30L
    filled    <- round(pct * bar_width)
    bar       <- paste0(
      "[", strrep("=", filled),
      if (filled < bar_width) ">" else "",
      strrep(" ", max(0L, bar_width - filled - 1L)),
      "]"
    )

    cat(sprintf(
      "\r%s %5.1f%%  checked: %d/%d  valid: %d  elapsed: %s  ETA: %s     ",
      bar, pct * 100, checked, total, n_valid,
      fmt_time(elapsed),
      if (is.na(eta)) "---" else fmt_time(eta)
    ))
    flush.console()
  }

  # ── 4. Iterate combinations with live progress ────────────────────────────
  for (k in sizes) {
    combos_k <- combn(covariates, k, simplify = FALSE)

    for (combo in combos_k) {
      sub_collinear <- collinear[combo, combo, drop = FALSE]
      if (!any(sub_collinear)) valid <- c(valid, list(combo))

      checked <- checked + 1L
      if (checked %% 500L == 0L || checked == total)
        draw_progress(checked, total, length(valid), start_time)
    }
  }

  # Final completed bar
  draw_progress(total, total, length(valid), start_time)
  cat("\n")

  elapsed_total <- proc.time()["elapsed"] - start_time
  message(sprintf(
    "Done. %d valid combination(s) from %d candidates | threshold = %.2f | time: %s",
    length(valid), length(covariates), threshold, fmt_time(elapsed_total)
  ))

  valid
}

#' Assemble candidate covariate names by category (contemporaneous, lagged, rolling)
#'
#' @param data         Data frame to search for matching column names.
#' @param groups       Character vector of categories to include: any of
#'                     "contemporaneous", "lagged", "rolling".
#' @param vars         Base contemporaneous variable names to look for verbatim.
#' @param lag_select   Optional integer vector restricting lagged columns to
#'                     specific lag values (matches "<var>_lag<k>"). If NULL,
#'                     all "_lag<digits>" columns are included.
#' @param roll_select  Optional character vector restricting rolling-window
#'                     columns to specific window labels (matches
#'                     "<var>_mean_<window>"). If NULL, all "_mean_<digits>mo"
#'                     columns are included.
#'
#' @return Character vector of unique column names found in `data` that match
#'         the requested groups/filters.
get_candidates <- function(data,
                           groups       = c("contemporaneous", "lagged", "rolling"),
                           vars         = c("temp_med_mean", "precip_med_mean",
                                            "rel_humid_med_mean", "thermal_range_mean",
                                            "rainy_days_mean", "enso", "iod", "pdo"),
                           lag_select   = NULL,
                           roll_select  = NULL
                          ) # single rolling window to keep
{
  all_cols <- names(data)
  selected <- character(0)

  if ("contemporaneous" %in% groups)
    selected <- c(selected, intersect(vars, all_cols))

  if ("lagged" %in% groups) {
    if (!is.null(lag_select)) {
      lag_pat  <- sprintf("_lag%d$", lag_select)
      selected <- c(selected, grep(lag_pat, all_cols, value = TRUE))
    } else {
      selected <- c(selected, grep("_lag\\d+$",    all_cols, value = TRUE))
    }
  }

  if ("rolling" %in% groups) {
    if (!is.null(roll_select)) {
      roll_pat <- sprintf("_mean_%s$", roll_select)
      selected <- c(selected, grep(roll_pat, all_cols, value = TRUE))
    } else {
      selected <- c(selected, grep("_mean_\\d+mo$", all_cols, value = TRUE))
    }
  }

  unique(selected)
}

#' Drop near-constant covariates (standard deviation below a threshold)
#'
#' @param data        Data frame containing the covariate columns.
#' @param covariates  Character vector of covariate names to check.
#' @param threshold   Minimum standard deviation required to keep a covariate
#'                    (default 0.01). Columns with sd <= threshold are dropped.
#'
#' @return Character vector of covariate names with low-variance columns removed.
filter_low_variance <- function(data, covariates, threshold = 0.01) {
  vars_sd <- sapply(covariates, \(v) sd(data[[v]], na.rm = TRUE))
  keep    <- names(vars_sd[vars_sd > threshold])
  dropped <- setdiff(covariates, keep)
  if (length(dropped) > 0)
    message("Dropped low-variance columns: ", paste(dropped, collapse = ", "))
  keep
}

#' Score a prediction tibble against observed values using scoringutils
#'
#' Converts the wide lower_*/pred/upper_* forecast tibble produced by
#' `fit_sarimax()` into the long quantile-forecast format `scoringutils`
#' expects, then computes the Weighted Interval Score plus point-forecast
#' error metrics and per-level empirical coverage.
#'
#' @param pred_df  Output tibble from `fit_sarimax()`: columns `pred`,
#'                 `lower_<level>`, `upper_<level>` for each level in `levels`.
#' @param actual   Numeric vector of observed case counts, same length/order
#'                 as `pred_df`.
#' @param levels   Numeric vector of prediction-interval coverage levels
#'                 present in `pred_df` (default c(50, 80, 90, 95)).
#'
#' @return A named list: `wis`, `mae`, `mse`, `rmse`, `mape`, and one
#'         `coverage_<level>` entry per level (empirical interval coverage).
# scoringutils expects a long data frame with quantile forecasts
compute_metrics <- function(
  pred_df,
  actual,
  levels = c(50, 80, 90, 95)
) {
  levels     <- sort(levels)
  lower_cols <- paste0("lower_", levels)
  upper_cols <- paste0("upper_", levels)
  alphas     <- 1 - levels / 100

  n <- nrow(pred_df)

  make_quantile_rows <- function(predicted, quantile_level) {
    tibble(
      time           = seq_len(n),
      observed       = actual,
      predicted      = predicted,
      quantile_level = quantile_level,
      model          = "sarimax"
    )
  }

  su <- bind_rows(
    # Lower bounds: alpha/2  (e.g. level=95 -> 0.025, level=50 -> 0.25)
    lapply(seq_along(levels), function(i)
      make_quantile_rows(pred_df[[lower_cols[i]]], alphas[i] / 2)),
    # Upper bounds: 1 - alpha/2
    lapply(seq_along(levels), function(i)
      make_quantile_rows(pred_df[[upper_cols[i]]], 1 - alphas[i] / 2)),
    # Median — always added; alpha/2 never equals 0.5 so no duplication risk
    make_quantile_rows(pred_df$pred, 0.5)
  ) |>
    arrange(time, quantile_level) |>
    as_forecast_quantile(
      observed       = "observed",
      predicted      = "predicted",
      quantile_level = "quantile_level",
      model          = "model"
    )

  sc <- withCallingHandlers(
    score(su) |> summarise_scores(by = "model"),
    warning = function(w) {
      if (grepl("interval_coverage", conditionMessage(w)))
        invokeRestart("muffleWarning")
    }
  )

  coverage <- setNames(
    sapply(seq_along(levels), function(i) {
      mean(actual >= pred_df[[lower_cols[i]]] &
           actual <= pred_df[[upper_cols[i]]])
    }),
    paste0("coverage_", levels)
  )

  y    <- actual
  yhat <- pred_df$pred

  c(
    list(
      wis  = sc$wis,
      mae  = mean(abs(y - yhat)),
      mse  = mean((y - yhat)^2),
      rmse = sqrt(mean((y - yhat)^2)),
      mape = mean(abs((y - yhat) / pmax(y, 1))) * 100
    ),
    as.list(coverage)
  )
}

# ── Main function ─────────────────────────────────────────────────────────────

#' Grid search over ARIMA orders and covariate formulas with CV evaluation
#'
#' @param data       Data frame with all covariates and train/target indicators.
#' @param formulas   List of formulas (from build_covariate_combinations output).
#' @param quantiles  Numeric vector of coverage levels, e.g. c(0.80, 0.95).
#' @param seasonal   Seasonal period list passed to fit_sarimax.
#' @param max_order  Named integer vector with upper bounds for each order
#'                   component: c(p=2, d=2, q=2, P=1, D=1, Q=1).
#' @param n_cores    Number of parallel workers for the order loop.
#'                   Defaults to all available cores minus 1.
#' @param table_path Optional: passed through to the per-fold
#'                   `forecast_covariates()` precompute step (see
#'                   `covariate_forecast_table_path()`) so the fold
#'                   covariate forecasts this grid search needs are looked
#'                   up in the persistent table first, and only fit live
#'                   (then saved) on a miss -- e.g. reusing forecasts a
#'                   prior `model_sel.r`, `fit.r`, or
#'                   `evaluate_covariate_forecasts.r` run already computed
#'                   for this exact unit/fold. Default `NULL` preserves the
#'                   original always-fit-live behavior exactly.
#' @param id         Required when `table_path` is supplied: this unit's id
#'                   (UF code or municipality geocode).
#'
#' @return A list with two tibbles:
#'   $predictions : one row per formula x order x train_id x date
#'   $metrics     : one row per formula x order x train_id, with all scores
run_grid_search <- function(data,
                            formulas,
                            levels = c(50, 80, 90, 95),
                            method = "CSS",
                            lambda = NULL,
                            optim.control = list(maxit = 500),   # more iterations
                            optim.method  = "BFGS",
                            seasonal  = list(order = c(1, 0, 1), period = 52),
                            max_order = c(p = 2, d = 2, q = 2,
                                          P = 1, D = 1, Q = 1),
                            fixed_order = F,
                            # FALSE searches d (and D) over the grid via CV WIS instead of
                            # locking it to a single KPSS/OCSB-recommended value: on short,
                            # irregular epidemic series the auto test can pick a worse d than
                            # what actually scores best out-of-sample.
                            fixed_stat_par = F,
                            bootstrap = TRUE,
                            npaths    = 1000,
                            n_cores   = max(1L, detectCores() - 1L),
                            table_path = NULL,
                            id         = NULL) {

  # ── 1. Build order grid ───────────────────────────────────────────────────
  # `if (fixed_order)` on a real named order vector (e.g. `c(p=1, d=1, ...)`,
  # as model_sel.r's pca = FALSE branch passes on its second run_grid_search()
  # call) is a length-6 condition -- R >= 4.3 makes that a hard error
  # ("the condition has length > 1"), where older R merely warned and used
  # the first element. isFALSE() correctly distinguishes "the default FALSE
  # sentinel" from "a real order vector was passed" regardless of R version.
  if (!isFALSE(fixed_order)) {
    # Single row grid from fixed order — no filtering needed
    order_grid <- data.frame(
      p = fixed_order["p"], d = fixed_order["d"], q = fixed_order["q"],
      P = fixed_order["P"], D = fixed_order["D"], Q = fixed_order["Q"],
      row.names = NULL
    )
  } else {
    required_names <- c("p", "d", "q", "P", "D", "Q")
    stopifnot(all(required_names %in% names(max_order)))
    stopifnot(all(max_order >= 0))

    if (fixed_stat_par) {
      diff_par <- determine_d(log1p(data$cases))
    }

    order_grid <- expand.grid(
      p = seq(0, max_order["p"]),
      d = if (fixed_stat_par) diff_par$d else seq(0, max_order["d"]),
      q = seq(0, max_order["q"]),
      P = seq(0, max_order["P"]),
      D = if (fixed_stat_par) diff_par$D else seq(0, max_order["D"]),
      Q = seq(0, max_order["Q"])
    ) |> filter(!(p == 0 & q == 0), !(P == 0 & Q == 0))
  }

  train_ids  <- paste0("train_", 1:4)
  n_orders   <- nrow(order_grid)
  n_formulas <- length(formulas)

  # ── 3. Build flat grid: every (formula_idx, order_idx) pair ───────────────
  job_grid <- expand.grid(
    fi = seq_len(n_formulas),
    oi = seq_len(n_orders)
  )
  total_jobs <- nrow(job_grid)

  message(sprintf(
    "Grid: %d formula(s) x %d order(s) x %d splits = %d total fits.",
    n_formulas, n_orders, length(train_ids),
    total_jobs * length(train_ids)
  ))

  # ── 4. Precompute formula metadata (outside workers) ──────────────────────
  formula_ids <- vapply(formulas, FUN.VALUE = character(1), function(f) {
    if (is.character(f)) {
      if (length(f) > 1) {
        paste(trimws(f), collapse = "+")
      } else { # length == 1
        single <- f[[1]]
        if (grepl("~", single)) {
          rhs <- sub(".*~", "", single)
          terms <- strsplit(rhs, "\\+")[[1]]
          paste(trimws(terms), collapse = "+")
        } else {
          paste(trimws(single), collapse = "+")
        }
      }
    } else if (inherits(f, "formula")) {
      paste(attr(stats::terms(f), "term.labels"), collapse = "+")
    } else {
      stop("Each element of 'formulas' must be a formula or character string/vector")
    }
  })
  # ── 4b. Precompute covariate forecasts once per fold, not once per grid
  #        cell. Every formula in `formulas` draws its regressors from the
  #        same small candidate universe (PCA components / surviving
  #        ocean-climate indices); forecasting each one once per fold here
  #        -- instead of once per (formula, order, fold) combination inside
  #        the parallel loop -- is what keeps genuinely forecasting
  #        covariates for every validation split computationally tractable.
  all_rhs_terms <- unique(unlist(lapply(formulas, parse_rhs_terms)))

  fold_windows <- setNames(lapply(train_ids, function(train_id) {
    target_id <- sub("train_", "target_", train_id)
    train_ew  <- data$epiweek[data[[train_id]]  == 1]
    target_ew <- data$epiweek[data[[target_id]] == 1]
    if (length(train_ew) == 0)  stop("No training rows found for ", train_id)
    if (length(target_ew) == 0) stop("No target rows found for ", target_id)
    list(
      train_start    = min(train_ew),
      train_end      = max(train_ew),
      forecast_start = min(target_ew),
      forecast_end   = max(target_ew)
    )
  }), train_ids)

  fold_covariate_forecasts <- if (length(all_rhs_terms) == 0) {
    setNames(vector("list", length(train_ids)), train_ids)
  } else {
    message(sprintf(
      "Precomputing covariate forecasts for %d candidate regressor(s) x %d fold(s)...",
      length(all_rhs_terms), length(train_ids)
    ))
    setNames(lapply(train_ids, function(train_id) {
      w <- fold_windows[[train_id]]
      forecast_covariates(
        data        = data,
        rhs_terms   = all_rhs_terms,
        train_start = w$train_start,
        train_end   = w$train_end,
        h           = length(enumerate_epiweeks(w$forecast_start, w$forecast_end)),
        table_path  = table_path,
        id          = id
      )
    }), train_ids)
  }

  # ── 5. Parallel backend ───────────────────────────────────────────────────
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)          # no progress_file to unlink

  clusterExport(cl, varlist = c(
    "fit_sarimax", "compute_metrics",
    "enumerate_epiweeks", "epiweek_to_date", "has_53_weeks", "parse_rhs_terms",
    "forecast_covariates", "forecast_covariate_series", "parse_derived_covariate",
    ".epiweek_to_date", ".has_53_weeks", ".next_epiweek", ".epiweek_seq_forward",
    "data", "formulas", "formula_ids",
    "order_grid", "job_grid",
    "fold_windows", "fold_covariate_forecasts",
    "seasonal", "levels", "train_ids",
    "levels", "method", "lambda", "optim.control", "optim.method",
    "bootstrap", "npaths"
  ), envir = environment())

  clusterEvalQ(cl, {
    library(forecast)
    library(dplyr)
    library(tidyr)
    library(scoringutils)
    library(pbapply)
  })

  message(sprintf("Running %d jobs on %d core(s)...", total_jobs, n_cores))
  start_time <- proc.time()["elapsed"]

  # ── 6. Single flat parallel loop over all (formula, order) pairs ──────────
  op <- pbapply::pboptions(type="timer") 
  pbapply::pboptions(op)
  all_results <- pbapply::pblapply(
    seq_len(total_jobs),
    function(job) {

      fi <- job_grid$fi[job]
      oi <- job_grid$oi[job]

      formula    <- formulas[[fi]]
      formula_id <- formula_ids[[fi]]
      og         <- order_grid[oi, ]
      ord        <- c(og$p, og$d, og$q)
      sea        <- list(order = c(og$P, og$D, og$Q), period = seasonal$period)
      order_str  <- sprintf("(%d,%d,%d)(%d,%d,%d)",
                            og$p, og$d, og$q, og$P, og$D, og$Q)

      split_results <- lapply(train_ids, function(train_id) {
        target_id <- sub("train_", "target_", train_id)
        actual    <- data$cases[data[[target_id]] == 1]
        w         <- fold_windows[[train_id]]

        # Capture WHY a fit failed, not just THAT it failed -- with every
        # failure previously collapsed to NULL, a run where every single fit
        # fails (as opposed to a handful of genuinely bad formula/order
        # combinations) gave no way to see the underlying Arima()/
        # compute_metrics() error short of re-running one fit manually
        # outside the cluster. fail_reason rides along in the fallback rows
        # below and gets summarized once, after the loop, exactly like
        # build_covariate_forecast_table.r's "Skip reasons:" summary.
        fail_reason <- NA_character_

        preds <- tryCatch(
          fit_sarimax(data,
                      formula        = formula,
                      train_start    = w$train_start,
                      train_end      = w$train_end,
                      forecast_start = w$forecast_start,
                      forecast_end   = w$forecast_end,
                      levels         = levels,
                      method         = method,
                      order          = ord,
                      seasonal       = sea,
                      lambda         = lambda,
                      optim.control  = optim.control,
                      optim.method   = optim.method,
                      bootstrap      = bootstrap,
                      npaths         = npaths,
                      covariate_forecasts = fold_covariate_forecasts[[train_id]]
                    ),
          error = function(e) {
            fail_reason <<- conditionMessage(e)
            NULL
          }
        )

        if (is.null(preds)) {
          # fail_reason already set by fit_sarimax()'s own tryCatch above.
        } else if (nrow(preds) != length(actual)) {
          fail_reason <- sprintf(
            "fit_sarimax() returned %d row(s), expected %d (target window mismatch)",
            nrow(preds), length(actual)
          )
        }
        failed  <- is.null(preds) || nrow(preds) != length(actual)
        metrics <- if (!failed) {
          tryCatch(
            as_tibble(compute_metrics(preds, actual)),
            error = function(e) {
              fail_reason <<- paste("compute_metrics() failed:", conditionMessage(e))
              NULL
            }
          )
        } else NULL
        failed <- failed || is.null(metrics)

        list(
          predictions = if (!failed) {
            preds |> mutate(formula_id = formula_id,
                            order      = order_str,
                            train_id   = train_id,
                            failed     = FALSE)
          } else {
            tibble(formula_id = formula_id, order = order_str,
                  train_id = train_id, failed = TRUE, fail_reason = fail_reason)
          },
          metrics = if (!failed) {
            metrics |> mutate(formula_id = formula_id,
                              order      = order_str,
                              train_id   = train_id,
                              failed     = FALSE)
          } else {
            tibble(formula_id = formula_id, order = order_str,
                  train_id = train_id, failed = TRUE, fail_reason = fail_reason)
          }
        )
      })

      split_results
    },
    cl = cl   # PSOCK cluster — pblapply handles distribution automatically
  )

  # ── 7. Progress summary ───────────────────────────────────────────────────
  elapsed_total <- proc.time()["elapsed"] - start_time
  fmt_time <- function(secs) {
    if (secs < 60)    return(sprintf("%.0fs", secs))
    if (secs < 3600)  return(sprintf("%.0fm %.0fs", secs %/% 60, secs %% 60))
    sprintf("%.0fh %.0fm", secs %/% 3600, (secs %% 3600) %/% 60)
  }

  # ── 8. Flatten and bind ───────────────────────────────────────────────────
  flat <- unlist(all_results, recursive = FALSE)

  predictions <- bind_rows(lapply(flat, `[[`, "predictions"))
  metrics     <- bind_rows(lapply(flat, `[[`, "metrics"))

  n_failed <- sum(metrics$failed, na.rm = TRUE)
  message(sprintf("Done. Total time: %s. Failed fits: %d / %d.",
                  fmt_time(elapsed_total), n_failed,
                  total_jobs * length(train_ids)))
  if (n_failed > 0 && "fail_reason" %in% names(metrics)) {
    # Printed here rather than left to each job's own message() so the
    # reasons are visible even though every job ran on a parallel worker
    # whose message() output is otherwise silently discarded (PSOCK workers'
    # outfile is /dev/null by default) -- see build_covariate_forecast_
    # table.r's identical "Skip reasons:" pattern.
    reason_counts <- sort(table(metrics$fail_reason[metrics$failed]), decreasing = TRUE)
    message("Failure reasons:")
    for (r in names(reason_counts)) {
      message(sprintf("  %d x  %s", reason_counts[[r]], r))
    }
  }

  list(predictions = predictions, metrics = metrics)
}

#' Reduce each lag/rolling family to its single most-predictive candidate
#'
#' For every base variable (after stripping "_lag<k>" / "_mean_<k>mo"
#' suffixes), keep only the one candidate column most correlated with
#' log1p(response).
#'
#' @param data        Data frame containing `covariates` and `response`.
#' @param covariates  Character vector of candidate covariate names, possibly
#'                    spanning several lag/rolling variants per base variable.
#' @param response    Name of the response column to correlate against
#'                    (default "cases"; correlated on the log1p scale).
#'
#' @return Character vector with one (best) covariate name per base variable.
select_best_per_variable <- function(data, covariates, response = "cases") {
  # Extract base variable name by stripping lag/rolling suffixes
  base_var <- gsub("_lag\\d+$|_mean_\\d+mo$", "", covariates)

  y <- log1p(data[[response]])

  # For each base variable, keep the candidate most correlated with response
  result <- tapply(covariates, base_var, function(candidates) {
    cors <- sapply(candidates, function(v)
      abs(cor(data[[v]], y, use = "pairwise.complete.obs")))
    candidates[which.max(cors)]
  })

  unname(unlist(result))
}

#' Drop covariates weakly correlated with the (log1p) response
#'
#' @param data        Data frame containing `covariates` and `response`.
#' @param covariates  Character vector of candidate covariate names.
#' @param response    Name of the response column (default "cases";
#'                    correlated on the log1p scale).
#' @param min_cor     Minimum absolute Pearson correlation required to keep a
#'                    covariate (default 0.1).
#'
#' @return Character vector of covariate names passing the correlation filter.
filter_by_correlation <- function(data, covariates,
                                  response  = "cases",
                                  min_cor   = 0.1) {
  y    <- log1p(data[[response]])
  cors <- sapply(covariates, function(v)
    abs(cor(data[[v]], y, use = "pairwise.complete.obs")))

  keep    <- names(cors[cors >= min_cor])
  dropped <- setdiff(covariates, keep)

  if (length(dropped) > 0)
    message(sprintf("Dropped %d weak predictor(s): %s",
                    length(dropped), paste(dropped, collapse = ", ")))
  keep
}

#' Drop climate indices (ENSO/IOD/PDO) redundant with local weather covariates
#'
#' Residualises log1p(cases) on the local-weather covariates, then keeps an
#' index only if its correlation with that residual still exceeds `threshold`
#' — i.e. the index explains variation not already captured by local weather.
#'
#' @param data        Data frame containing `covariates` and `cases`.
#' @param covariates  Character vector of candidate covariate names (mix of
#'                    local weather variables and climate indices).
#' @param indices     Names treated as climate indices (default
#'                    c("enso", "iod", "pdo")).
#' @param threshold   Minimum absolute residual correlation required to keep
#'                    an index (default 0.3).
#'
#' @return Character vector: all local-weather covariates plus any indices
#'         that passed the redundancy filter.
filter_redundant_indices <- function(data, covariates,
                                     indices   = c("enso", "iod", "pdo"),
                                     threshold = 0.3) {
  # Keep an index only if its partial correlation with response
  # (after removing linear effect of local weather) exceeds threshold
  local_weather <- setdiff(covariates, indices)
  present_indices <- intersect(indices, covariates)

  if (length(local_weather) == 0 || length(present_indices) == 0)
    return(covariates)

  y       <- log1p(data$cases)
  X_local <- as.matrix(data[, local_weather, drop = FALSE])

  # Residualise response on local weather
  y_resid <- residuals(lm(y ~ X_local))

  keep_indices <- Filter(function(idx) {
    abs(cor(data[[idx]], y_resid, use = "pairwise.complete.obs")) >= threshold
  }, present_indices)

  dropped <- setdiff(present_indices, keep_indices)
  if (length(dropped) > 0)
    message(sprintf("Dropped redundant index/indices: %s",
                    paste(dropped, collapse = ", ")))

  c(local_weather, keep_indices)
}

#' Parse a "(p,d,q)(P,D,Q)"-style order string back into numeric vectors
#'
#' @param order_str  Character string formatted like "(1,1,1)(1,0,1)", as
#'                   produced by `run_grid_search()`.
#'
#' @return A list with `order` (p,d,q) and `seasonal_order` (P,D,Q) integer
#'         vectors.
parse_order <- function(order_str) {
  nums <- as.integer(regmatches(order_str, gregexpr("[0-9]", order_str))[[1]])
  list(order = nums[1:3], seasonal_order = nums[4:6])
}

#' Recommend non-seasonal and seasonal differencing orders for a series
#'
#' @param y      Numeric vector or ts object (typically log1p(cases)).
#' @param max_d  Maximum non-seasonal differencing order to consider for the
#'               KPSS test (default 2).
#'
#' @return A list with recommended `d` (KPSS test, non-seasonal) and `D`
#'         (OCSB test, seasonal period 52) differencing orders.
determine_d <- function(y, max_d = 2) {
  d  <- ndiffs(y,  test = "kpss", max.d = max_d)
  D  <- nsdiffs(y, test = "ocsb", m = 52)
  message(sprintf("Recommended: d = %d, D = %d", d, D))
  list(d = d, D = D)
}

#' Fold-consistent PCA dimensionality reduction over candidate covariates
#'
#' Fits PCA separately within each fold\'s training rows (to avoid leakage),
#' chooses the number of components needed to explain `var_threshold` of
#' variance in the worst-case fold, then projects every fold\'s train+target
#' rows onto that common number of components.
#'
#' @param data           Data frame containing `candidates` and the
#'                       `train_cols` indicator columns.
#' @param candidates     Character vector of covariate names to reduce via PCA.
#' @param train_cols     Character vector of training-indicator column names,
#'                       one per CV fold (default paste0("train_", 1:4)).
#' @param var_threshold  Minimum cumulative explained-variance fraction used
#'                       to pick the number of components (default 0.90).
#'
#' @return A list with `data` (original data plus PC1..PCk score columns),
#'         `variables` (the new PC column names), and `var_table` (per-
#'         component and cumulative variance explained).
pca_all <- function(data,
                    candidates,
                    train_cols    = paste0("train_", 1:4),
                    var_threshold = 0.90) {

  mat    <- as.matrix(data[, candidates, drop = FALSE])
  n_rows <- nrow(data)

  # ── 1. Pre-compute n_comp consistently across all folds ───────────────────
  n_comps <- sapply(seq_along(train_cols), function(i) {
    train_rows <- data[[train_cols[i]]] == 1
    pca_fit    <- prcomp(mat[train_rows, , drop = FALSE],
                         center = TRUE, scale. = TRUE)
    var_exp    <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
    max(1L, which(var_exp >= var_threshold)[1])
  })

  # Use minimum across folds so all folds produce the same number of components
  n_comp     <- min(n_comps)
  comp_names <- paste0("PC", seq_len(n_comp))
  scores_out <- matrix(NA_real_, nrow = n_rows, ncol = n_comp,
                       dimnames = list(NULL, comp_names))

  message(sprintf("Using %d components (min across folds: %s).",
                  n_comp, paste(n_comps, collapse = ", ")))

  # ── 2. Fit per fold and write scores ──────────────────────────────────────
  pca_fit_last <- NULL

  for (i in seq_along(train_cols)) {
    train_col  <- train_cols[i]
    target_col <- sub("train_", "target_", train_col)
    train_rows <- data[[train_col]] == 1
    fold_rows  <- data[[train_col]] == 1 | data[[target_col]] == 1

    message(sprintf("Fold %d: fitting PCA on %d training rows.", i, sum(train_rows)))

    pca_fit <- prcomp(mat[train_rows, , drop = FALSE],
                      center = TRUE, scale. = TRUE)

    mat_scaled <- scale(mat[fold_rows, , drop = FALSE],
                        center = pca_fit$center,
                        scale  = pca_fit$scale)

    # Always take exactly n_comp components
    scores <- mat_scaled %*% pca_fit$rotation[, seq_len(n_comp), drop = FALSE]

    scores_out[fold_rows, ] <- scores
    pca_fit_last <- pca_fit
  }

  # ── 3. Variance table (based on consistent n_comp) ────────────────────────
  var_table <- tibble(
    component = comp_names,
    var_exp   = pca_fit_last$sdev[seq_len(n_comp)]^2 / sum(pca_fit_last$sdev^2),
    cum_var   = cumsum(pca_fit_last$sdev[seq_len(n_comp)]^2 /
                       sum(pca_fit_last$sdev^2))
  )
  print(var_table)

  list(
    data      = bind_cols(data, as_tibble(scores_out)),
    variables = comp_names,
    var_table = var_table
  )
}