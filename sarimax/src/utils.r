library(forecast)
library(dplyr)
library(tidyverse)
library(tidyr)
library(purrr)
library(lubridate)
library(scoringutils)
library(parallel)
library(pbapply)

# ── IMDC temporal split helpers ─────────────────────────────────────────────

#' Convert a YYYYWW epiweek integer to the Date of its opening Sunday (MMWR)
epiweek_to_date <- function(yw) {
  yr   <- yw %/% 100L
  wk   <- yw  %% 100L
  jan4 <- as.Date(paste0(yr, "-01-04"))
  dow_jan4  <- as.integer(format(jan4, "%w"))
  sunday_w1 <- jan4 - dow_jan4
  sunday_w1 + (wk - 1L) * 7L
}

#' Check whether a year contains epidemiological week 53.
has_53_weeks <- function(year) {
  epiweek_to_date(year * 100L + 53L) <
    epiweek_to_date((year + 1L) * 100L + 1L)
}

#' Enumerate all epiweeks (inclusive) between two YYYYWW integers.
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

#' Official retrospective validation splits for the 3rd IMDC.
imdc_validation_splits <- function() {
  tibble::tibble(
    split_id       = 1:4,
    train_id       = paste0("train_", 1:4),
    target_id      = paste0("target_", 1:4),
    cutoff         = c(202225L, 202325L, 202425L, 202525L),
    forecast_start = c(202241L, 202341L, 202441L, 202541L),
    forecast_end   = c(202340L, 202440L, 202540L, 202640L)
  )
}

#' Final 2026-2027 forecast split.
imdc_final_split <- function() {
  tibble::tibble(
    split_id       = 5L,
    train_id       = "final_train",
    target_id      = "final",
    cutoff         = 202625L,
    forecast_start = 202641L,
    forecast_end   = 202740L
  )
}

#' Circular distance between epidemiological week numbers.
#'
#' Week 53 is treated as adjacent to week 52 and week 1. This is used only
#' for the short seasonal smoothing window applied to future climatology.
epiweek_distance <- function(a, b) {
  a <- pmin(as.integer(a), 52L)
  b <- pmin(as.integer(b), 52L)
  d <- abs(a - b)
  pmin(d, 52L - d)
}

#' Smoothed seasonal climatology using only observations available at cutoff.
seasonal_climatology_value <- function(data,
                                       variable,
                                       target_epiweek,
                                       cutoff,
                                       climatology_years = 5,
                                       smooth_weeks = 2) {
  if (!variable %in% names(data)) return(NA_real_)

  x <- data[[variable]]
  ok <- data$epiweek <= cutoff & is.finite(x)
  hist <- data[ok, c("epiweek", variable), drop = FALSE]

  if (nrow(hist) == 0) return(NA_real_)

  target_week <- target_epiweek %% 100L
  hist_week   <- hist$epiweek %% 100L
  keep_week   <- epiweek_distance(hist_week, target_week) <= smooth_weeks
  cand        <- hist[keep_week, , drop = FALSE]

  if (nrow(cand) == 0) {
    return(mean(hist[[variable]], na.rm = TRUE))
  }

  if (is.finite(climatology_years) && climatology_years > 0) {
    cand_year <- cand$epiweek %/% 100L
    years_keep <- head(sort(unique(cand_year), decreasing = TRUE),
                       climatology_years)
    cand <- cand[cand_year %in% years_keep, , drop = FALSE]
  }

  mean(cand[[variable]], na.rm = TRUE)
}

#' Aggregate the official monthly climate forecast to one modeled unit.
prepare_unit_climate_forecast <- function(forecasting_climate,
                                          cutoff,
                                          unit_type,
                                          unit_id,
                                          geocode_map = NULL) {
  if (is.null(forecasting_climate) || nrow(forecasting_climate) == 0)
    return(tibble::tibble())

  fc <- forecasting_climate
  fc$reference_month <- as.Date(fc$reference_month)

  if (unit_type == "city") {
    fc <- fc[fc$geocode == as.integer(unit_id), , drop = FALSE]
  } else if (unit_type == "state") {
    if (is.null(geocode_map))
      stop("geocode_map is required for state-level climate forecasts.")

    codes <- geocode_map |>
      dplyr::filter(uf == unit_id) |>
      dplyr::distinct(geocode) |>
      dplyr::pull(geocode)

    fc <- fc[fc$geocode %in% codes, , drop = FALSE]
  } else {
    stop("unit_type must be 'state' or 'city'.")
  }

  cutoff_date <- epiweek_to_date(cutoff)
  fc <- fc[
    fc$reference_month <= cutoff_date &
      fc$forecast_months_ahead >= 0 &
      fc$forecast_months_ahead <= 6,
    ,
    drop = FALSE
  ]
  if (nrow(fc) == 0) return(tibble::tibble())

  # Use the most recent vintage that already existed at the simulated cutoff.
  latest_reference <- max(fc$reference_month, na.rm = TRUE)
  fc <- fc[fc$reference_month == latest_reference, , drop = FALSE]

  fc <- fc |>
    dplyr::mutate(
      target_month = lubridate::floor_date(
        reference_month %m+% lubridate::months(forecast_months_ahead),
        unit = "month"
      )
    ) |>
    dplyr::group_by(target_month) |>
    dplyr::summarise(
      temp_med_mean = mean(temp_med, na.rm = TRUE),
      rel_humid_med_mean = mean(rel_humid_med, na.rm = TRUE),
      .groups = "drop"
    )

  attr(fc, "reference_month") <- latest_reference
  fc
}

#' Build the leakage-free "as of cutoff" data view used by the SARIMAX.
#'
#' Observed climate and cases are retained only through `cutoff`. From the
#' next epiweek onward, future climate is first filled with a smoothed seasonal
#' climatology computed from allowed historical data. Where the official
#' Copernicus forecast vintage available at the cutoff provides temperature
#' and relative humidity, those values overwrite the climatology. The official
#' forecast is monthly and limited to its available horizon; after it ends the
#' climatology remains in place. Lagged and rolling weather features are then
#' rebuilt over this permitted trajectory.
make_split_data <- function(data,
                            cutoff,
                            forecast_start,
                            forecast_end,
                            forecasting_climate = NULL,
                            unit_type = c("state", "city"),
                            unit_id = NULL,
                            geocode_map = NULL,
                            climatology_years = 5,
                            smooth_weeks = 2) {
  unit_type <- match.arg(unit_type)
  stopifnot(is.data.frame(data))
  stopifnot(all(c("epiweek", "cases") %in% names(data)))

  data <- data |>
    dplyr::arrange(epiweek) |>
    dplyr::distinct(epiweek, .keep_all = TRUE)

  observed <- data
  min_epiweek <- min(observed$epiweek, na.rm = TRUE)

  all_epiweeks <- enumerate_epiweeks(min_epiweek, forecast_end)
  trajectory <- tibble::tibble(
    epiweek = all_epiweeks,
    date    = epiweek_to_date(all_epiweeks),
    year    = all_epiweeks %/% 100L
  )

  join_data <- observed |>
    dplyr::select(-dplyr::any_of(c("date", "year")))

  trajectory <- trajectory |>
    dplyr::left_join(join_data, by = "epiweek")

  # Keep actual target observations outside the model view for scoring only.
  forecast_epiweeks <- enumerate_epiweeks(forecast_start, forecast_end)
  actual <- tibble::tibble(epiweek = forecast_epiweeks) |>
    dplyr::left_join(
      observed |> dplyr::select(epiweek, cases),
      by = "epiweek"
    ) |>
    dplyr::pull(cases)

  weather_vars <- intersect(
    c(
      "temp_min_mean", "temp_med_mean", "temp_max_mean",
      "precip_min_mean", "precip_med_mean", "precip_max_mean",
      "pressure_min_mean", "pressure_med_mean", "pressure_max_mean",
      "rel_humid_min_mean", "rel_humid_med_mean", "rel_humid_max_mean",
      "thermal_range_mean", "rainy_days_mean"
    ),
    names(trajectory)
  )
  climate_index_vars <- intersect(c("enso", "iod", "pdo"), names(trajectory))
  future_idx <- trajectory$epiweek > cutoff

  # Remove all actually observed future values before constructing substitutes.
  if ("cases" %in% names(trajectory))
    trajectory$cases[future_idx] <- NA_real_

  # Weather after the cutoff is unknown. Start from a smoothed seasonal
  # climatology built only from years already available at the cutoff.
  future_rows <- which(future_idx)
  for (v in weather_vars) {
    trajectory[[v]][future_idx] <- NA_real_

    if (length(future_rows) > 0) {
      trajectory[[v]][future_rows] <- vapply(
        trajectory$epiweek[future_rows],
        function(ew) seasonal_climatology_value(
          data = observed,
          variable = v,
          target_epiweek = ew,
          cutoff = cutoff,
          climatology_years = climatology_years,
          smooth_weeks = smooth_weeks
        ),
        numeric(1)
      )
    }
  }

  # No future ENSO/IOD/PDO forecast is distributed with the challenge data.
  # Use a simple persistence assumption: the last index value known at the
  # cutoff is carried forward. This is causal and avoids observed future values.
  for (v in climate_index_vars) {
    trajectory[[v]][future_idx] <- NA_real_
    known <- observed[[v]][observed$epiweek <= cutoff &
                             is.finite(observed[[v]])]
    if (length(known) > 0)
      trajectory[[v]][future_idx] <- tail(known, 1)
  }

  # Use the official forecast only for variables that have direct equivalents
  # in the distributed forecast file. `precip_tot` is deliberately not renamed
  # to precip_med/precip_min/precip_max because those are different quantities.
  fc_unit <- prepare_unit_climate_forecast(
    forecasting_climate = forecasting_climate,
    cutoff = cutoff,
    unit_type = unit_type,
    unit_id = unit_id,
    geocode_map = geocode_map
  )

  selected_vintage <- as.Date(NA_character_)
  if (nrow(fc_unit) > 0) {
    selected_vintage <- attr(fc_unit, "reference_month")
    trajectory <- trajectory |>
      dplyr::mutate(forecast_month = lubridate::floor_date(date, "month")) |>
      dplyr::left_join(
        fc_unit |>
          dplyr::rename(
            temp_med_mean_forecast = temp_med_mean,
            rel_humid_med_mean_forecast = rel_humid_med_mean
          ),
        by = c("forecast_month" = "target_month")
      )

    use_fc <- trajectory$epiweek > cutoff

    if ("temp_med_mean" %in% names(trajectory)) {
      idx <- use_fc & is.finite(trajectory$temp_med_mean_forecast)
      trajectory$temp_med_mean[idx] <- trajectory$temp_med_mean_forecast[idx]
    }
    if ("rel_humid_med_mean" %in% names(trajectory)) {
      idx <- use_fc & is.finite(trajectory$rel_humid_med_mean_forecast)
      trajectory$rel_humid_med_mean[idx] <-
        trajectory$rel_humid_med_mean_forecast[idx]
    }

    trajectory <- trajectory |>
      dplyr::select(
        -forecast_month,
        -dplyr::any_of(c(
          "temp_med_mean_forecast",
          "rel_humid_med_mean_forecast"
        ))
      )
  }

  # Carry static/reference fields into rows created beyond the observed table.
  static_cols <- intersect(
    c("uf", "uf_code", "geocode", "koppen", "biome", "disease"),
    names(trajectory)
  )
  known_rows <- observed$epiweek <= cutoff
  for (v in static_cols) {
    vals <- observed[[v]][known_rows]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0)
      trajectory[[v]][is.na(trajectory[[v]])] <- vals[length(vals)]
  }

  # Population is not a candidate regressor in the current SARIMAX, but keep a
  # cutoff-safe value in future rows for completeness.
  if ("pop" %in% names(trajectory)) {
    known_pop <- observed$pop[observed$epiweek <= cutoff & is.finite(observed$pop)]
    if (length(known_pop) > 0)
      trajectory$pop[future_idx] <- tail(known_pop, 1)
  }

  # Rebuild all historical lag/rolling weather features over the allowed path.
  for (v in weather_vars) {
    trajectory[[paste0(v, "_lag4")]]  <- dplyr::lag(trajectory[[v]], 4)
    trajectory[[paste0(v, "_lag8")]]  <- dplyr::lag(trajectory[[v]], 8)
    trajectory[[paste0(v, "_lag12")]] <- dplyr::lag(trajectory[[v]], 12)
    trajectory[[paste0(v, "_lag16")]] <- dplyr::lag(trajectory[[v]], 16)

    trajectory[[paste0(v, "_mean_3mo")]] <-
      dplyr::lag(runner::mean_run(trajectory[[v]], k = 12, na_rm = TRUE))
    trajectory[[paste0(v, "_mean_6mo")]] <-
      dplyr::lag(runner::mean_run(trajectory[[v]], k = 24, na_rm = TRUE))
    trajectory[[paste0(v, "_mean_9mo")]] <-
      dplyr::lag(runner::mean_run(trajectory[[v]], k = 36, na_rm = TRUE))
    trajectory[[paste0(v, "_mean_12mo")]] <-
      dplyr::lag(runner::mean_run(trajectory[[v]], k = 48, na_rm = TRUE))
  }

  trajectory <- trajectory |>
    dplyr::mutate(
      is_train = epiweek <= cutoff & !is.na(cases),
      is_forecast = epiweek >= forecast_start & epiweek <= forecast_end
    )

  list(
    data = trajectory,
    actual = actual,
    cutoff = cutoff,
    forecast_start = forecast_start,
    forecast_end = forecast_end,
    forecast_vintage = selected_vintage
  )
}

#' Determine a fixed number of PCs using only the earliest training split.
determine_pca_components <- function(data,
                                     candidates,
                                     var_threshold = 0.90,
                                     max_k = 5,
                                     threshold_low_variance = 0.01) {
  train <- data[data$is_train, , drop = FALSE]
  candidates <- intersect(candidates, names(train))

  sds <- sapply(candidates, function(v) sd(train[[v]], na.rm = TRUE))
  keep <- candidates[is.finite(sds) & sds > threshold_low_variance]
  if (length(keep) == 0) return(0L)

  mat <- as.matrix(train[, keep, drop = FALSE])
  means <- colMeans(mat, na.rm = TRUE)
  for (j in seq_len(ncol(mat))) {
    mat[is.na(mat[, j]), j] <- means[j]
  }

  fit <- prcomp(mat, center = TRUE, scale. = TRUE)
  var_exp <- cumsum(fit$sdev^2) / sum(fit$sdev^2)
  n_comp <- max(1L, which(var_exp >= var_threshold)[1])
  min(as.integer(max_k), n_comp, ncol(mat))
}

#' Fit PCA on one split's training rows and project its permitted future path.
pca_transform_split <- function(split,
                                candidates,
                                n_comp,
                                threshold_low_variance = 0.01) {
  if (n_comp <= 0) return(split)

  d <- split$data
  train <- d[d$is_train, , drop = FALSE]
  candidates <- intersect(candidates, names(d))

  sds <- sapply(candidates, function(v) sd(train[[v]], na.rm = TRUE))
  keep <- candidates[is.finite(sds) & sds > threshold_low_variance]

  if (length(keep) < n_comp) {
    stop("Not enough non-constant covariates to construct ", n_comp,
         " principal components at cutoff ", split$cutoff, ".")
  }

  train_mat <- as.matrix(train[, keep, drop = FALSE])
  means <- colMeans(train_mat, na.rm = TRUE)
  for (j in seq_len(ncol(train_mat))) {
    train_mat[is.na(train_mat[, j]), j] <- means[j]
  }

  pca_fit <- prcomp(train_mat, center = TRUE, scale. = TRUE)

  all_mat <- as.matrix(d[, keep, drop = FALSE])
  for (j in seq_len(ncol(all_mat))) {
    all_mat[is.na(all_mat[, j]), j] <- means[j]
  }

  all_scaled <- scale(
    all_mat,
    center = pca_fit$center,
    scale = pca_fit$scale
  )

  scores <- all_scaled %*%
    pca_fit$rotation[, seq_len(n_comp), drop = FALSE]
  colnames(scores) <- paste0("PC", seq_len(n_comp))

  split$data <- dplyr::bind_cols(
    d,
    tibble::as_tibble(scores)
  )
  attr(split$data, "pca_fit") <- pca_fit
  split
}

#' Build the four leakage-free validation splits and the common model recipe.
#'
#' Candidate screening is performed using only the training rows available
#' inside each cutoff. Split 1 is used only to lock the common formula/PC
#' structure used for fair comparison across backtests. PCA loadings and
#' standardisation are re-fitted independently inside each cutoff.
prepare_imdc_model_data <- function(data,
                                    forecasting_climate,
                                    geocode_map,
                                    unit_type,
                                    unit_id,
                                    pca = TRUE,
                                    pca_var_threshold = 0.90,
                                    k = 5,
                                    threshold_low_variance = 0.01,
                                    threshold_cor = 0.6,
                                    min_cor = 0.1,
                                    max_size_covariates = 3,
                                    index_cor_threshold = 0.3,
                                    climatology_years = 5,
                                    smooth_weeks = 2) {
  specs <- imdc_validation_splits()

  splits <- lapply(seq_len(nrow(specs)), function(i) {
    sp <- specs[i, ]
    s <- make_split_data(
      data = data,
      cutoff = sp$cutoff,
      forecast_start = sp$forecast_start,
      forecast_end = sp$forecast_end,
      forecasting_climate = forecasting_climate,
      unit_type = unit_type,
      unit_id = unit_id,
      geocode_map = geocode_map,
      climatology_years = climatology_years,
      smooth_weeks = smooth_weeks
    )
    s$split_id  <- sp$split_id
    s$train_id  <- sp$train_id
    s$target_id <- sp$target_id
    s
  })

  schema_candidates <- get_candidates(splits[[1]]$data)

  # The screening operation is repeated independently inside every cutoff.
  # No target row is used to decide which raw weather variables enter that
  # split's PCA. Split 1 is used only to lock the common PC count/formula
  # structure, because it is the earliest information set and is safe for all
  # later backtests.
  split_candidates <- lapply(splits, function(split) {
    train <- split$data[split$data$is_train, , drop = FALSE]
    cand <- filter_low_variance(
      train, schema_candidates, threshold = threshold_low_variance
    )
    filter_by_correlation(
      train, cand, min_cor = min_cor
    )
  })

  pca_n_comp <- 0L
  kept_indices <- character(0)

  if (pca) {
    first_candidates <- split_candidates[[1]]
    first_pca_candidates <-
      first_candidates[!grepl("enso|iod|pdo", first_candidates)]
    first_indices <- intersect(
      c("enso", "iod", "pdo"),
      first_candidates
    )

    if (length(first_pca_candidates) == 0) {
      formulas <- list(reformulate(
        if (length(first_indices) > 0) first_indices else character(0),
        response = "cases"
      ))
      kept_indices <- first_indices
    } else {
      pca_n_comp <- determine_pca_components(
        splits[[1]]$data,
        candidates = first_pca_candidates,
        var_threshold = pca_var_threshold,
        max_k = k,
        threshold_low_variance = threshold_low_variance
      )

      splits <- lapply(seq_along(splits), function(i) {
        pca_candidates_i <-
          split_candidates[[i]][
            !grepl("enso|iod|pdo", split_candidates[[i]])
          ]

        pca_transform_split(
          splits[[i]],
          candidates = pca_candidates_i,
          n_comp = pca_n_comp,
          threshold_low_variance = threshold_low_variance
        )
      })

      pcs <- paste0("PC", seq_len(pca_n_comp))
      train1_pca <- splits[[1]]$data[
        splits[[1]]$data$is_train, , drop = FALSE
      ]

      if (length(first_indices) > 0) {
        retained <- filter_redundant_indices(
          data = train1_pca,
          covariates = c(pcs, first_indices),
          indices = first_indices,
          threshold = index_cor_threshold
        )
        kept_indices <- intersect(first_indices, retained)
      }

      formulas <- lapply(seq_len(pca_n_comp), function(i) {
        reformulate(pcs[1:i], response = "cases")
      })

      if (length(kept_indices) > 0) {
        formulas <- c(formulas, lapply(seq_len(pca_n_comp), function(i) {
          reformulate(c(pcs[1:i], kept_indices), response = "cases")
        }))
      }
    }
  } else {
    # Non-PCA formulas are locked from the earliest cutoff so the 2022
    # validation cannot be influenced by later data.
    train1 <- splits[[1]]$data[splits[[1]]$data$is_train, , drop = FALSE]
    candidates <- select_best_per_variable(
      train1,
      split_candidates[[1]]
    )

    if (length(candidates) == 0) {
      formulas <- list(reformulate(character(0), response = "cases"))
    } else {
      combos <- build_covariate_combinations(
        data = train1,
        covariates = candidates,
        threshold = threshold_cor,
        max_size = max_size_covariates
      )
      formulas <- lapply(combos, \(vars) reformulate(vars, response = "cases"))
    }
  }

  list(
    splits = splits,
    formulas = formulas,
    schema_candidates = schema_candidates,
    split_candidates = split_candidates,
    pca = pca,
    pca_n_comp = pca_n_comp,
    kept_indices = kept_indices,
    threshold_low_variance = threshold_low_variance,
    threshold_cor = threshold_cor,
    min_cor = min_cor,
    max_size_covariates = max_size_covariates,
    index_cor_threshold = index_cor_threshold,
    climatology_years = climatology_years,
    smooth_weeks = smooth_weeks
  )
}

#' Prepare the final 2026-2027 split using the same recipe selected in CV.
prepare_imdc_final_data <- function(data,
                                    recipe,
                                    forecasting_climate,
                                    geocode_map,
                                    unit_type,
                                    unit_id) {
  sp <- imdc_final_split()
  split <- make_split_data(
    data = data,
    cutoff = sp$cutoff,
    forecast_start = sp$forecast_start,
    forecast_end = sp$forecast_end,
    forecasting_climate = forecasting_climate,
    unit_type = unit_type,
    unit_id = unit_id,
    geocode_map = geocode_map,
    climatology_years = recipe$climatology_years,
    smooth_weeks = recipe$smooth_weeks
  )
  split$split_id  <- sp$split_id
  split$train_id  <- sp$train_id
  split$target_id <- sp$target_id

  if (isTRUE(recipe$pca) && recipe$pca_n_comp > 0) {
    train <- split$data[split$data$is_train, , drop = FALSE]
    candidates <- get_candidates(split$data)
    candidates <- filter_low_variance(
      train,
      candidates,
      threshold = recipe$threshold_low_variance
    )
    candidates <- filter_by_correlation(
      train,
      candidates,
      min_cor = recipe$min_cor
    )
    pca_candidates <- candidates[
      !grepl("enso|iod|pdo", candidates)
    ]

    split <- pca_transform_split(
      split,
      candidates = pca_candidates,
      n_comp = recipe$pca_n_comp,
      threshold_low_variance = recipe$threshold_low_variance
    )
  }
  split
}

#' Fit a SARIMAX model and return forecasts with prediction intervals
#'
#' @param data        A data frame containing all columns referenced in formula,
#'                    plus `epiweek`, `cases`, and train/target indicator columns.
#' @param formula     A formula with `cases` as LHS and covariates as RHS,
#'                    e.g. ~ temp_med_mean + precip_med_mean + enso
#' @param train_id    Deprecated. Non-NULL values are rejected to prevent the
#'                    legacy target-covariate leakage path.
#' @param levels      Prediction interval coverages (50, 80, 90, 95).
#' @param order       ARIMA (p,d,q) order. Default c(1,1,1).
#' @param seasonal    Seasonal ARIMA list, e.g. list(order=c(1,0,1), period=52).
#'                    Default list(order=c(1,0,1), period=52) for weekly epiweek data.
#'
#' @return A tibble with columns: date, pred, lower_*, upper_*
fit_sarimax <- function(data,
                        formula,
                        train_id,
                        method = "CSS-ML",
                        lambda = NULL,
                        optim.control = list(maxit = 500),   # more iterations
                        optim.method  = "BFGS",
                        levels = c(50, 80, 90, 95),
                        order      = c(1, 1, 1),
                        seasonal   = list(order = c(1, 0, 1), period = 52),
                        bootstrap  = TRUE,   # simulate forecast paths instead of
                        npaths     = 1000,   # assuming Gaussian normal-theory intervals
                        seed       = 123) { 

  stopifnot(is.data.frame(data))
  stopifnot(all(c("is_train", "is_forecast", "cases", "date") %in% names(data)))

  if (!is.null(train_id)) {
    stop(
      "The legacy train_id/target_id path is disabled because it can expose ",
      "observed target covariates. Build the split with make_split_data() first."
    )
  }

  # ── 1. Split train / target rows ─────────────────────────────────────────
  train_rows  <- data[data$is_train, , drop = FALSE]
  target_rows <- data[data$is_forecast, , drop = FALSE]

  if (nrow(train_rows) == 0) stop("No training rows in prepared split.")
  if (nrow(target_rows) == 0) stop("No forecast rows in prepared split.")

  # ── 2. Log-transform response (log1p to handle zero counts) ──────────────
  if (is.null(lambda)) {
    y <- log1p(train_rows$cases)
  } else {
    y <- train_rows$cases
  }
  

  # ── 3. Build & standardize regressor matrices ─────────────────────────────
  rhs_terms <- {
    if (is.null(formula)) {
      character(0)
    } else if (is.character(formula)) {
      if (length(formula) == 0 || identical(formula, "")) character(0)
      else {
        fo <- stats::as.formula(paste("~", paste(formula, collapse = "+")))
        attr(stats::terms(fo), "term.labels")
      }
    } else {
      attr(stats::terms(formula), "term.labels")
    }
  }

  if (!all(rhs_terms %in% names(data))) {
    stop("Missing regressor(s): ",
         paste(setdiff(rhs_terms, names(data)), collapse = ", "))
  }

  # Drop only training rows that still lack a selected regressor. This mostly
  # affects the first weeks of lagged/rolling features. Future xreg must be
  # complete because it was explicitly constructed by make_split_data().
  train_keep <- is.finite(train_rows$cases)
  if (length(rhs_terms) > 0) {
    train_keep <- train_keep &
      stats::complete.cases(train_rows[, rhs_terms, drop = FALSE])
  }
  train_rows <- train_rows[train_keep, , drop = FALSE]

  if (nrow(train_rows) == 0)
    stop("No complete training rows remain for the selected formula.")

  if (length(rhs_terms) > 0 &&
      !all(stats::complete.cases(target_rows[, rhs_terms, drop = FALSE]))) {
    stop("Prepared future regressors contain NA values.")
  }

  y <- if (is.null(lambda)) log1p(train_rows$cases) else train_rows$cases

  if (length(rhs_terms) == 0) {
    xreg_train  <- NULL
    xreg_future <- NULL
  } else {
    raw_train  <- as.matrix(train_rows[,  rhs_terms, drop = FALSE])
    raw_future <- as.matrix(target_rows[, rhs_terms, drop = FALSE])

    # Compute mean and sd from training data only (no data leakage)
    col_means <- colMeans(raw_train, na.rm = TRUE)
    col_sds   <- apply(raw_train, 2, sd, na.rm = TRUE)

    # Avoid division by zero for constant columns
    col_sds[!is.finite(col_sds) | col_sds == 0] <- 1

    xreg_train  <- scale(raw_train,  center = col_means, scale = col_sds)
    xreg_future <- scale(raw_future, center = col_means, scale = col_sds)
  }

  # ── 4. Fit SARIMAX — capture warnings ─────────────────────────────────────
  y_ts <- ts(y, frequency = seasonal$period)
  fit_warnings <- character(0)

  fit <- withCallingHandlers(
    Arima(y_ts,
          order    = order,
          seasonal = seasonal,
          xreg     = xreg_train,
          method   = method,
          lambda   = lambda,
          optim.control = optim.control,
          optim.method  = optim.method
        ),
    error = function(e) stop("Arima() failed: ", conditionMessage(e)),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # ── 5. Detect NaN standard errors ─────────────────────────────────────────
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
      paste0(
        "NaN standard errors for: ",
        paste(nan_se_params, collapse = ", "),
        ". Prediction intervals may be unreliable. ",
        "Consider reducing model order or using auto.arima()."
      )
    )
  }

  # ── 6. Forecast & back-transform ──────────────────────────────────────────
  h      <- nrow(target_rows)
  set.seed(seed)

  sim_paths <- replicate(
    npaths,
    as.numeric(
      simulate(
        fit,
        nsim = h,
        xreg = xreg_future,
        future = TRUE,
        bootstrap = bootstrap
      )
    )
  )

  # ── 8. Assemble output tibble ─────────────────────────────────────────────
  sim_paths <- matrix(sim_paths, nrow = h, ncol = npaths)

  if (is.null(lambda)) {
    sim_paths <- expm1(sim_paths)
  }

  sim_paths <- pmax(sim_paths, 0)

  probs <- c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  qmat <- t(apply(
    sim_paths,
    1,
    stats::quantile,
    probs = probs,
    na.rm = TRUE,
    names = FALSE,
    type = 8
  ))
  

  for (lv in levels) {
    lv_char <- paste0(lv, "%")
    if (is.null(lambda)) {
      out[[paste0("lower_", lv)]] <- bt(as.numeric(fc$lower[, lv_char]))
      out[[paste0("upper_", lv)]] <- bt(as.numeric(fc$upper[, lv_char]))
    } else {
      out[[paste0("lower_", lv)]] <- as.numeric(fc$lower[, lv_char])
      out[[paste0("upper_", lv)]] <- as.numeric(fc$upper[, lv_char])
    }
  }


  out <- tibble::tibble(
    date     = as.Date(target_rows$date),
    lower_95 = qmat[, 1],
    lower_90 = qmat[, 2],
    lower_80 = qmat[, 3],
    lower_50 = qmat[, 4],
    pred     = qmat[, 5],
    upper_50 = qmat[, 6],
    upper_80 = qmat[, 7],
    upper_90 = qmat[, 8],
    upper_95 = qmat[, 9]
  ) |>
    dplyr::mutate(
      dplyr::across(
        c(pred, dplyr::starts_with("lower_"), dplyr::starts_with("upper_")),
        \(x) pmax(x, 0)
      )
    )
  
  # ── 9. Attach warnings as attribute ───────────────────────────────────────

  attr(out, "warnings") <- if (length(fit_warnings) > 0) fit_warnings else NULL
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
  if (length(covariates) == 0) return(character(0))
  vars_sd <- sapply(covariates, \(v) sd(data[[v]], na.rm = TRUE))
  keep    <- names(vars_sd[is.finite(vars_sd) & vars_sd > threshold])
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

  keep <- is.finite(actual)
  keep <- keep &
    stats::complete.cases(
      pred_df[, c("pred", lower_cols, upper_cols), drop = FALSE]
    )

  if (!any(keep))
    stop("No observed target rows are available for scoring.")

  pred_df <- pred_df[keep, , drop = FALSE]
  actual  <- actual[keep]
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
                            n_cores   = max(1L, detectCores() - 1L)) {

  stopifnot(is.list(splits), length(splits) > 0)

  # ── 1. Build order grid ───────────────────────────────────────────────────
  if (fixed_order) {
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
      first_train <- splits[[1]]$data |>
        dplyr::filter(is_train)
      diff_par <- determine_d(log1p(first_train$cases))
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
    n_formulas, n_orders, length(splits),
    total_jobs * length(splits)
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
  # ── 5. Parallel backend ───────────────────────────────────────────────────
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)          # no progress_file to unlink

  clusterExport(cl, varlist = c(
    "fit_sarimax", "compute_metrics",
    "splits", "formulas", "formula_ids",
    "order_grid", "job_grid",
    "seasonal", "levels", "method", "lambda",
    "optim.control", "optim.method",
    "bootstrap", "npaths", "seed"
  ), envir = environment())

  clusterEvalQ(cl, {
    library(forecast)
    library(dplyr)
    library(tidyr)
    library(scoringutils)
  })

  message(sprintf("Running %d jobs on %d core(s)...", total_jobs, n_cores))
  start_time <- proc.time()["elapsed"]

  # ── 6. Single flat parallel loop over all (formula, order) pairs ──────────
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

      split_results <- lapply(seq_along(splits), function(si) {
        split <- splits[[si]]

        preds <- tryCatch(
          fit_sarimax(
            split$data,
            formula = formula,
            levels = levels,
            method = method,
            order = ord,
            seasonal = sea,
            lambda = lambda,
            optim.control = optim.control,
            optim.method = optim.method,
            bootstrap = bootstrap,
            npaths = npaths,
            seed = seed + si
          ),
          error = function(e) NULL
        )

        failed <- is.null(preds) || nrow(preds) != length(split$actual)
        metrics <- if (!failed) {
          tryCatch(
            as_tibble(compute_metrics(preds, split$actual, levels = levels)),
            error = function(e) NULL
          )
        } else NULL

        list(
          predictions = if (!failed) {
            preds |>
              mutate(
                formula_id = formula_id,
                order = order_str,
                train_id = split$train_id,
                target_id = split$target_id,
                failed = FALSE
              )
          } else {
            tibble(
              formula_id = formula_id,
              order = order_str,
              train_id = split$train_id,
              target_id = split$target_id,
              failed = TRUE
            )
          },
          metrics = if (!is.null(metrics)) {
            metrics |>
              mutate(
                formula_id = formula_id,
                order = order_str,
                train_id = split$train_id,
                target_id = split$target_id,
                failed = FALSE
              )
          } else {
            tibble(
              formula_id = formula_id,
              order = order_str,
              train_id = split$train_id,
              target_id = split$target_id,
              failed = TRUE
            )
          }
        )
      })

      split_results
    },
    cl = cl   # PSOCK cluster — pblapply handles distribution automatically
  )

  # ── 7. Progress summary ───────────────────────────────────────────────────
  elapsed_total <- proc.time()["elapsed"] - start_time

  # ── 8. Flatten and bind ───────────────────────────────────────────────────
  flat <- unlist(all_results, recursive = FALSE)

  predictions <- bind_rows(lapply(flat, `[[`, "predictions"))
  metrics     <- bind_rows(lapply(flat, `[[`, "metrics"))

  n_failed <- sum(metrics$failed, na.rm = TRUE)
  message(sprintf(
    "Done. Total time: %.1f min. Failed fits: %d / %d.",
    elapsed_total / 60, n_failed,
    total_jobs * length(splits)
  ))

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
  if (length(covariates) == 0) return(character(0))
  # Extract base variable name by stripping lag/rolling suffixes
  base_var <- gsub("_lag\\d+$|_mean_\\d+mo$", "", covariates)

  y <- log1p(data[[response]])

  # For each base variable, keep the candidate most correlated with response
  result <- tapply(covariates, base_var, function(candidates) {
    cors <- sapply(candidates, function(v)
      abs(cor(data[[v]], y, use = "pairwise.complete.obs")))
    cors[!is.finite(cors)] <- -Inf
    if (all(cors == -Inf)) return(character(0))
    candidates[which.max(cors)]
  })
  
  unique(unname(unlist(result)))
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
  if (length(covariates) == 0) return(character(0))
  y    <- log1p(data[[response]])
  cors <- sapply(covariates, function(v)
    abs(cor(data[[v]], y, use = "pairwise.complete.obs")))

  keep    <- names(cors[is.finite(cors) & cors >= min_cor])
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

  needed <- c("cases", local_weather, present_indices)
  d <- data[stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
  if (nrow(d) < 3)
    return(covariates)

  y       <- log1p(data$cases)
  X_local <- as.matrix(data[, local_weather, drop = FALSE])

  # Residualise response on local weather
  y_resid <- residuals(lm(y ~ X_local))

  keep_indices <- Filter(function(idx) {
    r <- abs(cor(d[[idx]], y_resid))
    is.finite(r) && r >= threshold
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