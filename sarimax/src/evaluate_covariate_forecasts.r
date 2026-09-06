source("sarimax/src/utils.r")  # also sources forecast_covariates.r
source("sarimax/src/validate_challenge_rules.r")

# ── Why this file exists ─────────────────────────────────────────────────────
# fit_sarimax() (sarimax/src/utils.r) and run_grid_search() never see the real
# future values of a formula's covariates -- they call forecast_covariates()
# to produce genuine, model-based forecasts instead (see root README.md,
# Section 6, "Data Usage Restriction"). That keeps the pipeline ex-ante
# compliant, but it also means the disease forecast is only as good as those
# covariate forecasts are. This script answers the separate question: how
# good ARE they? For every validation split whose target window has already
# happened (splits 1-3 always; split 4 partially, once the season is
# underway), the real values are sitting right there in the processed data
# files -- so this re-runs forecast_covariates() exactly as fit.r would, and
# scores the result against what actually occurred.
#
# This is a diagnostic/reporting script, not part of the submission
# pipeline: nothing here feeds back into model_sel.r or fit.r.

# ── One unit, one metric table ──────────────────────────────────────────────

#' Score forecast_covariates() against real observed values for one unit
#'
#' For each requested split, re-forecasts `rhs_terms` over that split's
#' target window using only data through that split's own train_end (exactly
#' as `fit_sarimax()` does), then compares the forecast against whichever
#' target-window epiweeks are actually present in the data. A split with no
#' observed target-window data yet (most commonly split 4, before the season
#' it forecasts has happened) is skipped with a message, not an error.
#'
#' @param disease    "dengue" or "chikungunya".
#' @param level      "state" or "city".
#' @param id         UF code or municipality geocode.
#' @param rhs_terms  Character vector of covariates to check. Default: this
#'                    unit's contemporaneous candidates from `get_candidates()`
#'                    (the 5 local-weather aggregates plus enso/iod/pdo) --
#'                    the raw climate signal `pca_all()` builds PCs from.
#'                    Pass `groups = c("contemporaneous","lagged","rolling")`
#'                    candidates (state level only) to also check lagged/
#'                    rolling-window columns, which `forecast_covariates()`
#'                    forecasts no differently from any other series.
#' @param splits     Integer vector of which of the 4 validation splits to
#'                    check (default 1:4).
#' @param baseline   Logical: also score a naive seasonal-climatology
#'                    forecast (the training window's mean value for each
#'                    epiweek-of-year) alongside the model forecast, so
#'                    `skill_rmse` (1 - rmse_model/rmse_climatology; positive
#'                    means the model forecast beats naive seasonality) means
#'                    something (default TRUE).
#' @param table_path  Path to the persistent covariate-forecast table (see
#'                    `covariate_forecast_table_path()`), shared with
#'                    `model_sel.r` and `fit.r` -- default reuses whatever
#'                    they already computed for this unit/fold instead of
#'                    re-fitting. Pass `NULL` to force a fresh, uncached
#'                    evaluation (e.g. to double check the table isn't
#'                    masking a regression).
#'
#' @return A tibble, one row per (split, covariate): `n_weeks` compared,
#'   `mae`, `rmse`, `mape`, `bias`, `cor`, and if `baseline = TRUE`,
#'   `rmse_climatology` and `skill_rmse`. Empty (0-row) if no split had any
#'   observed data to check against.
evaluate_unit_covariate_forecasts <- function(disease, level, id,
                                               rhs_terms  = NULL,
                                               splits     = 1:4,
                                               baseline   = TRUE,
                                               table_path = covariate_forecast_table_path(disease, level)) {
  d <- load_unit_data(disease = disease, level = level, id = id)

  if (is.null(rhs_terms)) rhs_terms <- get_candidates(d, groups = "contemporaneous")
  missing_terms <- setdiff(rhs_terms, names(d))
  if (length(missing_terms) > 0) {
    stop("evaluate_unit_covariate_forecasts(): data is missing column(s): ",
         paste(missing_terms, collapse = ", "))
  }

  unit_col <- if (level == "state") "state" else "city"

  per_split <- lapply(splits, function(s) {
    win <- if (s == 4) {
      expected_challenge_window("validation", test_id = 4, data_start = min(d$epiweek))
    } else {
      expected_challenge_window("validation", test_id = s, data_start = min(d$epiweek))
    }

    fc_epiweeks <- enumerate_epiweeks(win$forecast_start, win$forecast_end)
    observed    <- d[d$epiweek %in% fc_epiweeks, ]

    if (nrow(observed) == 0) {
      message(sprintf("[%s | split %d] No observed data yet for %d..%d -- skipping.",
                      id, s, win$forecast_start, win$forecast_end))
      return(NULL)
    }

    h  <- length(fc_epiweeks)
    fc <- forecast_covariates(
      data        = d,
      rhs_terms   = rhs_terms,
      train_start = win$train_start,
      train_end   = win$train_end,
      h           = h,
      # Reuses whatever model_sel.r/fit.r already computed for this exact
      # unit/fold via the shared persistent table; pass table_path = NULL
      # explicitly to force a fresh, uncached evaluation instead.
      table_path  = table_path,
      id          = id
    )

    keep     <- fc_epiweeks %in% observed$epiweek
    fc       <- fc[keep, , drop = FALSE]
    observed <- observed[match(fc_epiweeks[keep], observed$epiweek), ]

    if (baseline) {
      train_rows  <- d[d$epiweek >= win$train_start & d$epiweek <= win$train_end, ]
      woy_train   <- train_rows$epiweek %% 100
      woy_target  <- observed$epiweek   %% 100
    }

    per_covariate <- lapply(rhs_terms, function(term) {
      pred <- fc[, term]
      obs  <- observed[[term]]
      err  <- pred - obs

      out <- tibble(
        split     = s,
        covariate = term,
        n_weeks   = length(obs),
        mae       = mean(abs(err)),
        rmse      = sqrt(mean(err^2)),
        mape      = mean(abs(err) / pmax(abs(obs), 1e-6)) * 100,
        bias      = mean(err),
        cor       = if (length(obs) > 1) suppressWarnings(stats::cor(pred, obs)) else NA_real_
      )

      if (baseline) {
        clim_by_woy <- tapply(train_rows[[term]], woy_train, mean, na.rm = TRUE)
        clim_pred   <- unname(clim_by_woy[as.character(woy_target)])
        clim_pred[is.na(clim_pred)] <- mean(train_rows[[term]], na.rm = TRUE)
        rmse_clim   <- sqrt(mean((clim_pred - obs)^2))
        out$rmse_climatology <- rmse_clim
        out$skill_rmse       <- if (rmse_clim > 0) 1 - out$rmse / rmse_clim else NA_real_
      }
      out
    })

    bind_rows(per_covariate)
  })

  out <- bind_rows(per_split)
  if (nrow(out) == 0) return(out)
  out |> mutate(!!unit_col := id, .before = 1)
}

# ── Every unit in one (disease, level) combination ──────────────────────────

#' Per-unit worker for evaluate_all_units_covariate_forecasts()
#'
#' Factored out to a top-level, named function purely so it can be shipped
#' to parallel workers by name via `clusterExport()` -- the same reason
#' `run_grid_search()` (`src/utils.r`) keeps its per-job logic addressable
#' rather than as an inline closure. Not meant to be called directly.
eval_one_unit_covariate_forecasts <- function(id, disease, level, ...) {
  tryCatch(
    evaluate_unit_covariate_forecasts(disease = disease, level = level, id = id, ...),
    error = function(e) {
      message(sprintf("[%s] evaluate_unit_covariate_forecasts() failed: %s", id, conditionMessage(e)))
      NULL
    }
  )
}

#' Run `evaluate_unit_covariate_forecasts()` over every unit of a level
#'
#' Unit list defaults to whichever units have a `best_wis_<disease>_all_
#' <states|cities>.csv` row (i.e. have been through `model_sel.r`); pass
#' `ids` explicitly to check units before/without running model selection.
#'
#' Each unit's evaluation is independent, CPU-bound work -- one
#' `tbats()`/`stlf()`/`Arima()` fit per covariate per split, via
#' `forecast_covariate_series()` -- so units are fanned out across a local
#' PSOCK cluster by default, exactly the way `run_grid_search()` parallelizes
#' across (formula, order) jobs. With the default 8 contemporaneous
#' covariates x 4 splits, that's 32 model fits per unit before any
#' parallelism, which is why a full run across e.g. 27 states is worth
#' spreading across cores rather than running one unit at a time.
#'
#' @param n_cores  Number of local worker processes (default
#'                  `detectCores() - 1`, capped at `length(ids)` -- no point
#'                  starting more workers than units to evaluate). Pass
#'                  `n_cores = 1` to force a plain sequential loop instead
#'                  (no cluster spun up at all) -- useful for debugging a
#'                  single failing unit, or on a machine where cluster
#'                  start-up overhead isn't worth it for just one or two
#'                  units.
#'
#' @return The row-bound per-unit results (possibly empty). Also written to
#'   `sarimax/results/metrics/covariate_forecast_eval_<disease>_<level>.csv`
#'   when non-empty.
evaluate_all_units_covariate_forecasts <- function(disease, level, ids = NULL,
                                                    n_cores = max(1L, detectCores() - 1L),
                                                    ...) {
  unit_col <- if (level == "state") "state" else "city"

  if (is.null(ids)) {
    best_wis_path <- file.path(
      "sarimax/results/metrics",
      paste0("best_wis_", disease, "_all_", if (level == "state") "states" else "cities", ".csv")
    )
    if (!file.exists(best_wis_path)) {
      message(sprintf(
        "[%s | %s] %s not found and no `ids` given -- skipping.",
        disease, level, best_wis_path
      ))
      return(invisible(tibble()))
    }
    ids <- read_csv(best_wis_path, show_col_types = FALSE)[[unit_col]]
  }

  n_cores <- max(1L, min(n_cores, length(ids)))

  # Sequential (n_cores == 1, e.g. explicitly requested, or only one unit to
  # evaluate) skips cluster start-up entirely; pblapply(cl = NULL) just runs
  # a plain, progress-barred lapply() in that case.
  cl <- NULL
  if (n_cores > 1) {
    message(sprintf("[%s | %s] Evaluating %d unit(s) on %d core(s)...",
                     disease, level, length(ids), n_cores))

    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)

    clusterExport(cl, varlist = c(
      "eval_one_unit_covariate_forecasts", "evaluate_unit_covariate_forecasts",
      "load_unit_data", "get_candidates",
      "forecast_covariates", "forecast_covariate_series", "parse_derived_covariate",
      "enumerate_epiweeks", "expected_challenge_window",
      "epiweek_to_date", "has_53_weeks",
      ".epiweek_to_date", ".has_53_weeks", ".next_epiweek", ".epiweek_seq_forward",
      "covariate_forecast_table_path", "get_or_compute_root_forecast",
      ".read_covariate_forecast_table", ".write_covariate_forecast_table", ".with_covariate_forecast_table_lock",
      ".covariate_forecast_table_cols"
    ), envir = environment())

    clusterEvalQ(cl, {
      library(forecast)
      library(dplyr)
      library(readr)
      library(tibble)
      library(tidyr)
    })

    op <- pbapply::pboptions(type = "timer")
    pbapply::pboptions(op)
  } else {
    message(sprintf("[%s | %s] Evaluating %d unit(s) sequentially (n_cores = 1)...",
                     disease, level, length(ids)))
  }

  results <- pbapply::pblapply(ids, eval_one_unit_covariate_forecasts,
                                disease = disease, level = level, ..., cl = cl)

  out <- bind_rows(results)
  if (nrow(out) > 0) {
    write_csv(
      out,
      file.path("sarimax/results/metrics", paste0("covariate_forecast_eval_", disease, "_", level, ".csv"))
    )
  }
  out
}

# ── Human-readable summary ──────────────────────────────────────────────────

#' Aggregate a covariate-forecast evaluation table into one row per covariate
#'
#' Averages across every unit and split scored, so you can see at a glance
#' which covariates `forecast_covariates()` forecasts well versus which ones
#' it's barely beating (or losing to) naive seasonal climatology.
#'
#' @param eval_df  Output of `evaluate_unit_covariate_forecasts()` or
#'                 `evaluate_all_units_covariate_forecasts()`.
#'
#' @return A tibble, one row per covariate, ordered by mean `skill_rmse`
#'   descending (if present) or mean `rmse` ascending otherwise.
summarize_covariate_forecast_eval <- function(eval_df) {
  if (nrow(eval_df) == 0) {
    message("summarize_covariate_forecast_eval(): nothing to summarize (empty input).")
    return(eval_df)
  }

  has_skill <- "skill_rmse" %in% names(eval_df)

  summary_df <- eval_df |>
    group_by(covariate) |>
    summarise(
      n_unit_splits = n(),
      mean_mae      = mean(mae, na.rm = TRUE),
      mean_rmse     = mean(rmse, na.rm = TRUE),
      mean_mape     = mean(mape, na.rm = TRUE),
      mean_bias     = mean(bias, na.rm = TRUE),
      mean_cor      = mean(cor, na.rm = TRUE),
      mean_skill_rmse = if (has_skill) mean(skill_rmse, na.rm = TRUE) else NA_real_,
      pct_beats_climatology = if (has_skill) mean(skill_rmse > 0, na.rm = TRUE) * 100 else NA_real_,
      .groups = "drop"
    )

  if (has_skill) {
    summary_df |> arrange(desc(mean_skill_rmse))
  } else {
    summary_df |> arrange(mean_rmse)
  }
}

# ── Run for every disease x level, if that level has been model-selected ───
covariate_eval_dengue_state      <- evaluate_all_units_covariate_forecasts("dengue",      "state")
covariate_eval_chikungunya_state <- evaluate_all_units_covariate_forecasts("chikungunya", "state")
covariate_eval_dengue_city       <- evaluate_all_units_covariate_forecasts("dengue",      "city")
covariate_eval_chikungunya_city  <- evaluate_all_units_covariate_forecasts("chikungunya", "city")

for (nm in c("covariate_eval_dengue_state", "covariate_eval_chikungunya_state",
             "covariate_eval_dengue_city",  "covariate_eval_chikungunya_city")) {
  df <- get(nm)
  if (nrow(df) > 0) {
    message("\n== ", nm, " summary (mean across all units/splits) ==")
    print(summarize_covariate_forecast_eval(df), n = Inf)
  }
}
