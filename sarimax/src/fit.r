source("sarimax/src/utils.r")  # also sources forecast_covariates.r
source("sarimax/src/validate_challenge_rules.r")

# ── Submission window config ────────────────────────────────────────────────
# Change ONLY this block to switch which window the train_4/target_4 splits
# below fit and forecast (they read SUBMISSION_WINDOW, not literal numbers).
# validate_epiweek_window() aborts immediately -- before any of the long
# per-unit fitting loops run -- if the resulting bounds don't match what the
# IMDC rules allow for the chosen phase; see sarimax/src/validate_challenge_rules.r
# and root README.md, Section 6 ("Data Usage Restriction").
#
# Currently set to reproduce Validation Test 4 (2025-2026 season, due July 1).
# For the Forecast Phase (2026-2027 season, due September 10, after the July 31
# data update), set instead:
#   SUBMISSION_PHASE        <- "forecast"
#   SUBMISSION_TEST_ID       <- NULL
#   SUBMISSION_FORECAST_YEAR <- 2026
SUBMISSION_PHASE         <- "validation"
SUBMISSION_TEST_ID       <- 4
SUBMISSION_FORECAST_YEAR <- NULL

SUBMISSION_WINDOW <- expected_challenge_window(
  phase         = SUBMISSION_PHASE,
  test_id       = SUBMISSION_TEST_ID,
  forecast_year = SUBMISSION_FORECAST_YEAR
)
validate_epiweek_window(
  train_start    = SUBMISSION_WINDOW$train_start,
  train_end      = SUBMISSION_WINDOW$train_end,
  forecast_start = SUBMISSION_WINDOW$forecast_start,
  forecast_end   = SUBMISSION_WINDOW$forecast_end,
  phase          = SUBMISSION_PHASE,
  test_id        = SUBMISSION_TEST_ID,
  forecast_year  = SUBMISSION_FORECAST_YEAR
)
message(
  "Submission window: ", SUBMISSION_WINDOW$label, " -- train ",
  SUBMISSION_WINDOW$train_start, "..", SUBMISSION_WINDOW$train_end,
  ", forecast ", SUBMISSION_WINDOW$forecast_start, "..", SUBMISSION_WINDOW$forecast_end
)

#' Resolve the train/forecast epiweek window for one of the four validation
#' splits used throughout this file. Split 4 uses SUBMISSION_WINDOW (see the
#' config block above, since that is the one split whose window changes
#' depending on which phase this script is currently producing); splits 1-3
#' are the three already-elapsed retrospective seasons, which always use the
#' canonical challenge window for their own year regardless of
#' SUBMISSION_WINDOW.
#'
#' @param train_id  One of "train_1", "train_2", "train_3", "train_4".
#' @param data      The unit/disease data frame (used only to find the
#'                  earliest epiweek actually present, e.g. chikungunya's
#'                  later start).
#'
#' @return A list with train_start, train_end, forecast_start, forecast_end.
window_for_train_id <- function(train_id, data) {
  test_id <- as.integer(sub("train_", "", train_id))
  if (test_id == 4) return(SUBMISSION_WINDOW)
  expected_challenge_window("validation", test_id = test_id, data_start = min(data$epiweek))
}

# ── Candidate covariate selection + PCA, shared by every unit/level/disease ─
#
# Identical steps `run_model_selection()` (sarimax/src/model_sel.r) uses to
# build the covariate space a formula's PCs are drawn from -- kept in sync
# here so the formula/order picked by model selection resolves to the same
# regressor columns at final-fit time.
prepare_unit_covariates <- function(data) {
  train_rows <- data$train_1 == 1

  candidates <- get_candidates(data)
  candidates <- filter_low_variance(data, candidates, threshold = 0.01)
  candidates <- filter_by_correlation(data[train_rows, ], candidates, min_cor = 0.1)

  pca_result <- pca_all(
    data          = data,
    candidates    = candidates[!grepl("enso|iod|pdo", candidates)],
    var_threshold = 0.9
  )
  pca_result$data
}

# ── Fit + forecast every validation split for one unit (state or city) ─────
#
# Loads the unit's data, rebuilds its PCA covariate space, then fits the
# challenge's four train/target splits with the formula/order
# `run_model_selection()` selected as best for this unit. If that fit
# produces an actionable warning (a NaN standard error, an internal
# `auto.arima()` fallback, or an outright fit failure) on any split, it is
# retried with the next-best formula/order from this unit's full
# leaderboard (`metrics_all_formulas_<disease>_<id>.csv`), and so on, until
# a clean fit is found or the leaderboard is exhausted. This generalizes
# what used to be a chikungunya-only retry loop to every disease and level,
# since a warning is just as much a reason to prefer a different candidate
# for dengue (or a municipality) as it is for chikungunya.
#
# @param disease      "dengue" or "chikungunya".
# @param level        "state" or "city".
# @param id           UF code or municipality geocode.
# @param best_wis_df  This level's leaderboard-winners table (one row per
#                      unit), as written by `run_model_selection()` to
#                      `best_wis_<disease>_all_<states|cities>.csv`.
# @param unit_col     Name of `best_wis_df`'s unit-identifier column
#                      ("state" or "city").
# @param metrics_dir  Directory holding each unit's full formula/order
#                      leaderboard (default "sarimax/results/metrics").
# @param preds_dir    Directory predictions are written to (default
#                      "sarimax/results/preds").
#
# Each split's `fit_sarimax()` call passes this unit's persistent
# covariate-forecast table (`covariate_forecast_table_path()`) so a root
# series' forecast already computed by `model_sel.r`'s grid search for this
# exact unit/fold is reused here instead of being re-fit from scratch, and
# vice versa on a later `model_sel.r` re-run.
#
# @return A list: `id`, `resolved` (logical), `used_row` (1 = the top pick
#   worked untouched), `warnings` (named list by target_id), and `preds` (the
#   four splits' forecasts, row-bound, with a `target_id` and unit_col
#   column attached).
fit_unit_predictions <- function(disease, level, id, best_wis_df, unit_col,
                                  metrics_dir = "sarimax/results/metrics",
                                  preds_dir   = "sarimax/results/preds") {

  message(sprintf("[%s | %s | %s] Loading and preparing covariates...", disease, level, id))
  d <- load_unit_data(disease = disease, level = level, id = id)
  d <- prepare_unit_covariates(d)

  unit_row <- best_wis_df[best_wis_df[[unit_col]] == id, ]
  candidate_rows <- tibble(formula_id = unit_row$formula_id, order = unit_row$order)

  metrics_file <- file.path(metrics_dir, paste0("metrics_all_formulas_", disease, "_", id, ".csv"))
  if (file.exists(metrics_file)) {
    metrics_df <- read_csv(metrics_file, show_col_types = FALSE)
    # Rows 2.. are fallback candidates, already ranked by CV score.
    if (nrow(metrics_df) > 1) {
      candidate_rows <- bind_rows(candidate_rows, metrics_df[-1, c("formula_id", "order")])
    }
  }

  train_ids  <- paste0("train_", 1:4)
  target_ids <- paste0("target_", 1:4)

  resolved       <- FALSE
  first_attempt  <- NULL   # kept as the fallback if nothing ever resolves
  final_preds    <- NULL
  final_warnings <- list()
  used_row       <- 1L

  for (row_i in seq_len(nrow(candidate_rows))) {
    best_formula <- candidate_rows$formula_id[row_i]
    best_order   <- candidate_rows$order[row_i]
    ord          <- parse_order(best_order)

    if (row_i > 1) {
      message(sprintf("[%s] Trying formula row %d of %d: %s",
                      id, row_i, nrow(candidate_rows), best_formula))
    }

    attempt_warnings <- list()
    attempt_failed   <- FALSE

    pred_target <- lapply(seq_along(train_ids), function(j) {
      target_id <- target_ids[j]
      win <- window_for_train_id(train_ids[j], d)
      fit <- tryCatch(
        fit_sarimax(
          data           = d,
          formula        = best_formula,
          train_start    = win$train_start,
          train_end      = win$train_end,
          forecast_start = win$forecast_start,
          forecast_end   = win$forecast_end,
          order          = c(ord$order[1], ord$order[2], ord$order[3]),
          seasonal       = list(order = c(ord$seasonal_order[1], ord$seasonal_order[2], ord$seasonal_order[3]), period = 52),
          optim.method   = "BFGS",
          # Reuse whatever root-covariate forecasts model_sel.r (or an
          # earlier fit.r/evaluate_covariate_forecasts.r run) already
          # computed for this exact unit/fold, instead of re-fitting
          # tbats/stlf/Arima from scratch -- see covariate_forecast_table_path().
          table_path     = covariate_forecast_table_path(disease, level),
          id             = id
        ),
        error = function(e) {
          attempt_failed <<- TRUE
          message(sprintf("[%s | %s] fit_sarimax() failed: %s", id, target_id, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(fit)) return(NULL)

      w <- attr(fit, "warnings")
      if (!is.null(w)) {
        attempt_warnings[[target_id]] <<- w
        message(sprintf("[%s | %s] %d warning(s):\n%s",
                        id, target_id, length(w), paste0("  - ", w, collapse = "\n")))
      }
      fit
    })
    names(pred_target) <- target_ids

    if (row_i == 1) first_attempt <- list(preds = pred_target, warnings = attempt_warnings)

    actionable <- attempt_failed || any(unlist(lapply(attempt_warnings, function(w) {
      any(grepl("auto.arima|unreliable|failed", w, ignore.case = TRUE))
    })))

    if (!actionable) {
      resolved       <- TRUE
      used_row       <- row_i
      final_preds    <- pred_target
      final_warnings <- attempt_warnings
      if (row_i > 1) message(sprintf("[%s] Resolved with formula row %d: %s", id, row_i, best_formula))
      break
    }
    message(sprintf("[%s] Formula row %d still has warnings; trying next.", id, row_i))
  }

  if (!resolved) {
    message(sprintf(
      "[%s] All %d formula(s) exhausted; could not resolve warnings. Keeping the top-pick's forecast.",
      id, nrow(candidate_rows)
    ))
    final_preds    <- first_attempt$preds
    final_warnings <- first_attempt$warnings
  }

  for (j in seq_along(target_ids)) {
    if (is.null(final_preds[[j]])) {
      message(sprintf("[%s | %s] No forecast produced -- skipping this target's CSV.", id, target_ids[j]))
      next
    }
    write_csv(
      final_preds[[j]],
      file.path(preds_dir, paste0("pred_", disease, "_", id, "_", target_ids[j], ".csv"))
    )
  }

  list(
    id          = id,
    resolved    = resolved,
    used_row    = used_row,
    used_formula = candidate_rows$formula_id[used_row],
    used_order   = candidate_rows$order[used_row],
    warnings    = final_warnings,
    preds       = bind_rows(final_preds, .id = "target_id") |> mutate(!!unit_col := id)
  )
}

# ── Run every unit in one (disease, level)'s leaderboard ────────────────────
#
# Skips gracefully (with a message) if `run_model_selection()` hasn't been
# run yet for this level -- e.g. today, only the state-level leaderboards
# exist, so the city-level calls below are no-ops until model_sel.r has been
# run for cities too.
run_all_units <- function(disease, level) {
  unit_col <- if (level == "state") "state" else "city"
  best_wis_path <- file.path(
    "sarimax/results/metrics",
    paste0("best_wis_", disease, "_all_", if (level == "state") "states" else "cities", ".csv")
  )
  if (!file.exists(best_wis_path)) {
    message(sprintf(
      "[%s | %s] %s not found -- skipping (run model_sel.r for this disease/level first).",
      disease, level, best_wis_path
    ))
    return(invisible(NULL))
  }

  best_wis_df <- read_csv(best_wis_path, show_col_types = FALSE)
  ids <- best_wis_df[[unit_col]]

  results <- lapply(ids, function(id) {
    fit_unit_predictions(
      disease     = disease,
      level       = level,
      id          = id,
      best_wis_df = best_wis_df,
      unit_col    = unit_col
    )
  })
  names(results) <- ids

  resolved_formulas <- lapply(results, function(r) {
    if (r$used_row > 1) {
      tibble(id = r$id, formula_id = r$used_formula, order = r$used_order, row = r$used_row) |>
        rename(!!unit_col := id)
    } else {
      NULL
    }
  }) |> bind_rows()

  if (nrow(resolved_formulas) > 0) {
    write_csv(
      resolved_formulas,
      file.path("sarimax/results/metrics", paste0("resolved_formulas_", disease, "_", level, ".csv"))
    )
  }

  unresolved <- vapply(results, function(r) !r$resolved, logical(1))
  if (any(unresolved)) {
    message(sprintf(
      "[%s | %s] %d unit(s) never fully resolved warnings; kept their top-pick forecast: %s",
      disease, level, sum(unresolved), paste(ids[unresolved], collapse = ", ")
    ))
  }

  bind_rows(lapply(results, `[[`, "preds"))
}

# ── Run every disease x level combination ───────────────────────────────────
preds_dengue_state      <- run_all_units("dengue",      "state")
preds_chikungunya_state <- run_all_units("chikungunya", "state")
preds_dengue_city       <- run_all_units("dengue",      "city")
preds_chikungunya_city  <- run_all_units("chikungunya", "city")
