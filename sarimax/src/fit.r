source("sarimax/src/utils.r")  # also sources forecast_covariates.r
source("sarimax/src/validate_challenge_rules.r")

# ── Forecast Phase config ────────────────────────────────────────────────────
# The four validation splits (train_1/target_1 .. train_4/target_4) are
# always fixed by the challenge rules and never need any configuration here
# -- Validation Test 4 is always the 2025-2026 season (train through EW25
# 2025, forecast EW41 2025..EW40 2026), computed the same way as tests 1-3
# (see window_for_split() below). It is NOT affected by anything in this
# block.
#
# The Forecast Phase (challenge rules, Section 5.2) is a genuinely separate,
# later window: the 2026-2027 season, train through EW25 2026, forecast
# EW41 2026..EW40 2027. It is a DIFFERENT window from Validation Test 4, not
# the same slot under another name -- fitting/submitting it never touches
# target_4, and vice versa; each gets its own file
# (pred_<disease>_<id>_target_forecast.csv vs. ..._target_4.csv). Update
# SUBMISSION_FORECAST_YEAR each year when a new Forecast Phase opens.
#
# validate_epiweek_window() aborts immediately -- before any of the long
# per-unit fitting loops run -- if FORECAST_WINDOW doesn't match what the
# IMDC rules allow; see sarimax/src/validate_challenge_rules.r and root
# README.md, Section 6 ("Data Usage Restriction").
SUBMISSION_FORECAST_YEAR <- 2026

FORECAST_WINDOW <- expected_challenge_window(
  phase         = "forecast",
  forecast_year = SUBMISSION_FORECAST_YEAR
)
validate_epiweek_window(
  train_start    = FORECAST_WINDOW$train_start,
  train_end      = FORECAST_WINDOW$train_end,
  forecast_start = FORECAST_WINDOW$forecast_start,
  forecast_end   = FORECAST_WINDOW$forecast_end,
  phase          = "forecast",
  forecast_year  = SUBMISSION_FORECAST_YEAR
)
message(
  "Forecast Phase window: ", FORECAST_WINDOW$label, " -- train ",
  FORECAST_WINDOW$train_start, "..", FORECAST_WINDOW$train_end,
  ", forecast ", FORECAST_WINDOW$forecast_start, "..", FORECAST_WINDOW$forecast_end
)

#' Check whether processed_data actually covers FORECAST_WINDOW's train_end
#' yet, against one representative unit.
#'
#' Unlike the four validation splits (already-elapsed seasons, guaranteed
#' complete), FORECAST_WINDOW's train_end (EW25 of SUBMISSION_FORECAST_YEAR)
#' may be beyond the most recent real data actually on disk -- fitting on
#' data that stops short of that cutoff wouldn't error, it would silently
#' train on fewer weeks than the window claims. Called once, right before
#' the Forecast Phase driver calls at the bottom of this file (never before
#' the validation driver calls above, which don't depend on this at all),
#' against dengue/state/SP -- always present once model_sel.r has run at
#' the official submission level. Non-fatal if the check itself can't run
#' (e.g. that file doesn't exist yet in a fresh checkout) -- just a warning
#' to verify manually in that case, so a missing reference file never blocks
#' the fits that don't need it.
check_forecast_phase_data_ready <- function() {
  tryCatch({
    d_check   <- load_unit_data(disease = "dengue", level = "state", id = "SP")
    readiness <- validate_forecast_phase_data(d_check, forecast_year = SUBMISSION_FORECAST_YEAR)
    if (!readiness$data_ready) {
      stop(sprintf(
        "Forecast Phase data not ready: max epiweek available is %d, need through %d (EW25 %d). ",
        readiness$max_epiweek_available, readiness$needed_through_epiweek, SUBMISSION_FORECAST_YEAR
      ), "Update processed_data (data_prep/agg_data_uf.r) before running the Forecast Phase fits.", call. = FALSE)
    }
    message(sprintf(
      "Forecast Phase data readiness check (dengue | state | SP): max epiweek %d >= needed %d -- OK.",
      readiness$max_epiweek_available, readiness$needed_through_epiweek
    ))
  }, error = function(e) {
    if (grepl("^Forecast Phase data not ready", conditionMessage(e))) stop(e)
    message(
      "Could not run the Forecast Phase data-readiness check (", conditionMessage(e),
      ") -- proceeding without it; verify manually that processed_data covers through EW25 ",
      SUBMISSION_FORECAST_YEAR, " before trusting the Forecast Phase fits."
    )
  })
}

# ── Which target window(s) to actually fit this run ─────────────────────────
# By default every run refits and (re)writes predictions for all four
# validation splits (1, 2, 3, 4) -- always the same four fixed,
# already-elapsed (or currently in progress) seasons, per the challenge
# rules. The Forecast Phase is a separate, fifth thing entirely (see the
# config block above) that must be requested explicitly by name.
#
# TARGET_TO_FIT accepts either a plain numeric subset of c(1, 2, 3, 4), or a
# friendlier string spec (see resolve_target_spec() below):
#   "validation"                     -- all four validation splits (the full
#                                        4-fold retrospective backtest)
#   "validation_1".."validation_4"   -- just that one validation split
#   "forecast"                       -- the Forecast Phase submission (its
#                                        own target_forecast.csv -- a
#                                        DIFFERENT window and file from
#                                        Validation Test 4, never confused
#                                        with or overwriting it)
#
# Restricting this skips fitting (and leaves untouched on disk) whatever
# isn't included -- e.g. TARGET_TO_FIT <- "forecast" fits/writes only
# target_forecast.csv and leaves any existing target_1..4 CSVs alone.
# Validation Test 4 (target_4.csv) and the Forecast Phase (target_forecast.csv)
# are independent slots that can both exist on disk at the same time --
# fitting one never touches the other, and there is no config to switch
# between them (unlike the two phases sharing a slot, which is NOT how this
# works).
TARGET_TO_FIT <- "validation"

#' Resolve a target spec into the target_<n>/target_forecast split(s) to
#' fit/submit
#'
#' Duplicated (not sourced) between fit.r and sub_pred.r -- same reasoning as
#' resolve_case_file() in sub_pred.r: each file is meant to be runnable on
#' its own. Keep the two copies in sync if this changes.
#'
#' @param target  NULL (-> all four validation splits), a numeric subset of
#'   c(1,2,3,4), or one of "validation" (-> all four), "validation_1"..
#'   "validation_4" (-> that one split), or "forecast" (-> the separate
#'   Forecast Phase submission, its own "forecast" slot -- never 4).
#'
#' @return list(targets = <sorted integer vector, subset of 1:4, OR the
#'   single string "forecast">, label = <NULL, or a human-readable string
#'   for a single-target spec>).
resolve_target_spec <- function(target) {
  if (is.null(target)) target <- 1:4

  if (is.numeric(target)) {
    targets <- sort(unique(as.integer(target)))
    stopifnot(
      "target must be a non-empty subset of c(1, 2, 3, 4)" =
        length(targets) > 0 && all(targets %in% 1:4)
    )
    return(list(targets = targets, label = NULL))
  }

  if (!is.character(target) || length(target) != 1) {
    stop(
      "target must be NULL, a numeric subset of c(1,2,3,4), or one of ",
      "\"validation\", \"validation_1\"..\"validation_4\", \"forecast\"."
    )
  }

  spec <- tolower(trimws(target))

  if (spec == "validation") {
    return(list(targets = 1:4, label = NULL))
  }
  if (spec == "forecast") {
    return(list(targets = "forecast", label = "final forecast"))
  }
  m <- regmatches(spec, regexec("^validation_?([1-4])$", spec))[[1]]
  if (length(m) == 2) {
    n <- as.integer(m[2])
    return(list(targets = n, label = sprintf("Validation Test %d", n)))
  }

  stop(
    "Unrecognized target spec \"", target, "\". Use a numeric subset of ",
    "c(1,2,3,4), \"validation\", \"validation_1\"..\"validation_4\", or \"forecast\"."
  )
}

# Resolve + sanity-check TARGET_TO_FIT immediately -- before any of the long
# per-unit fitting loops below run -- exactly like validate_epiweek_window()
# above. Aborts right away on an unrecognized spec instead of partway
# through a long batch.
invisible(resolve_target_spec(TARGET_TO_FIT))

#' Resolve the train/forecast epiweek window for one fit/forecast split.
#'
#' Splits 1-4 are the four validation tests -- always the same fixed,
#' already-elapsed (or currently in progress) seasons, per the challenge
#' rules, computed identically regardless of anything else in this file.
#' "forecast" is the separate Forecast Phase window (FORECAST_WINDOW, see
#' the config block above) -- a different season/window from Validation
#' Test 4, never derived from it.
#'
#' @param split  One of 1, 2, 3, 4 (or "1".."4"), or the string "forecast".
#' @param data   The unit/disease data frame (used only to find the
#'               earliest epiweek actually present, e.g. chikungunya's
#'               later start).
#'
#' @return A list with train_start, train_end, forecast_start, forecast_end.
window_for_split <- function(split, data) {
  if (identical(split, "forecast")) return(FORECAST_WINDOW)
  expected_challenge_window("validation", test_id = as.integer(split), data_start = min(data$epiweek))
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
  pca_candidates <- candidates[!grepl("enso|iod|pdo", candidates)]

  # No covariates survived filtering (more common for municipalities, which
  # lack the lag/rolling-mean variants state files carry -- a smaller
  # candidate set to begin with) -- pca_all() calls prcomp()/svd() on a
  # zero-column matrix and errors ("a dimension is zero"). Mirrors the same
  # guard model_sel.r's own PCA branch already has (see its `pca_all()`
  # call): skip PCA and hand the data back untouched. Whatever formula
  # model selection picked for this unit is necessarily covariate-free
  # (e.g. `cases ~ 1`) if it also saw zero candidates here -- model
  # selection and this function build the same covariate space from the
  # same filtering steps.
  if (length(pca_candidates) == 0) {
    message("No covariates survived filtering -- skipping PCA (0 candidates).")
    return(data)
  }

  pca_result <- pca_all(
    data          = data,
    candidates    = pca_candidates,
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
# @param target       Which split(s) to fit and write a prediction CSV for
#                      (default NULL -> all four validation splits). Either
#                      a numeric subset of c(1,2,3,4), or a friendlier
#                      string spec -- "validation" (all four),
#                      "validation_1".."validation_4" (just that one split),
#                      or "forecast" (the separate Forecast Phase window --
#                      its own target_forecast.csv, never target_4) -- see
#                      `resolve_target_spec()` above. Splits not included
#                      are skipped entirely -- neither fit nor overwritten
#                      -- so a prior run's CSV for a skipped target is left
#                      exactly as it was.
#
# Each split's `fit_sarimax()` call passes this unit's persistent
# covariate-forecast table (`covariate_forecast_table_path()`) so a root
# series' forecast already computed by `model_sel.r`'s grid search for this
# exact unit/fold is reused here instead of being re-fit from scratch, and
# vice versa on a later `model_sel.r` re-run.
#
# @return A list: `id`, `resolved` (logical), `used_row` (1 = the top pick
#   worked untouched), `warnings` (named list by target_id), and `preds` (the
#   fitted splits' forecasts, row-bound, with a `target_id` and unit_col
#   column attached).
fit_unit_predictions <- function(disease, level, id, best_wis_df, unit_col,
                                  metrics_dir = "sarimax/results/metrics",
                                  preds_dir   = "sarimax/results/preds",
                                  target      = NULL) {

  spec       <- resolve_target_spec(target)
  targets    <- spec$targets
  target_ids <- paste0("target_", targets)

  message(sprintf(
    "[%s | %s | %s] Loading and preparing covariates... (target(s): %s%s)",
    disease, level, id, paste(target_ids, collapse = ", "),
    if (!is.null(spec$label)) paste0(" -- ", spec$label) else ""
  ))
  d <- load_unit_data(disease = disease, level = level, id = id)
  d <- prepare_unit_covariates(d)

  unit_row <- best_wis_df[best_wis_df[[unit_col]] == id, ]
  candidate_rows <- tibble(formula_id = unit_row$formula_id, order = unit_row$order)

  metrics_file <- file.path(metrics_dir, paste0("metrics_all_formulas_", disease, "_", id, ".csv"))
  if (file.exists(metrics_file)) {
    metrics_df <- read_leaderboard_csv(metrics_file)
    # Rows 2.. are fallback candidates, already ranked by CV score.
    if (nrow(metrics_df) > 1) {
      candidate_rows <- bind_rows(candidate_rows, metrics_df[-1, c("formula_id", "order")])
    }
  }

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

    pred_target <- lapply(seq_along(targets), function(j) {
      target_id <- target_ids[j]
      win <- window_for_split(targets[j], d)
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
#
# @param target  Which split(s) to fit for every unit in this
#                 (disease, level) -- passed straight through to
#                 `fit_unit_predictions()` (see its own `target` doc above):
#                 NULL/a numeric subset of c(1,2,3,4), or a friendly string
#                 spec ("validation", "validation_1".."validation_4",
#                 "forecast"). Resolved once here rather than once per unit.
run_all_units <- function(disease, level, target = NULL) {
  spec <- resolve_target_spec(target)

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

  best_wis_df <- read_leaderboard_csv(best_wis_path)
  ids <- best_wis_df[[unit_col]]

  results <- lapply(ids, function(id) {
    fit_unit_predictions(
      disease     = disease,
      level       = level,
      id          = id,
      best_wis_df = best_wis_df,
      unit_col    = unit_col,
      target      = spec$targets
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
# target = TARGET_TO_FIT applies the config block at the top of this file
# uniformly to every (disease, level) combination. Call run_all_units()
# directly with your own target spec (e.g. run_all_units("dengue", "state",
# "validation")) to override this for a one-off run instead.
preds_dengue_state      <- run_all_units("dengue",      "state", TARGET_TO_FIT)
preds_chikungunya_state <- run_all_units("chikungunya", "state", TARGET_TO_FIT)
preds_dengue_city       <- run_all_units("dengue",      "city",  TARGET_TO_FIT)
preds_chikungunya_city  <- run_all_units("chikungunya", "city",  TARGET_TO_FIT)

# ── Forecast Phase (challenge rules, Section 5.2) ───────────────────────────
# A separate, additional run for every disease x level combination -- NOT
# part of TARGET_TO_FIT above, and never a substitute for it. Fits
# FORECAST_WINDOW (the 2026-2027 season, configured via
# SUBMISSION_FORECAST_YEAR near the top of this file) and writes each unit's
# own target_forecast.csv, alongside its target_1..4.csv from the validation
# run above -- neither run touches or overwrites the other's files.
check_forecast_phase_data_ready()
preds_dengue_state_forecast      <- run_all_units("dengue",      "state", "forecast")
preds_chikungunya_state_forecast <- run_all_units("chikungunya", "state", "forecast")
preds_dengue_city_forecast       <- run_all_units("dengue",      "city",  "forecast")
preds_chikungunya_city_forecast  <- run_all_units("chikungunya", "city",  "forecast")
