source("sarimax/src/utils.r")   # also sources forecast_covariates.r
source("sarimax/src/validate_challenge_rules.r")

# ── Why this file exists ─────────────────────────────────────────────────────
# forecast_covariates()/get_or_compute_root_forecast() (sarimax/src/
# forecast_covariates.r) already cache every root covariate forecast they
# compute into a persistent table (processed_data/covariate_forecasts_
# <disease>_<state|city>.csv), so model_sel.r, fit.r, and
# evaluate_covariate_forecasts.r automatically reuse whatever was fit before
# for a given (unit, fold) instead of re-fitting. But that cache only ever
# gets warmed lazily, as a side effect of running the rest of the pipeline --
# the first person/run to touch a given unit/fold pays the fitting cost.
#
# This script does the same job eagerly and exhaustively: for every unit
# under processed_data/ (not just the ones that have already made it into a
# best_wis_*.csv leaderboard) and every validation split, it forecasts each
# root covariate once and lets it land in the shared table. Run it once
# up front (e.g. after a data refresh) and every later model_sel.r/fit.r/
# evaluate_covariate_forecasts.r call for any unit becomes a pure table hit.
#
# This is a precompute/warm-up utility, not part of the submission pipeline
# itself -- nothing here fits a disease model or touches best_wis_*.csv.

# ── 1. Enumerate every unit under processed_data/ ───────────────────────────

#' List every unit id with a processed data file for one disease/level
#'
#' State-level files are `processed_data/<disease>/<disease>_<UF>_agg.csv.gz`;
#' municipality-level files are `processed_data/<disease>/sel_cities/
#' <disease>_<geocode>.csv.gz`. This just lists whichever directory applies
#' and extracts the id from the filename -- it does not open the files.
#'
#' @param disease "dengue" or "chikungunya".
#' @param level   "state" or "city".
#' @param processed_data_dir Root of the processed_data tree (default
#'   "processed_data", i.e. run this from the repo root).
#'
#' @return Character vector of unit ids (UF codes or geocodes), possibly
#'   empty if the expected directory doesn't exist or has no matching files.
enumerate_processed_units <- function(disease, level = c("state", "city"),
                                       processed_data_dir = "processed_data") {
  level <- match.arg(level)

  if (level == "state") {
    dir_path <- file.path(processed_data_dir, disease)
    pattern  <- paste0("^", disease, "_(.+)_agg\\.csv\\.gz$")
  } else {
    dir_path <- file.path(processed_data_dir, disease, "sel_cities")
    pattern  <- paste0("^", disease, "_(.+)\\.csv\\.gz$")
  }

  if (!dir.exists(dir_path)) {
    message(sprintf("enumerate_processed_units(): %s does not exist -- returning 0 units.", dir_path))
    return(character(0))
  }

  files <- list.files(dir_path, pattern = pattern)
  ids   <- sub(pattern, "\\1", files)
  sort(unique(ids))
}

# ── 2. One (unit, split) precompute task ────────────────────────────────────

# The 8 "root" candidate series (see get_candidates()'s `contemporaneous`
# group): everything else forecast_covariates() is ever asked for (the
# _lag{4,8,12,16} and _mean_{3,6,9,12}mo derived columns; see
# parse_derived_covariate()) is a deterministic transform of one of these,
# so precomputing forecasts for just these 8 per unit/fold is sufficient to
# make every derived-column forecast in forecast_covariates() a table hit
# too.
default_root_vars <- function() {
  c("temp_med_mean", "precip_med_mean", "rel_humid_med_mean",
    "thermal_range_mean", "rainy_days_mean", "enso", "iod", "pdo")
}

#' Precompute and cache root-covariate forecasts for one unit's one split
#'
#' Loads the unit's data, resolves that split's train/forecast window via
#' `expected_challenge_window()`, and calls `forecast_covariates()` with
#' `rhs_terms` set to whichever of `vars` are present in this unit's data --
#' `forecast_covariates()` does the actual fitting/caching via
#' `get_or_compute_root_forecast()`, so this function's only job is to ask
#' for the right window. Units with too little history for a given split
#' (e.g. a municipality whose data starts after that split's train_start)
#' are skipped with a message, not an error, so one bad unit/fold doesn't
#' abort a whole precompute run.
#'
#' @param id         Unit id (UF code or municipality geocode).
#' @param disease    "dengue" or "chikungunya".
#' @param level      "state" or "city".
#' @param split      Which of the 4 validation splits (1-4) to precompute.
#' @param vars       Root covariates to precompute (default `default_root_vars()`).
#' @param table_path Persistent table path (default `covariate_forecast_table_path(disease, level)`).
#'
#' @return A list with `ok` (logical) and `reason` (`NA_character_` on
#'   success, a short description of why otherwise) -- returned as a value
#'   rather than left as only a `message()` call, because `message()`
#'   output from a parallel worker is silently discarded by default (PSOCK
#'   workers' `outfile` is `/dev/null` unless told otherwise), which would
#'   otherwise leave a `n_cores > 1` run with no way to see *why* tasks
#'   were skipped -- only that they were.
precompute_unit_split <- function(id, disease, level, split,
                                   vars       = default_root_vars(),
                                   table_path = covariate_forecast_table_path(disease, level)) {
  tryCatch({
    d <- load_unit_data(disease = disease, level = level, id = id)

    win <- expected_challenge_window("validation", test_id = split, data_start = min(d$epiweek))

    root_vars <- intersect(vars, names(d))
    if (length(root_vars) == 0) {
      reason <- "none of the requested root vars are present"
      message(sprintf("[%s | %s | %s | split %d] %s -- skipping.",
                       disease, level, id, split, reason))
      return(list(ok = FALSE, reason = reason))
    }

    h <- length(enumerate_epiweeks(win$forecast_start, win$forecast_end))

    forecast_covariates(
      data        = d,
      rhs_terms   = root_vars,
      train_start = win$train_start,
      train_end   = win$train_end,
      h           = h,
      table_path  = table_path,
      id          = id
    )
    list(ok = TRUE, reason = NA_character_)
  }, error = function(e) {
    reason <- conditionMessage(e)
    message(sprintf(
      "[%s | %s | %s | split %d] Precompute failed, skipping: %s",
      disease, level, id, split, reason
    ))
    list(ok = FALSE, reason = reason)
  })
}

# Top-level named worker (not an anonymous closure) so it and everything it
# calls can be handed to clusterExport() by name -- same reasoning as
# eval_one_unit_covariate_forecasts() in evaluate_covariate_forecasts.r.
precompute_one_task <- function(task, disease, level, vars, table_path) {
  precompute_unit_split(
    id = task$id, disease = disease, level = level, split = task$split,
    vars = vars, table_path = table_path
  )
}

# ── 3. Bulk driver: every unit x every split, optionally in parallel ────────

#' Precompute and cache root-covariate forecasts for every unit under
#' processed_data/, across every requested split
#'
#' @param disease  "dengue" or "chikungunya".
#' @param level    "state" or "city".
#' @param ids      Unit ids to precompute. Default: every unit with a file
#'                 under `processed_data/` (see `enumerate_processed_units()`)
#'                 -- not just units that already appear in a
#'                 `best_wis_*.csv` leaderboard, unlike
#'                 `evaluate_all_units_covariate_forecasts()`.
#' @param splits   Integer vector of validation splits to precompute
#'                 (default 1:4).
#' @param vars     Root covariates to precompute (default `default_root_vars()`).
#' @param n_cores  Worker processes for a real PSOCK cluster; `1` runs
#'                 sequentially. Default: all but one detected core.
#'
#' @return Invisibly, a tibble with one row per (id, split) task, a logical
#'   `ok` column (`TRUE` = precomputed or already cached, `FALSE` =
#'   skipped), and a `reason` column (`NA` on success, else why it was
#'   skipped) -- printed as a grouped summary at the end regardless of
#'   `n_cores`, so failures are visible even when the per-task `message()`s
#'   themselves were swallowed by a parallel worker's suppressed output.
build_covariate_forecast_table <- function(disease, level, ids = NULL,
                                            splits  = 1:4,
                                            vars    = default_root_vars(),
                                            n_cores = max(1L, detectCores() - 1L)) {
  if (is.null(ids)) ids <- enumerate_processed_units(disease, level)

  if (length(ids) == 0) {
    message(sprintf("[%s | %s] No units found under processed_data/ -- nothing to precompute.",
                     disease, level))
    return(invisible(tibble(id = character(0), split = integer(0),
                             ok = logical(0), reason = character(0))))
  }

  table_path <- covariate_forecast_table_path(disease, level)
  tasks <- expand.grid(id = ids, split = splits, stringsAsFactors = FALSE)
  tasks <- tasks[order(tasks$id, tasks$split), , drop = FALSE]
  task_list <- lapply(seq_len(nrow(tasks)), function(i) as.list(tasks[i, ]))

  n_cores <- max(1L, min(n_cores, length(task_list)))

  message(sprintf(
    "[%s | %s] Precomputing %d unit(s) x %d split(s) = %d task(s) into %s ...",
    disease, level, length(ids), length(splits), length(task_list), table_path
  ))

  cl <- NULL
  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)

    clusterExport(cl, varlist = c(
      "precompute_one_task", "precompute_unit_split", "default_root_vars",
      "load_unit_data", "forecast_covariates", "forecast_covariate_series",
      "parse_derived_covariate", "enumerate_epiweeks", "expected_challenge_window",
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
    })
  }

  results <- pbapply::pblapply(task_list, precompute_one_task,
                                disease = disease, level = level,
                                vars = vars, table_path = table_path, cl = cl)

  tasks$ok     <- vapply(results, function(r) isTRUE(r$ok), logical(1))
  tasks$reason <- vapply(results, function(r) {
    if (is.null(r$reason) || is.na(r$reason)) NA_character_ else r$reason
  }, character(1))

  n_done <- sum(tasks$ok)
  message(sprintf(
    "[%s | %s] Done: %d/%d task(s) precomputed or already cached, %d skipped.",
    disease, level, n_done, nrow(tasks), nrow(tasks) - n_done
  ))
  if (n_done < nrow(tasks)) {
    # Printed here (not just left to each task's own message()) so the
    # reasons are visible even under n_cores > 1, where a worker's
    # message() output is otherwise silently discarded -- see the
    # precompute_unit_split() doc comment.
    reason_counts <- sort(table(tasks$reason[!tasks$ok]), decreasing = TRUE)
    message(sprintf("[%s | %s] Skip reasons:", disease, level))
    for (r in names(reason_counts)) {
      message(sprintf("  %d x  %s", reason_counts[[r]], r))
    }
  }
  invisible(as_tibble(tasks))
}

# ── Run for every disease x level ────────────────────────────────────────
precompute_dengue_state      <- build_covariate_forecast_table("dengue",      "state")
precompute_chikungunya_state <- build_covariate_forecast_table("chikungunya", "state")
precompute_dengue_city       <- build_covariate_forecast_table("dengue",      "city")
precompute_chikungunya_city  <- build_covariate_forecast_table("chikungunya", "city")
