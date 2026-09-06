# sarimax/src/validate_challenge_rules.r
#
# Validates that this pipeline's train/target epiweek windows follow the
# IMDC 2026 data-usage rule (challenge rules, Section 5 -- "Validation
# Phase" / "Forecast Phase"):
#
#   Validation Test i (i = 1..4): forecast weekly cases for the
#   (2021+i)-(2022+i) season [EW41 (2021+i) .. EW40 (2022+i)], using data
#   from EW01 2010 to EW25 (2021+i).
#
#   Forecast Phase: forecast weekly cases for the forecast_year/
#   (forecast_year+1) season [EW41 forecast_year .. EW40 (forecast_year+1)],
#   using all available data from EW01 2010 to EW25 forecast_year.
#
# See root README.md, Section 6 ("Data Usage Restriction") for how this
# repo otherwise documents the rule, and sarimax/src/fit.r for where the
# train_4/target_4 window is actually fit (the unified fit_sarimax() in
# sarimax/src/utils.r, called via window_for_train_id()).
#
# This script is self-contained (base R only -- no tidyverse dependency) so
# it can be sourced and run on its own, in a fresh R session, or in CI:
#
#     source("sarimax/src/validate_challenge_rules.r")
#     run_challenge_validation("dengue")
#     run_challenge_validation("chikungunya")
#
# It does NOT fit any models and does NOT modify any files.

# -- epiweek helpers (kept in sync with the copies in sarimax/src/fit.r) ----
# Duplicated here (rather than sourced from fit.r) because fit.r also
# contains top-level code that runs the full fitting/forecasting pipeline as
# a side effect of being sourced. If these helpers are ever promoted into
# sarimax/src/utils.r (which only defines functions), this block can be
# deleted in favor of `source("sarimax/src/utils.r")`.

epiweek_to_date <- function(yw) {
  yr   <- yw %/% 100L
  wk   <- yw  %% 100L
  jan4 <- as.Date(paste0(yr, "-01-04"))
  dow_jan4  <- as.integer(format(jan4, "%w"))
  sunday_w1 <- jan4 - dow_jan4
  sunday_w1 + (wk - 1L) * 7L
}

has_53_weeks <- function(year) {
  epiweek_to_date(year * 100L + 53L) < epiweek_to_date((year + 1L) * 100L + 1L)
}

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

# -- 1. What the rules require -----------------------------------------------

#' Expected train/forecast epiweek window for one challenge phase
#'
#' @param phase          "validation" or "forecast".
#' @param test_id        Integer 1..4, required when phase == "validation".
#' @param forecast_year  Integer year (e.g. 2026), required when
#'                       phase == "forecast".
#' @param data_start     Earliest epiweek (YYYYWW) actually present in the
#'                       dataset being checked. The rules allow data back to
#'                       EW01 2010, but a disease whose reporting only
#'                       starts later (e.g. chikungunya, tracked from 2014)
#'                       legitimately has a later train_start -- the
#'                       expected train_start is
#'                       max(201001, data_start), not a hard 201001.
#'                       Default 201001 (i.e. assume data goes back that
#'                       far unless told otherwise).
#'
#' @return A list with train_start, train_end, forecast_start, forecast_end
#'   (YYYYWW integers) and a human-readable label.
expected_challenge_window <- function(phase = c("validation", "forecast"),
                                       test_id = NULL,
                                       forecast_year = NULL,
                                       data_start = 201001L) {
  phase <- match.arg(phase)

  if (phase == "validation") {
    if (is.null(test_id) || !(test_id %in% 1:4)) {
      stop("test_id must be one of 1, 2, 3, 4 when phase == 'validation'")
    }
    season_year <- 2021L + as.integer(test_id)
    label <- sprintf(
      "Validation Test %d (%d-%d season)", test_id, season_year, season_year + 1L
    )
  } else {
    if (is.null(forecast_year)) {
      stop("forecast_year is required when phase == 'forecast'")
    }
    season_year <- as.integer(forecast_year)
    label <- sprintf(
      "Forecast Phase (%d-%d season)", season_year, season_year + 1L
    )
  }

  list(
    train_start    = max(201001L, as.integer(data_start)),
    train_end      = season_year * 100L + 25L,
    forecast_start = season_year * 100L + 41L,
    forecast_end   = (season_year + 1L) * 100L + 40L,
    label          = label
  )
}

# -- 2. Guard for an explicit fit_sarimax() call -----------------------------

#' Validate an explicit train/forecast epiweek window against the rules
#'
#' Meant to guard a `fit_sarimax()` call in sarimax/src/fit.r: pass
#' the exact train_start/train_end/forecast_start/forecast_end you are
#' about to fit with, and this checks them against what the rules require
#' for that phase/test. With `stop_on_fail = TRUE` (default) it aborts with
#' an informative error instead of letting an out-of-window forecast get
#' written/submitted; pass FALSE to just get the result back (e.g. for a
#' report).
#'
#' @return Invisibly, a list with `valid`, `label`, `expected`, `actual`,
#'   and `mismatches` (character vector, empty if valid).
validate_epiweek_window <- function(train_start, train_end,
                                     forecast_start, forecast_end,
                                     phase = c("validation", "forecast"),
                                     test_id = NULL,
                                     forecast_year = NULL,
                                     data_start = 201001L,
                                     stop_on_fail = TRUE) {
  phase <- match.arg(phase)
  expected <- expected_challenge_window(phase, test_id, forecast_year, data_start)
  actual <- list(
    train_start = as.integer(train_start), train_end = as.integer(train_end),
    forecast_start = as.integer(forecast_start), forecast_end = as.integer(forecast_end)
  )

  fields <- c("train_start", "train_end", "forecast_start", "forecast_end")
  mismatches <- character(0)
  for (f in fields) {
    if (!identical(actual[[f]], as.integer(expected[[f]]))) {
      mismatches <- c(mismatches, sprintf(
        "%-15s expected %d, got %d", paste0(f, ":"), expected[[f]], actual[[f]]
      ))
    }
  }

  result <- list(
    valid      = length(mismatches) == 0,
    label      = expected$label,
    expected   = expected,
    actual     = actual,
    mismatches = mismatches
  )

  if (!result$valid && stop_on_fail) {
    stop(sprintf(
      "Epiweek window does not match the IMDC rules for %s:\n  %s",
      expected$label, paste(mismatches, collapse = "\n  ")
    ), call. = FALSE)
  }

  invisible(result)
}

# -- 3. Check the train_i/target_i flag columns in a processed data file ----

#' Validate one state/disease file's train_<test_id>/target_<test_id> flags
#'
#' Checks that `train_<test_id>` is TRUE for exactly the epiweeks the rules
#' allow (data_start..EW25 of the test's season, no gaps, nothing beyond),
#' and that `target_<test_id>` is TRUE only for epiweeks inside
#' EW41..EW40 of that season. `target_<test_id>` is allowed to be an
#' incomplete *prefix* of that window (real case counts for future weeks of
#' an ongoing season don't exist yet) but must never skip a week in the
#' middle or contain a week outside the allowed window.
#'
#' @param data         Data frame with an epiweek column and
#'                     `train_<test_id>`/`target_<test_id>` boolean/0-1
#'                     columns (e.g. one row of
#'                     processed_data/<disease>/<disease>_<UF>_agg.csv.gz).
#' @param test_id      Integer 1..4.
#' @param epiweek_col  Name of the epiweek column (default "epiweek").
#'
#' @return A one-row data.frame summarizing pass/fail for the train and
#'   target windows.
validate_train_target_flags <- function(data, test_id, epiweek_col = "epiweek") {
  if (!(test_id %in% 1:4)) stop("test_id must be one of 1, 2, 3, 4")
  train_col  <- paste0("train_",  test_id)
  target_col <- paste0("target_", test_id)
  needed <- c(epiweek_col, train_col, target_col)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    stop("data is missing column(s): ", paste(missing_cols, collapse = ", "))
  }

  ew <- as.integer(data[[epiweek_col]])
  is_true <- function(x) x %in% c(TRUE, 1, "1", "TRUE", "true")

  actual_train  <- sort(ew[is_true(data[[train_col]])])
  actual_target <- sort(ew[is_true(data[[target_col]])])

  expected <- expected_challenge_window(
    "validation", test_id = test_id, data_start = min(ew)
  )
  expected_train  <- enumerate_epiweeks(expected$train_start, expected$train_end)
  expected_target <- enumerate_epiweeks(expected$forecast_start, expected$forecast_end)

  train_missing <- setdiff(expected_train, actual_train)
  train_extra   <- setdiff(actual_train, expected_train)

  target_extra <- setdiff(actual_target, expected_target)
  target_is_valid_prefix <- length(actual_target) == 0 ||
    (length(actual_target) <= length(expected_target) &&
       identical(actual_target, expected_target[seq_along(actual_target)]))

  data.frame(
    test_id           = test_id,
    label             = expected$label,
    train_ok          = length(train_missing) == 0 && length(train_extra) == 0,
    train_missing_n   = length(train_missing),
    train_extra_n     = length(train_extra),
    target_ok         = length(target_extra) == 0 && target_is_valid_prefix,
    target_extra_n    = length(target_extra),
    target_n          = length(actual_target),
    target_expected_n = length(expected_target),
    target_complete   = length(actual_target) == length(expected_target),
    stringsAsFactors  = FALSE
  )
}

#' Run validate_train_target_flags() for all 4 validation tests
validate_all_tests <- function(data, epiweek_col = "epiweek") {
  do.call(rbind, lapply(1:4, function(i) {
    validate_train_target_flags(data, i, epiweek_col)
  }))
}

# -- 4. Run the check across every state file for a disease ------------------

#' Run the train_i/target_i validity check across every processed state file
#'
#' @param disease         "dengue" or "chikungunya".
#' @param processed_dir   Root of the processed_data/ tree (default
#'                        "processed_data" -- run from the repo root).
#' @param states          Character vector of UF codes to check; NULL
#'                        (default) checks every
#'                        `<disease>_<UF>_agg.csv.gz` file found.
#'
#' @return A data.frame, one row per state x test_id, with a `passed`
#'   summary column. Also prints a one-line summary and any failing rows.
run_challenge_validation <- function(disease, processed_dir = "processed_data", states = NULL) {
  dir <- file.path(processed_dir, disease)
  pattern <- paste0("^", disease, "_[A-Z]{2}_agg\\.csv\\.gz$")
  files <- list.files(dir, pattern = pattern, full.names = TRUE)

  if (!is.null(states)) {
    state_of <- sub(paste0("^", disease, "_(.*)_agg\\.csv\\.gz$"), "\\1", basename(files))
    files <- files[state_of %in% states]
  }
  if (length(files) == 0) {
    stop("No '", disease, "_<UF>_agg.csv.gz' files found under ", dir)
  }

  rows <- lapply(files, function(f) {
    st <- sub(paste0("^", disease, "_(.*)_agg\\.csv\\.gz$"), "\\1", basename(f))
    d  <- read.csv(gzfile(f), stringsAsFactors = FALSE)
    res <- validate_all_tests(d)
    cbind(disease = disease, state = st, res, stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)
  out$passed <- out$train_ok & out$target_ok

  n_fail <- sum(!out$passed)
  message(sprintf(
    "[%s] checked %d state file(s) x 4 validation tests = %d rows, %d failed.",
    disease, length(files), nrow(out), n_fail
  ))
  if (n_fail > 0) {
    print(out[!out$passed, c("state", "test_id", "label", "train_ok", "target_ok",
                              "train_missing_n", "train_extra_n", "target_extra_n")])
  }

  out
}

# -- 5. Is the data actually updated for the Forecast Phase yet? -------------

#' Check whether a data file has been updated far enough for the Forecast
#' Phase (Section 5.2: data update on July 31, using data through EW25 of
#' forecast_year). This does NOT check that fit.r's hardcoded window has
#' been updated to match -- use validate_epiweek_window() for that.
#'
#' @return A one-row data.frame: max epiweek available vs. the epiweek
#'   needed, and whether the data is ready.
validate_forecast_phase_data <- function(data, forecast_year, epiweek_col = "epiweek") {
  ew <- as.integer(data[[epiweek_col]])
  needed <- as.integer(forecast_year) * 100L + 25L
  data.frame(
    forecast_year          = forecast_year,
    max_epiweek_available  = max(ew),
    needed_through_epiweek = needed,
    data_ready             = max(ew) >= needed
  )
}
