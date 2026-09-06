# sarimax/src/forecast_covariates.r
#
# Produces genuine, out-of-sample forecasts for the SARIMAX regressor
# columns (PCA components and/or surviving ocean-climate indices) needed to
# forecast the EW41..EW40 target window in sarimax/src/fit.r's
# fit_sarimax_epiweek().
#
# Why this exists: the real future values of these covariates are not known
# at submission time (the weather/ocean indices for the target season
# haven't happened yet), and the IMDC rules require the whole forecast --
# not just the case-count model -- to be produced using only data available
# up to EW25 of the training year (see root README.md, Section 6, "Data
# Usage Restriction"). Earlier versions of this pipeline substituted the
# previous year's realized value at the same epiweek as a stand-in; that is
# simple but is not itself a forecast, and reusing a single past year's
# realized trajectory (rather than an estimate built from the full history)
# risks looking like it's smuggling in real observed data rather than a
# genuine out-of-sample prediction.
#
# An earlier version of this module fit each regressor its own univariate
# seasonal time-series model (forecast::tbats(), falling back to stlf() and
# a differenced AR(1)-with-drift Arima()). At the scale of the real
# pipeline -- every unit x every validation split x every root covariate --
# that turned out to be the wrong tool: tbats()/stlf() routinely took many
# seconds per series on real (decade-plus) histories, and a handful of
# slow/non-converging fits were enough to make
# build_covariate_forecast_table.r's precompute run mostly time out (see
# the note on .with_covariate_forecast_table_lock() below for the specific
# failure mode that produced). It also wasn't buying much: weather and
# ocean-climate covariates are strongly seasonal and only weakly driven by
# their own most recent trajectory, so a full ARIMA-family fit was mostly
# rediscovering the seasonal cycle a much simpler estimate already captures.
#
# This module now instead uses a seasonal climatology: for each future
# epiweek, the mean of that exact calendar week (week-of-year) over the
# most recent few years of the series' own history (see
# forecast_covariate_series() below) -- computed only from observations up
# to the same train_end cutoff the disease model itself respects, so it's
# still a genuine out-of-sample estimate, just a much simpler and far more
# robust one. This has no numerical-convergence failure mode, is orders of
# magnitude faster, and no longer needs the `forecast` package at all for
# covariate forecasting (only base R).

#' Round-trip an MMWR/Brazilian epiweek (YYYYWW) to/from its opening Sunday,
#' and step forward one epiweek at a time
#'
#' Self-contained duplicate of the same logic in `sarimax/src/utils.r`
#' (`epiweek_to_date()` / `has_53_weeks()` / `enumerate_epiweeks()`), kept
#' here under dotted names so this file's forecasting functions stay usable
#' -- and testable -- without sourcing the rest of the pipeline, matching
#' this file's existing "depends on nothing else in this repo" design.
.epiweek_to_date <- function(yw) {
  yr   <- yw %/% 100L
  wk   <- yw  %% 100L
  jan4 <- as.Date(paste0(yr, "-01-04"))
  dow_jan4  <- as.integer(format(jan4, "%w"))  # %w: 0 = Sunday
  sunday_w1 <- jan4 - dow_jan4                 # Sunday on or before Jan 4
  sunday_w1 + (wk - 1L) * 7L
}

#' Whether `year` has 53 epiweeks (MMWR/Brazilian calendar) rather than 52
.has_53_weeks <- function(year) {
  .epiweek_to_date(year * 100L + 53L) < .epiweek_to_date((year + 1L) * 100L + 1L)
}

#' The epiweek immediately after `yw`, correctly rolling over the year
#' boundary and accounting for 52- vs 53-week years
.next_epiweek <- function(yw) {
  yr <- yw %/% 100L
  wk <- yw  %% 100L
  n_weeks <- if (.has_53_weeks(yr)) 53L else 52L
  if (wk < n_weeks) yr * 100L + (wk + 1L) else (yr + 1L) * 100L + 1L
}

#' The `h` epiweeks immediately following `start_yw`, in order
.epiweek_seq_forward <- function(start_yw, h) {
  out <- integer(h)
  cur <- start_yw
  for (i in seq_len(h)) {
    cur <- .next_epiweek(cur)
    out[i] <- cur
  }
  out
}

#' Forecast a single covariate forward via seasonal (epiweek-of-year)
#' climatology
#'
#' For each future epiweek, averages that exact calendar week (week-of-
#' year) across the most recent `n_years` years found in `y`'s own history
#' -- e.g. the forecast for epiweek 202642 is the mean of week 42 across
#' however many of the last 5 years of `y` have a week 42.
#'
#' @param y         Numeric vector: the covariate's history, in the same
#'                   order as `epiweeks`, up to (and not beyond) the
#'                   allowed training cutoff.
#' @param epiweeks  Integer vector (YYYYWW), same length as `y`. Need not
#'                   be gap-free or presorted -- both are handled.
#' @param h         Integer: number of steps ahead to forecast.
#' @param n_years   How many of the most recent years with a matching
#'                   week-of-year to average per forecast step (default 5).
#'                   Fewer are used if that many aren't available. NA
#'                   occurrences are dropped rather than propagated; if the
#'                   target week-of-year has no usable (non-NA) match
#'                   anywhere in the history -- only possible for week 53,
#'                   which not every year has -- week 52 is used instead,
#'                   and if a series has no usable data at all, every step
#'                   falls back to that series' overall mean.
#'
#' @return Numeric vector of length `h`: the point forecast for each of the
#'   `h` epiweeks immediately following `max(epiweeks)`.
forecast_covariate_series <- function(y, epiweeks, h, n_years = 5) {
  y        <- as.numeric(y)
  epiweeks <- as.integer(epiweeks)
  if (length(y) != length(epiweeks)) {
    stop("forecast_covariate_series(): 'y' and 'epiweeks' must be the same length.")
  }
  if (length(y) == 0) {
    stop("forecast_covariate_series(): 'y' has no observations to forecast from.")
  }

  ord      <- order(epiweeks)
  y        <- y[ord]
  epiweeks <- epiweeks[ord]
  week_of  <- epiweeks %% 100L

  target_epiweeks <- .epiweek_seq_forward(max(epiweeks), h)
  overall_mean    <- mean(y, na.rm = TRUE)

  candidates_for <- function(wk) {
    idx <- week_of == wk
    yw  <- epiweeks[idx]
    val <- y[idx]
    keep <- !is.na(val)
    list(yw = yw[keep], y = val[keep])
  }

  vapply(target_epiweeks, function(target_yw) {
    wk   <- target_yw %% 100L
    cand <- candidates_for(wk)

    if (length(cand$y) == 0 && wk == 53L) {
      cand <- candidates_for(52L)
    }
    if (length(cand$y) == 0) {
      return(overall_mean)
    }

    top <- order(cand$yw, decreasing = TRUE)[seq_len(min(n_years, length(cand$y)))]
    mean(cand$y[top])
  }, numeric(1))
}

#' Recognize a lag/rolling-mean covariate as a deterministic transform of
#' one "base" (contemporaneous) series
#'
#' `data_prep/agg_data_uf.r` derives, for every base weather aggregate
#' (e.g. `temp_med_mean`), eight further columns that are pure, fixed
#' functions of that ONE series: `<base>_lag{4,8,12,16}` (`dplyr::lag(x, n)`
#' -- a backward shift) and `<base>_mean_{3,6,9,12}mo` (a trailing
#' `runner::mean_run()` rolling mean over 12/24/36/48 weeks, itself lagged
#' one more week to avoid leakage: `dplyr::lag(mean_run(x, k), 1)`). They
#' are not independent covariates -- they're different fixed views of one
#' underlying series -- so `forecast_covariates()` uses this to forecast
#' each base series only once and mechanically re-derive the rest (see its
#' own doc comment for why fitting all nine independently is wrong, not
#' just wasteful).
#'
#' @param term  A candidate column name.
#'
#' @return NULL if `term` doesn't match the naming convention above
#'   (contemporaneous base series, a PCA component, ENSO/IOD/PDO, or
#'   anything else -- all forecast directly, unchanged). Otherwise a list
#'   with `base` (the underlying series' column name), `kind` (`"lag"` or
#'   `"roll"`), and either `n` (weeks to shift back) or `k` (trailing
#'   window width in weeks, i.e. months * 4).
parse_derived_covariate <- function(term) {
  m_lag <- regmatches(term, regexec("^(.+)_lag(4|8|12|16)$", term))[[1]]
  if (length(m_lag) == 3) {
    return(list(base = m_lag[2], kind = "lag", n = as.integer(m_lag[3])))
  }
  m_roll <- regmatches(term, regexec("^(.+)_mean_(3|6|9|12)mo$", term))[[1]]
  if (length(m_roll) == 3) {
    return(list(base = m_roll[2], kind = "roll", k = as.integer(m_roll[3]) * 4L))
  }
  NULL
}

#' Path to the persistent covariate-forecast table for one (disease, level)
#'
#' A single long-format CSV shared by every unit of that disease/level --
#' `model_sel.r`, `fit.r`, and `evaluate_covariate_forecasts.r` all read
#' from (and, on a miss, write to) the same file via this one naming
#' convention, so a forecast computed by any one of them is available to
#' the others without recomputing. Lives under `processed_data/` (not
#' `sarimax/results/`) alongside the disease/climate input tables it's
#' built from, since -- once populated -- it's read the same way those are:
#' a precomputed input the modeling scripts just load, not a per-run output.
#'
#' @param disease  "dengue" or "chikungunya".
#' @param level    "state" or "city".
#' @param dir      Base directory (default "processed_data").
#'
#' @return The file path (not guaranteed to exist yet).
covariate_forecast_table_path <- function(disease, level, dir = "processed_data") {
  file.path(dir, paste0("covariate_forecasts_", disease, "_", level, ".csv"))
}

#' Column schema of the persistent covariate-forecast table
#'
#' One row per (unit, root series, fold window, forecast step). Kept as a
#' named constant so every reader/writer agrees on column names/order.
.covariate_forecast_table_cols <- c(
  "id", "root_var", "train_start", "train_end", "h", "n_years", "step", "forecast_value"
)

#' Read the persistent covariate-forecast table, or an empty one if it
#' doesn't exist yet -- or can't be read cleanly right now
#'
#' The empty-table shape is also what's returned if `table_path` exists but
#' fails to parse (e.g. `read.csv()`'s "no lines available in input" on a
#' zero-byte file). That specific failure is expected, not exceptional: the
#' lock-free fast-path read in `get_or_compute_root_forecast()` reads this
#' file without acquiring the write lock (deliberately -- see the note on
#' `.with_covariate_forecast_table_lock()` for why holding a lock across
#' anything slow is the thing we're avoiding here), so it can observe the
#' table mid-write by another worker. `get_or_compute_root_forecast()`'s
#' write path writes to a temp file and renames it into place precisely so
#' a reader never sees a *partially written* file -- but the very first
#' write of a brand-new table still has one unavoidable instant where the
#' path exists-but-is-empty (between the OS creating the directory entry
#' and the rename completing), and on a slow/networked filesystem (this
#' pipeline's real deployment target is a OneDrive-synced folder) that
#' instant is wide enough for many simultaneous workers to actually hit it.
#' Treating any read failure here as "nothing usable yet" -- rather than
#' propagating the error -- is safe because the caller always falls back to
#' computing fresh (and the locked write below re-checks for a hit before
#' persisting, so no work or data is lost, just occasionally redone).
.read_covariate_forecast_table <- function(table_path) {
  empty <- as.data.frame(setNames(
    replicate(length(.covariate_forecast_table_cols), character(0), simplify = FALSE),
    .covariate_forecast_table_cols
  ), stringsAsFactors = FALSE)

  if (!file.exists(table_path)) {
    return(empty)
  }
  tryCatch(
    utils::read.csv(table_path, stringsAsFactors = FALSE, colClasses = "character"),
    error = function(e) empty
  )
}

#' Write `tbl` to `table_path` without ever exposing a reader to a
#' partially-written file
#'
#' Writes to a uniquely-named temp file in the same directory (so the
#' rename below stays on one filesystem/volume) and renames it into place.
#' `file.rename()` is atomic on POSIX even when the destination already
#' exists; on Windows it can fail in that case, so this falls back to
#' removing the old file first (a much narrower race than writing in place
#' ever was -- see `.read_covariate_forecast_table()` -- and one that
#' degrades to a harmless re-read-as-empty rather than a parse error).
.write_covariate_forecast_table <- function(tbl, table_path) {
  dir.create(dirname(table_path), showWarnings = FALSE, recursive = TRUE)
  tmp_path <- file.path(
    dirname(table_path),
    sprintf(".%s.tmp-%d-%s", basename(table_path), Sys.getpid(),
            paste(sample(c(letters, 0:9), 8, replace = TRUE), collapse = ""))
  )
  utils::write.csv(tbl, tmp_path, row.names = FALSE)
  if (!isTRUE(file.rename(tmp_path, table_path))) {
    unlink(table_path)
    if (!isTRUE(file.rename(tmp_path, table_path))) {
      # Last resort: still correct (we hold the write lock throughout), just
      # no longer atomic with respect to a lock-free reader.
      file.copy(tmp_path, table_path, overwrite = TRUE)
      unlink(tmp_path)
    }
  }
}

#' Acquire/release a simple, dependency-free cross-process file lock
#'
#' `dir.create()` is atomic on every filesystem this pipeline runs on: only
#' one caller ever succeeds in creating a given directory, so it doubles as
#' a mutex without needing a locking package. Used to serialize read-modify
#' -write access to one shared covariate-forecast table when multiple
#' parallel workers (e.g. `evaluate_all_units_covariate_forecasts()`'s
#' cluster) might otherwise race to append to it at the same time.
#'
#' IMPORTANT: only ever wrap the table's own read/write I/O in this lock,
#' never a compute step. An earlier version of `get_or_compute_root_forecast()`
#' held this lock for its *entire* body, including the covariate model fit
#' itself -- fine with the synthetic test data used during development
#' (fits took milliseconds), but on the real pipeline's much longer
#' histories, a single slow fit (tbats()/stlf() could genuinely take tens
#' of seconds) monopolized the lock and serialized every other worker
#' behind it; anything still waiting past `timeout` errored out instead of
#' ever getting a turn. That's what turned "run this for every unit" into
#' "1-2 out of 104 tasks complete, the rest time out". Keeping this lock's
#' critical section to just the read/write below (compute happens outside
#' it, see `get_or_compute_root_forecast()`) means lock hold time is
#' bounded by CSV I/O, not by how long a fit takes.
.with_covariate_forecast_table_lock <- function(table_path, expr_fn, timeout = 60) {
  lock_path <- paste0(table_path, ".lock")
  start <- Sys.time()
  while (!dir.create(lock_path, showWarnings = FALSE, recursive = TRUE)) {
    if (as.numeric(Sys.time() - start, units = "secs") > timeout) {
      stop("Could not acquire lock on ", table_path, " within ", timeout, "s ",
           "(stale lock? remove ", lock_path, " manually if no process is using it).")
    }
    Sys.sleep(0.1)
  }
  on.exit(unlink(lock_path, recursive = TRUE), add = TRUE)
  expr_fn()
}

#' Look up one unit's root-series forecast in the persistent table, or
#' compute it and persist it if it isn't there yet
#'
#' This is the single place both the on-demand path (inside
#' `forecast_covariates()`, when a formula needs a root series mid-fit) and
#' the bulk precompute driver (`build_covariate_forecast_table.r`) go
#' through, so a table built once by one and filled in on-demand by the
#' other are always in the exact same format.
#'
#' @param y            The root series' history (already sliced to
#'                      train_start..train_end, in epiweek order).
#' @param epiweeks,h,n_years  As in `forecast_covariate_series()`.
#' @param table_path   Path from `covariate_forecast_table_path()`.
#' @param id           This unit's id (UF code or municipality geocode) --
#'                      part of the lookup key, since the table spans every
#'                      unit of a disease/level.
#' @param root_var     The root series' column name (part of the lookup key).
#' @param train_start,train_end  Integer YYYYWW fold bounds -- part of the
#'                      lookup key, so a fold window shifting (e.g. the
#'                      submission window advancing as the season
#'                      progresses) is naturally treated as new work rather
#'                      than silently reusing a forecast for a different
#'                      window. NOTE: if historical weather/case data for
#'                      already-cached epiweeks is later corrected upstream,
#'                      that alone does NOT invalidate a cached row (the key
#'                      doesn't fingerprint the historical values) -- delete
#'                      the affected row(s), or the whole table, to force a
#'                      refresh after a real data correction.
#'
#' @return Numeric vector of length `h`.
get_or_compute_root_forecast <- function(y, epiweeks, h, table_path, id, root_var,
                                          train_start, train_end, n_years = 5) {
  key_match <- function(tbl) {
    tbl$id          == as.character(id) &
      tbl$root_var    == root_var &
      tbl$train_start == as.character(train_start) &
      tbl$train_end   == as.character(train_end) &
      tbl$h           == as.character(h) &
      tbl$n_years     == as.character(n_years)
  }

  # Lock-free read first: a cache hit needs no lock at all, and -- per the
  # note on .with_covariate_forecast_table_lock() -- nothing that can take
  # a nontrivial amount of time should ever happen while the lock is held.
  # Worst case on a miss is a small amount of duplicated compute below,
  # resolved by the re-check inside the write's own (brief) lock.
  tbl <- .read_covariate_forecast_table(table_path)
  hit <- tbl[key_match(tbl), , drop = FALSE]
  if (nrow(hit) == h) {
    return(as.numeric(hit$forecast_value)[order(as.integer(hit$step))])
  }

  fc <- forecast_covariate_series(y = y, epiweeks = epiweeks, h = h, n_years = n_years)

  .with_covariate_forecast_table_lock(table_path, function() {
    tbl <- .read_covariate_forecast_table(table_path)
    hit <- tbl[key_match(tbl), , drop = FALSE]
    if (nrow(hit) == h) {
      return(invisible(NULL))  # another worker already wrote this exact key
    }

    new_rows <- data.frame(
      id             = as.character(id),
      root_var       = root_var,
      train_start    = as.character(train_start),
      train_end      = as.character(train_end),
      h              = as.character(h),
      n_years        = as.character(n_years),
      step           = as.character(seq_len(h)),
      forecast_value = as.character(fc),
      stringsAsFactors = FALSE
    )
    # Drop any stale partial rows for this exact key before appending (e.g.
    # a previous run was interrupted mid-write) so re-runs stay consistent.
    stale <- key_match(tbl)
    tbl <- tbl[!stale, , drop = FALSE]
    tbl <- rbind(tbl, new_rows)

    .write_covariate_forecast_table(tbl, table_path)
    invisible(NULL)
  })

  fc
}

#' Forecast every regressor column a fit_sarimax_epiweek() call needs
#'
#' A naive implementation would fit `forecast_covariate_series()`
#' independently to every column in `rhs_terms`. That's correct for genuinely
#' independent series (PCA components, ENSO/IOD/PDO, a base weather
#' aggregate requested on its own) -- but state-level data also carries
#' lag/rolling-mean columns that are deterministic transforms of a single
#' base series (see `parse_derived_covariate()`), and forecasting those
#' independently is not just 9x more model fits than necessary, it's wrong:
#' `<base>_lag4`'s true value at forecast step h is exactly `<base>`'s own
#' value at step h-4 (already-known history for h<=4), so modeling it
#' separately throws away exact information and can produce a trajectory
#' that contradicts the base series' own forecast trend. This function
#' instead forecasts each distinct base series exactly once, appends that
#' forecast to the real historical tail, and re-derives every lag/rolling
#' column from that one combined series with the same `lag()`/trailing-mean
#' arithmetic `agg_data_uf.r` used to build them from real data -- so a
#' `<base>_lag4` forecast for h<=4 is the literal historical value, and for
#' h>4 is mechanically tied to `<base>`'s own forecast rather than an
#' independently-fit guess.
#'
#' NOTE: this does not (yet) reach into PCA components -- PC1/PC2 are still
#' forecast independently of each other even though both are fixed linear
#' combinations of the same underlying weather history. See the comment in
#' `sarimax/src/utils.r` above `pca_all()`/`run_grid_search()` if that gap
#' is ever worth closing (it would mean threading the fitted PCA rotation
#' through to here and reprojecting forecasted raw covariates, rather than
#' forecasting PC scores directly).
#'
#' @param data         The full state/disease data frame passed to
#'                      `fit_sarimax_epiweek()` (must contain `epiweek` and
#'                      every column named in `rhs_terms`, plus -- for any
#'                      derived `rhs_terms` entry -- that entry's base
#'                      column, even if the base itself isn't in
#'                      `rhs_terms`).
#' @param rhs_terms     Character vector of column names to forecast -- the
#'                      SARIMAX formula's regressors (e.g. `c("PC1", "PC2",
#'                      "iod")`, or, at state level, raw names like
#'                      `c("temp_med_mean_lag4", "temp_med_mean_mean_6mo")`).
#' @param train_start  Integer YYYYWW: first epiweek used to fit each
#'                      covariate's own model (inclusive) -- pass the same
#'                      value used for the disease model's own training
#'                      window, so the covariate forecasts respect the same
#'                      cutoff.
#' @param train_end    Integer YYYYWW: last epiweek used to fit each
#'                      covariate's own model (inclusive).
#' @param h            Integer: number of steps (weeks) ahead to forecast --
#'                      pass `n_fc` (the length of the target epiweek range).
#' @param epiweek_col  Name of the epiweek column (default "epiweek").
#' @param n_years      Passed through to `forecast_covariate_series()`
#'                      (default 5): how many recent years' worth of each
#'                      calendar week to average per forecast step.
#' @param table_path   Optional: path to a persistent covariate-forecast
#'                      table (see `covariate_forecast_table_path()`). When
#'                      supplied (together with `id`), each root series
#'                      needed is looked up there first -- and only fit
#'                      live (then saved back for next time) on a miss --
#'                      instead of always re-fitting. Default `NULL`
#'                      preserves the original always-fit-live behavior
#'                      exactly.
#' @param id           Required when `table_path` is supplied: this unit's
#'                      id (UF code or municipality geocode), since the
#'                      table spans every unit of a disease/level.
#'
#' @return A numeric matrix with `h` rows and `length(rhs_terms)` columns,
#'   column names equal to `rhs_terms` and in the same order, suitable for
#'   `scale()`-ing and passing to `forecast::forecast(..., xreg = )`.
forecast_covariates <- function(data, rhs_terms, train_start, train_end, h,
                                 epiweek_col = "epiweek", n_years = 5,
                                 table_path = NULL, id = NULL) {
  if (!is.null(table_path) && is.null(id)) {
    stop("forecast_covariates(): 'id' is required when 'table_path' is supplied.")
  }
  missing_cols <- setdiff(c(epiweek_col, rhs_terms), names(data))
  if (length(missing_cols) > 0) {
    stop("forecast_covariates(): data is missing column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  hist_rows <- data[[epiweek_col]] >= train_start & data[[epiweek_col]] <= train_end
  if (sum(hist_rows) == 0) {
    stop("forecast_covariates(): no rows found for epiweeks ",
         train_start, "-", train_end)
  }
  # Explicit epiweek-order sort (rather than trusting caller order) so the
  # lag/rolling reconstruction below always shifts/averages along the true
  # time axis.
  hist_data <- data[hist_rows, , drop = FALSE]
  hist_data <- hist_data[order(hist_data[[epiweek_col]]), , drop = FALSE]

  # ── Classify requested terms and find every base series that needs its
  #    own model fit: one per distinct base, whether requested directly or
  #    only needed to reconstruct a derived (lag/rolling) term. ────────────
  derived_spec <- setNames(lapply(rhs_terms, parse_derived_covariate), rhs_terms)
  is_derived   <- !vapply(derived_spec, is.null, logical(1))

  base_needed <- unique(c(
    rhs_terms[!is_derived],
    vapply(derived_spec[is_derived], `[[`, character(1), "base")
  ))
  missing_bases <- setdiff(base_needed, names(hist_data))
  if (length(missing_bases) > 0) {
    stop("forecast_covariates(): derived covariate(s) requested but their ",
         "base column(s) are missing from data: ", paste(missing_bases, collapse = ", "))
  }

  base_forecasts <- setNames(
    lapply(base_needed, function(b) {
      if (is.null(table_path)) {
        forecast_covariate_series(y = hist_data[[b]], epiweeks = hist_data[[epiweek_col]],
                                   h = h, n_years = n_years)
      } else {
        get_or_compute_root_forecast(
          y = hist_data[[b]], epiweeks = hist_data[[epiweek_col]], h = h, n_years = n_years,
          table_path = table_path, id = id, root_var = b,
          train_start = train_start, train_end = train_end
        )
      }
    }),
    base_needed
  )

  out <- matrix(NA_real_, nrow = h, ncol = length(rhs_terms),
                dimnames = list(NULL, rhs_terms))

  for (term in rhs_terms) {
    spec <- derived_spec[[term]]

    if (is.null(spec)) {
      out[, term] <- base_forecasts[[term]]
      next
    }

    base_hist <- hist_data[[spec$base]]
    combined  <- c(base_hist, base_forecasts[[spec$base]])
    H         <- length(base_hist)

    if (spec$kind == "lag") {
      if (H < spec$n) {
        stop(sprintf(
          "forecast_covariates(): only %d week(s) of training history available, ",
          H), sprintf("not enough to reconstruct '%s' (needs %d week(s) of lookback).",
                       term, spec$n))
      }
      # <base>_lag<n> at absolute position t is base(t - n); positions
      # H+1..H+h are the forecast horizon within `combined`.
      out[, term] <- combined[(H + 1):(H + h) - spec$n]
    } else { # "roll"
      if (H < spec$k) {
        stop(sprintf(
          "forecast_covariates(): only %d week(s) of training history available, ",
          H), sprintf("not enough to reconstruct '%s' (needs a %d-week trailing window).",
                       term, spec$k))
      }
      # <base>_mean_Xmo at absolute position t is mean(base((t-k)..(t-1)))
      # -- a k-week trailing window ending the week before t (see
      # `agg_data_uf.r`'s `dplyr::lag(mean_run(., k), 1)`).
      out[, term] <- vapply((H + 1):(H + h), function(t) {
        mean(combined[(t - spec$k):(t - 1)])
      }, numeric(1))
    }
  }

  out
}
