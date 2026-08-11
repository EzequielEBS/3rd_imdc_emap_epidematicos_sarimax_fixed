#' Submit forecasts to the Mosqlimate platform
#'
#' Reads the best forecasts produced by `fit.r` (one CSV per state x target
#' window under `sarimax/results/preds/`) and uploads them to the
#' Mosqlimate Predictions Registry through the `mosqlient` Python package
#' (accessed from R via `reticulate`).
#'
#' Credentials: the Mosqlimate API key is NEVER hardcoded here. It is read
#' from the `MOSQLIMATE_API_KEY` environment variable, which should be
#' defined in a local, git-ignored `.env` file as:
#'   MOSQLIMATE_API_KEY=your-key-here
#' Get your key from the "Auth" section of your Mosqlimate profile
#' (https://mosqlimate.org/). Never commit a real key to a public repository.

library(reticulate)
library(data.table)
library(dotenv)
library(readr)
library(dplyr)

py_config()
py_require(c("epiweeks", "python-dotenv", "mosqlient"))

# Load MOSQLIMATE_API_KEY from a local .env file (ignored by git, see .gitignore).
# Falls back to whatever is already in the environment if no .env is present.
if (file.exists(".env")) dotenv::load_dot_env(".env")

required_prediction_columns <- c(
  "date",
  "lower_95",
  "lower_90",
  "lower_80",
  "lower_50",
  "pred",
  "upper_50",
  "upper_80",
  "upper_90",
  "upper_95"
)

# Official IMDC prediction windows. Used only to validate file completeness
# immediately before upload.
prediction_set_specs <- tibble::tribble(
  ~prediction_set, ~start_epiweek, ~end_epiweek,
  "target_1", 202241L, 202340L,
  "target_2", 202341L, 202440L,
  "target_3", 202441L, 202540L,
  "target_4", 202541L, 202640L,
  "final",    202641L, 202740L
)

epiweek_to_date_submission <- function(yw) {
  yr   <- yw %/% 100L
  wk   <- yw  %% 100L
  jan4 <- as.Date(paste0(yr, "-01-04"))
  sunday_w1 <- jan4 - as.integer(format(jan4, "%w"))
  sunday_w1 + (wk - 1L) * 7L
}

has_53_weeks_submission <- function(year) {
  epiweek_to_date_submission(year * 100L + 53L) <
    epiweek_to_date_submission((year + 1L) * 100L + 1L)
}

enumerate_epiweeks_submission <- function(start, end) {
  out <- integer(0)
  yr <- start %/% 100L
  wk <- start %% 100L
  end_year <- end %/% 100L
  end_week <- end %% 100L

  repeat {
    out <- c(out, yr * 100L + wk)
    if (yr == end_year && wk == end_week) break
    max_week <- if (has_53_weeks_submission(yr)) 53L else 52L
    if (wk < max_week) {
      wk <- wk + 1L
    } else {
      yr <- yr + 1L
      wk <- 1L
    }
  }
  out
}

#' Validate a prediction file against the IMDC weekly submission rules.
validate_imdc_prediction <- function(pred, prediction_set = NULL) {
  missing_cols <- setdiff(required_prediction_columns, names(pred))
  if (length(missing_cols) > 0) {
    stop(
      "Prediction is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  pred <- pred |>
    select(all_of(required_prediction_columns)) |>
    mutate(date = as.Date(date)) |>
    arrange(date)

  if (any(is.na(pred$date)))
    stop("Prediction contains invalid dates.")

  # Sunday = 0 under %w.
  if (any(as.integer(format(pred$date, "%w")) != 0))
    stop("All prediction dates must be Sundays.")

  if (nrow(pred) > 1 && any(diff(pred$date) != 7))
    stop("Prediction dates must be continuous weekly dates with no gaps.")

  numeric_cols <- setdiff(required_prediction_columns, "date")
  if (any(!is.finite(as.matrix(pred[, numeric_cols]))))
    stop("Prediction contains NA/NaN/Inf values.")

  if (any(as.matrix(pred[, numeric_cols]) < 0))
    stop("All predictions and interval bounds must be non-negative.")

  nested_ok <- with(
    pred,
    lower_95 <= lower_90 &
      lower_90 <= lower_80 &
      lower_80 <= lower_50 &
      lower_50 <= pred &
      pred <= upper_50 &
      upper_50 <= upper_80 &
      upper_80 <= upper_90 &
      upper_90 <= upper_95
  )

  if (!all(nested_ok))
    stop("Prediction intervals are not nested around the median.")

  if (!is.null(prediction_set)) {
    spec <- prediction_set_specs |>
      filter(.data$prediction_set == !!prediction_set)

    if (nrow(spec) != 1)
      stop("Unknown prediction set: ", prediction_set)

    expected_epiweeks <- enumerate_epiweeks_submission(
      spec$start_epiweek,
      spec$end_epiweek
    )
    expected_dates <- epiweek_to_date_submission(expected_epiweeks)

    if (!identical(pred$date, expected_dates)) {
      stop(
        prediction_set,
        " does not contain the complete official EW41-EW40 weekly window."
      )
    }
  }

  # Keep the API payload explicit and portable across reticulate versions.
  pred |>
    mutate(date = format(date, "%Y-%m-%d"))
}

#' Submit a batch of state-level forecasts for a single disease to Mosqlimate
#'
#' @param disease Character. ICD-10 code submitted to Mosqlimate
#'   (e.g. "A90" for dengue, "A92.0" for chikungunya).
#' @param disease_label Character. Lowercase label used in file paths
#'   (e.g. "dengue", "chikungunya"). Defaults to `disease` if not given.
#' @param commit Character. Commit hash of the model version used to
#'   generate these predictions.
#' @param repository Character. GitHub repository in "owner/repo" form.
#' @param states Character vector of two-letter state (UF) codes to submit.
#'   Defaults to all 26 states + DF.
#' @param n_targets Integer. Number of target windows per state
#'   (target_1..target_n). Defaults to 4.
#' @param processed_data_dir Character. Base directory containing the
#'   aggregated case-count CSVs, expected at
#'   `<processed_data_dir>/<disease_label>/<disease_label>_<state>_agg.csv.gz`.
#' @param preds_dir Character. Directory containing prediction CSVs produced
#'   by `fit.r`, expected at
#'   `<preds_dir>/pred_<disease_label>_<state>_target_<i>.csv`.
#' @param case_definition Character. Case definition reported to Mosqlimate
#'   (e.g. "probable").
#' @param adm_level Integer. Administrative level of the prediction
#'   (1 = state/UF).
#' @param adm_2 Optional. Sub-state admin code, NULL for state-level.
#' @param published Logical. Whether the prediction should be published.
#' @param api_key Character. Mosqlimate API key. Defaults to the
#'   `MOSQLIMATE_API_KEY` environment variable.
#' @param on_duplicate Character. What to do when Mosqlimate reports that a
#'   prediction already exists for this model/commit/state/target:
#'   "skip" (default, log and move on), "error" (stop the whole run), or
#'   "overwrite" (delete the existing prediction via `mosq$remove_prediction()`
#'   if available, then re-upload). "overwrite" requires the mosqlient
#'   version installed to expose a removal function; if it doesn't, this
#'   falls back to "skip" with a warning.
#' @param verbose Logical. Print progress messages. Defaults to TRUE.
#'
#' @return Invisibly, a data.table logging each state/target pair and its
#'   status ("submitted", "duplicate_skipped", "duplicate_overwritten",
#'   "missing_case_file", "missing_pred_file", or "error").
submit_mosqlimate_forecasts <- function(
  disease,
  commit,
  repository,
  disease_label,
  unit_type = c("state", "city"),
  units,
  prediction_sets = paste0("target_", 1:4),
  processed_data_dir = "processed_data",
  preds_dir = "sarimax/results/preds",
  case_definition = "probable",
  published = TRUE,
  api_key = Sys.getenv("MOSQLIMATE_API_KEY"),
  on_duplicate = c("skip", "error"),
  verbose = TRUE
) {

  unit_type <- match.arg(unit_type)
  on_duplicate <- match.arg(on_duplicate)

  if (identical(api_key, "")) {
    stop(
      "MOSQLIMATE_API_KEY is not set. Create a .env file in the repository ",
      "root with a line `MOSQLIMATE_API_KEY=your-key-here` (see sarimax/README.md), ",
      "or pass `api_key` explicitly."
    )
  }

  if (!grepl("^[0-9a-fA-F]{40}$", commit)) {
    stop("`commit` must be the full 40-character commit hash used for this run.")
  }

  mosq <- import("mosqlient")
  log_rows <- list()

  for (unit in units) {
    if (unit_type == "state") {
      case_file <- file.path(
        processed_data_dir,
        disease_label,
        paste0(disease_label, "_", unit, "_agg.csv.gz")
      )
      adm_level <- 1L
    } else {
      case_file <- file.path(
        processed_data_dir,
        disease_label,
        "sel_cities",
        paste0(disease_label, "_", unit, ".csv.gz")
      )
      adm_level <- 2L
    }

    if (!file.exists(case_file)) {
      warning("Skipping ", unit, ": processed data file not found: ", case_file)
      next
    }

    d <- read_csv(case_file, show_col_types = FALSE)
    adm_1 <- as.integer(d$uf_code[1])
    adm_2 <- if (unit_type == "city") as.integer(unit) else NULL

    for (prediction_set in prediction_sets) {
      pred_file <- file.path(
        preds_dir,
        paste0(
          "pred_", disease_label, "_", unit, "_",
          prediction_set, ".csv"
        )
      )

      if (!file.exists(pred_file)) {
        warning(
          "Skipping ", unit, " ", prediction_set,
          ": prediction file not found: ", pred_file
        )
        next
      }

      pred <- read_csv(pred_file, show_col_types = FALSE)
      pred <- validate_imdc_prediction(pred, prediction_set = prediction_set)

      description <- paste0(
        tools::toTitleCase(disease_label),
        " | ",
        if (unit_type == "state") "state " else "city ",
        unit,
        " | ",
        prediction_set,
        " | 3rd IMDC 2026"
      )

      if (verbose) {
        message(
          "Submitting ",
          disease_label, " | ", unit_type, " ", unit,
          " | ", prediction_set
        )
      }

      status <- tryCatch({
        mosq$upload_prediction(
          api_key = api_key,
          repository = repository,
          description = description,
          commit = commit,
          disease = disease,
          case_definition = case_definition,
          adm_level = adm_level,
          adm_0 = "BRA",
          adm_1 = adm_1,
          adm_2 = adm_2,
          published = published,
          prediction = pred
        )
        "submitted"
      }, error = function(e) {
        msg <- conditionMessage(e)
        is_duplicate <- grepl("Duplication found", msg, fixed = TRUE)

        if (is_duplicate && on_duplicate == "skip") {
          if (verbose) message("  -> duplicate, skipping")
          return("duplicate_skipped")
        }

        if (verbose) message("  -> error: ", msg)
        if (on_duplicate == "error") stop(e)
        "error"
      })

      log_rows[[length(log_rows) + 1]] <- data.table(
        disease = disease_label,
        unit_type = unit_type, unit = as.character(unit),
        prediction_set = prediction_set, commit = commit, status = status
      )
    }
  }

  if (length(log_rows) == 0)
    return(invisible(data.table()))

  result <- rbindlist(log_rows)

  if (verbose) {
    message("\nSummary:")
    print(result[, .N, by = status])
  }

  invisible(result)
}

# ── Example usage ─────────────────────────────────────────────────────────
# Fill in the commit hashes and repository name for the run being submitted,
# then call the function once per disease.

repository <- "EzequielEBS/3rd_imdc_emap_epidematicos_sarimax_state/fixed"

states <- c(
  "AC", "AL", "AM", "AP", "BA", "CE", "DF", "GO", "MA",
  "MG", "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ",
  "RN", "RO", "RR", "RS", "SC", "SE", "SP", "TO"
)

cities_dengue <- c(
  2931350, 2933307, 2302503, 3119401, 3549805,
  3541406, 1200401, 1200203, 1716109, 4113700,
  4103701, 4104808, 5201405, 5102637, 5215231
)

cities_chikungunya <- c(
  2211001, 2931350, 3143302, 3119401, 1721000,
  1716109, 4104808, 4219507, 5103403, 5102637
)

# Dengue state: four validations + final target.
submit_mosqlimate_forecasts(
  disease       = "A90",
  disease_label = "dengue",
  commit        = "372294033bf4aa2586eb4d770ef75acf7dec21a3",
  repository    = repository,
  unit_type     = "state",
  units         = states,
  prediction_sets = c(paste0("target_", 1:4), "final")
)

# Dengue cities: four validations + final target.
#submit_mosqlimate_forecasts(
#  disease       = "A90",
#  disease_label = "dengue",
#  commit        = ??,
#  repository    = repository,
#  unit_type     = "city",
#  units         = cities_dengue,
#  prediction_sets = c(paste0("target_", 1:4), "final")
#)

# Chikungunya state: current instructions define the four validations.
submit_mosqlimate_forecasts(
  disease       = "A92.0",
  disease_label = "chikungunya",
  commit        = "ec240514fbcdffa0bfdaadcdc65c056e6a80aff8",
  repository    = repository,
  unit_type     = "state",
  units         = states,
  prediction_sets = paste0("target_", 1:4)
)

# Chikungunya cities: current instructions define the four validations.
#submit_mosqlimate_forecasts(
#  disease       = "A92.0",
#  disease_label = "chikungunya",
#  commit        = ??,
#  repository    = repository,
#  unit_type     = "city",
#  units         = cities_chikungunya,
#  prediction_sets = paste0("target_", 1:4)
#)