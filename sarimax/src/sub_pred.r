#' Submit forecasts to the Mosqlimate platform
#'
#' Reads the best forecasts produced by `fit.r` (one CSV per unit x target
#' window under `sarimax/results/preds/`) and uploads them to the
#' Mosqlimate Predictions Registry through the `mosqlient` Python package
#' (accessed from R via `reticulate`). Works at both geographic levels this
#' pipeline produces forecasts for -- state (`adm_level = 1`, IBGE UF code)
#' and municipality (`adm_level = 2`, IBGE municipality geocode) -- since
#' `fit.r` already writes both under the same `pred_<disease>_<id>_<target_id>.csv`
#' naming convention (`id` a UF code or a geocode); this file previously only
#' wired up the state-level path even though `adm_2` already existed as a
#' parameter for exactly this purpose.
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

py_config()
py_require(c("epiweeks", "python-dotenv", "mosqlient"))

# Load MOSQLIMATE_API_KEY from a local .env file (ignored by git, see .gitignore).
# Falls back to whatever is already in the environment if no .env is present.
if (file.exists(".env")) dotenv::load_dot_env(".env")

#' Resolve the case-count file path for one unit, mirroring `load_unit_data()`'s
#' two naming conventions (`sarimax/src/utils.r`) without pulling in that
#' file's much heavier dependency chain (forecast_covariates.r, scoringutils,
#' etc.) -- sub_pred.r is meant to be runnable on its own, right after fit.r,
#' with nothing else sourced.
#'
#' @param level Either "state" or "city".
resolve_case_file <- function(level, disease_label, id, processed_data_dir) {
  if (level == "state") {
    file.path(processed_data_dir, disease_label, paste0(disease_label, "_", id, "_agg.csv.gz"))
  } else {
    file.path(processed_data_dir, disease_label, "sel_cities", paste0(disease_label, "_", id, ".csv.gz"))
  }
}

#' Resolve a target spec into the target_<n>/target_forecast split(s) to
#' submit
#'
#' Duplicated (not sourced) from fit.r's copy -- same reasoning as
#' resolve_case_file() above: this file is meant to be runnable on its own.
#' Keep the two copies in sync if this changes.
#'
#' @param target  NULL (-> all four validation splits, no label), a numeric
#'   subset of c(1,2,3,4) (no label), or one of "validation" (-> c(1,2,3,4),
#'   no label -- submitting all four doesn't call for one shared label),
#'   "validation_1".."validation_4" (-> that one split, label
#'   "Validation Test <n>"), "forecast" (-> the separate Forecast Phase
#'   submission, its own "forecast" slot -- never 4 -- label
#'   "final forecast").
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

#' Submit a batch of forecasts for a single disease/level to Mosqlimate
#'
#' @param disease Character. ICD-10 code submitted to Mosqlimate
#'   (e.g. "A90" for dengue, "A92.0" for chikungunya).
#' @param disease_label Character. Lowercase label used in file paths
#'   (e.g. "dengue", "chikungunya"). Defaults to `disease` if not given.
#' @param commit Character. Commit hash of the model version used to
#'   generate these predictions.
#' @param repository Character. GitHub repository in "owner/repo" form.
#' @param level Character. "state" (adm_level 1, UF-level forecasts) or
#'   "city" (adm_level 2, municipality-level forecasts). Defaults to
#'   "state".
#' @param ids Character vector of unit identifiers to submit -- two-letter
#'   UF codes when `level = "state"`, IBGE municipality geocodes when
#'   `level = "city"`. Defaults to all 26 states + DF for `level = "state"`;
#'   for `level = "city"`, defaults to every municipality listed in
#'   `<metrics_dir>/best_wis_<disease_label>_all_cities.csv` (the same
#'   leaderboard `fit.r`'s `run_all_units()` reads), erroring if that file
#'   doesn't exist yet and `ids` wasn't supplied explicitly.
#' @param target Which window(s) to submit per unit. Defaults to NULL (all
#'   four validation splits). Either a numeric subset of c(1,2,3,4), or a
#'   friendlier string spec -- "validation" (all four), "validation_1"..
#'   "validation_4" (just that one split), or "forecast" (the separate
#'   Forecast Phase submission) -- resolved via `resolve_target_spec()`
#'   above. "forecast" reads `pred_<disease_label>_<id>_target_forecast.csv`
#'   -- a different file from `..._target_4.csv` (Validation Test 4) -- so
#'   submitting one never touches or duplicates the other.
#' @param target_label Optional character. Cosmetic only: if supplied,
#'   overrides the label a string `target` spec would otherwise supply on
#'   its own (e.g. "forecast" already implies "final forecast";
#'   "validation_4" already implies "Validation Test 4") -- used in each
#'   upload's `description` in place of "target_<n>" so Mosqlimate
#'   submission descriptions say what the target actually represents
#'   instead of just its slot number. Doesn't affect which file is read or
#'   where it's uploaded. Only meaningful when `target` selects a single
#'   split.
#' @param processed_data_dir Character. Base directory containing the
#'   aggregated case-count CSVs. State files are expected at
#'   `<processed_data_dir>/<disease_label>/<disease_label>_<id>_agg.csv.gz`;
#'   municipality files at
#'   `<processed_data_dir>/<disease_label>/sel_cities/<disease_label>_<id>.csv.gz`.
#' @param preds_dir Character. Directory containing prediction CSVs produced
#'   by `fit.r`, expected at
#'   `<preds_dir>/pred_<disease_label>_<id>_<target_id>.csv`.
#' @param metrics_dir Character. Directory holding `best_wis_..._all_cities.csv`,
#'   used to auto-discover `ids` when `level = "city"` and `ids` isn't given.
#' @param case_definition Character. Case definition reported to Mosqlimate
#'   (e.g. "probable").
#' @param published Logical. Whether the prediction should be published.
#' @param api_key Character. Mosqlimate API key. Defaults to the
#'   `MOSQLIMATE_API_KEY` environment variable.
#' @param on_duplicate Character. What to do when Mosqlimate reports that a
#'   prediction already exists for this model/commit/unit/target:
#'   "skip" (default, log and move on), "error" (stop the whole run), or
#'   "overwrite" (delete the existing prediction via `mosq$remove_prediction()`
#'   if available, then re-upload). "overwrite" requires the mosqlient
#'   version installed to expose a removal function; if it doesn't, this
#'   falls back to "skip" with a warning.
#' @param verbose Logical. Print progress messages. Defaults to TRUE.
#'
#' @return Invisibly, a data.table logging each unit/target pair and its
#'   status ("submitted", "duplicate_skipped", "duplicate_overwritten",
#'   "missing_case_file", "missing_pred_file", or "error").
submit_mosqlimate_forecasts <- function(disease,
                                         commit,
                                         repository,
                                         disease_label = disease,
                                         level = c("state", "city"),
                                         ids = NULL,
                                         target = NULL,
                                         target_label = NULL,
                                         processed_data_dir = "processed_data",
                                         preds_dir = "sarimax/results/preds",
                                         metrics_dir = "sarimax/results/metrics",
                                         case_definition = "probable",
                                         published = TRUE,
                                         api_key = Sys.getenv("MOSQLIMATE_API_KEY"),
                                         on_duplicate = c("skip", "error", "overwrite"),
                                         verbose = TRUE) {

  level        <- match.arg(level)
  on_duplicate <- match.arg(on_duplicate)
  adm_level    <- if (level == "state") 1L else 2L

  spec    <- resolve_target_spec(target)
  targets <- spec$targets
  if (is.null(target_label)) target_label <- spec$label

  if (!is.null(target_label) && length(targets) > 1) {
    stop(
      "target_label is only meaningful with a single target -- got ",
      length(targets), " (", paste(targets, collapse = ", "), "). ",
      "Call once per target if you need a distinct label for each."
    )
  }

  if (identical(api_key, "")) {
    stop(
      "MOSQLIMATE_API_KEY is not set. Create a .env file in the repository ",
      "root with a line `MOSQLIMATE_API_KEY=your-key-here` (see sarimax/README.md), ",
      "or pass `api_key` explicitly."
    )
  }

  if (is.null(ids)) {
    if (level == "state") {
      ids <- c("AC", "AL", "AM", "AP", "BA", "CE", "DF", "GO", "MA",
               "MG", "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ",
               "RN", "RO", "RR", "RS", "SC", "SE", "SP", "TO")
    } else {
      best_wis_path <- file.path(metrics_dir, paste0("best_wis_", disease_label, "_all_cities.csv"))
      if (!file.exists(best_wis_path)) {
        stop(
          "`ids` not supplied and ", best_wis_path, " not found -- run ",
          "model_sel.r (and then fit.r) for city-level ", disease_label,
          " first, or pass `ids` explicitly."
        )
      }
      ids <- as.character(read_csv(best_wis_path, show_col_types = FALSE)$city)
    }
  }

  mosq <- import("mosqlient")
  target_ids <- paste0("target_", sort(unique(targets)))

  log_rows <- list()

  for (id in ids) {

    case_file <- resolve_case_file(level, disease_label, id, processed_data_dir)
    if (!file.exists(case_file)) {
      warning("Skipping ", id, ": case-count file not found at ", case_file)
      next
    }
    d <- read_csv(case_file, show_col_types = FALSE)
    adm_1 <- d$uf_code[1]
    adm_2 <- if (level == "state") NULL else id

    for (target_id in target_ids) {
      pred_file <- file.path(
        preds_dir,
        paste0("pred_", disease_label, "_", id, "_", target_id, ".csv")
      )
      if (!file.exists(pred_file)) {
        warning("Skipping ", id, " ", target_id, ": prediction file not found at ", pred_file)
        next
      }
      pred <- read_csv(pred_file, show_col_types = FALSE)
      target_desc <- if (!is.null(target_label)) target_label else target_id
      description <- paste0(
        tools::toTitleCase(disease_label), " prediction for ",
        if (level == "state") paste0("state ", id) else paste0("municipality ", id),
        " and target ", target_desc
      )

      if (verbose) message("Submitting ", disease_label, " | ", level, " | ", id, " | ", target_id)

      do_upload <- function() {
        mosq$upload_prediction(
          api_key = api_key,
          repository = repository,
          description = description,
          commit = commit,
          disease = disease,
          case_definition = case_definition,
          adm_level = adm_level,
          adm_1 = adm_1,
          adm_2 = adm_2,
          published = published,
          prediction = pred
        )
      }

      status <- tryCatch({
        do_upload()
        "submitted"
      }, error = function(e) {
        msg <- conditionMessage(e)
        is_duplicate <- grepl("Duplication found", msg, fixed = TRUE)

        if (!is_duplicate) {
          if (verbose) message("  -> error: ", msg)
          return("error")
        }

        if (on_duplicate == "error") {
          stop(e)
        } else if (on_duplicate == "skip") {
          if (verbose) message("  -> duplicate, skipping")
          return("duplicate_skipped")
        } else { # overwrite
          if (verbose) message("  -> duplicate, attempting overwrite")
          removed <- tryCatch({
            if ("remove_prediction" %in% names(mosq)) {
              mosq$remove_prediction(
                api_key = api_key,
                repository = repository,
                commit = commit,
                disease = disease,
                adm_level = adm_level,
                adm_1 = adm_1,
                adm_2 = adm_2
              )
              TRUE
            } else {
              FALSE
            }
          }, error = function(e2) {
            if (verbose) message("  -> removal failed: ", conditionMessage(e2))
            FALSE
          })

          if (!removed) {
            warning(
              "Could not remove existing prediction for ", id, " ", target_id,
              " (mosqlient has no removal function, or removal failed). ",
              "Delete it manually on mosqlimate.org if you need to overwrite. Skipping."
            )
            return("duplicate_skipped")
          }

          re_status <- tryCatch({
            do_upload()
            "duplicate_overwritten"
          }, error = function(e3) {
            if (verbose) message("  -> re-upload failed: ", conditionMessage(e3))
            "error"
          })
          return(re_status)
        }
      })

      log_rows[[length(log_rows) + 1]] <- data.table(
        disease = disease_label, level = level, id = id, target = target_id,
        commit = commit, status = status
      )
    }
  }

  result <- rbindlist(log_rows)
  if (verbose && nrow(result) > 0) {
    message("\nSummary:")
    print(result[, .N, by = status])
  }

  invisible(result)
}

# ── Example usage ─────────────────────────────────────────────────────────
# Fill in the commit hashes and repository name for the run being submitted,
# then call the function once per (disease, level) combination you want to
# submit. State-level is the official IMDC submission level; municipality-
# level submission is available here too (adm_level = 2) but check the
# current challenge rules before including it in an official entry -- see
# sarimax/README.md's note on the focal municipalities being, historically,
# a supplementary/exploratory level rather than a submitted one.

repository <- "EzequielEBS/3rd_imdc_emap_epidematicos_sarimax_fixed"

# Dengue -- states
submit_mosqlimate_forecasts(
  disease       = "A90",
  disease_label = "dengue",
  commit        = "45830c1e76b9645aef7a011e4cf2f7d0fe6c02dd",
  repository    = repository,
  level         = "state",
  on_duplicate  = "skip",
  target        = "validation"
)

# Dengue -- municipalities
submit_mosqlimate_forecasts(
  disease       = "A90",
  disease_label = "dengue",
  commit        = "372294033bf4aa2586eb4d770ef75acf7dec21a3",
  repository    = repository,
  level         = "city",
  on_duplicate  = "skip",
  target        = "validation"
)

# Chikungunya -- states
submit_mosqlimate_forecasts(
  disease       = "A92.0",
  disease_label = "chikungunya",
  commit        = "ec240514fbcdffa0bfdaadcdc65c056e6a80aff8",
  repository    = repository,
  level         = "state",
  on_duplicate  = "skip",
  target        = "validation"
)

# Chikungunya -- municipalities
submit_mosqlimate_forecasts(
  disease       = "A92.0",
  disease_label = "chikungunya",
  commit        = "ec240514fbcdffa0bfdaadcdc65c056e6a80aff8",
  repository    = repository,
  level         = "city",
  on_duplicate  = "skip",
  target        = "validation"
)

# ── Submitting the Forecast Phase ────────────────────────────────────────────
# The Forecast Phase (challenge rules, Section 5.2) is a separate submission
# from the four validation splits above -- a different season/window,
# produced by fit.r's own separate `TARGET_TO_FIT <- "forecast"` run, and
# read from its own `..._target_forecast.csv` file (never `..._target_4.csv`,
# which stays Validation Test 4 and is untouched by this). Pass
# target = "forecast" to submit just that:
#
submit_mosqlimate_forecasts(
  disease       = "A90",
  disease_label = "dengue",
  commit        = "45830c1e76b9645aef7a011e4cf2f7d0fe6c02dd",
  repository    = repository,
  level         = "state",
  target        = "forecast",
  on_duplicate  = "skip"
)
