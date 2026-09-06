source("sarimax/src/utils.r")

#' Build a `cases ~ ...` formula, tolerating an empty term list
#'
#' `reformulate()` errors ("'termlabels' must be a character vector of length
#' at least one") when `vars` is empty instead of producing a valid
#' intercept-only formula. That happens for real: a state/city can have every
#' local-weather covariate AND every climate index filtered out, leaving no
#' terms at all. Falls back to `cases ~ 1` in that case so `run_model_selection()`
#' still fits a plain (covariate-free) SARIMA instead of crashing.
safe_reformulate <- function(vars, response = "cases") {
  if (length(vars) == 0) {
    as.formula(paste(response, "~ 1"))
  } else {
    reformulate(vars, response = response)
  }
}

#' Run full covariate/order model selection for one state or city
#'
#' End-to-end pipeline for a single state (or city) and disease: loads the
#' aggregated covariate file, filters candidate covariates (low-variance,
#' weak-correlation), then either (a) collapses each base local-weather
#' variable's lag/rolling family down to its single best-correlated variant
#' via `select_best_per_variable()` and builds covariate combinations
#' (`pca = FALSE`), or (b) applies that same `select_best_per_variable()`
#' collapse before reducing the resulting one-column-per-variable set via
#' fold-consistent PCA (`pca_all()`) and testing whether ENSO/IOD/PDO add
#' signal beyond the retained PCs (`filter_redundant_indices()`)
#' (`pca = TRUE`, the default). Both branches apply the same collapse step
#' first so a base variable's own lag/rolling transforms -- which are fixed
#' transforms of one series, not genuinely distinct series (see
#' `src/forecast_covariates.r`) -- are never fed into either covariate
#' combinations or the PCA input matrix as if they were independent
#' variables; only the single most-predictive variant per base variable
#' reaches that stage. Then, either way, this builds a set of candidate
#' formulas, grid-searches SARIMAX (p,d,q)(P,D,Q) orders via
#' `run_grid_search()` ranked by mean cross-validated `metric` (default WIS),
#' and writes the best model\'s per-split metrics and the full formula/order
#' leaderboard to `sarimax/results/metrics/`.
#'
#' @param state                   State abbreviation to process (mutually
#'                                exclusive with `city`); skipped if already
#'                                in `concluded_states$state`.
#' @param concluded_states        Data frame of already-processed states
#'                                (column `state`), used to skip reruns.
#' @param city                    City code to process (mutually exclusive
#'                                with `state`); skipped if already in
#'                                `concluded_cities$city`.
#' @param concluded_cities        Data frame of already-processed cities
#'                                (column `city`), used to skip reruns.
#' @param disease                 "dengue" or "chikungunya" (default "dengue");
#'                                selects the input file path.
#' @param metric                  Name of the metric column (from
#'                                `compute_metrics()`) used to rank candidate
#'                                formulas/orders (default "wis").
#' @param sample_size             When `pca = FALSE`, number of candidate
#'                                formulas randomly sampled for the order
#'                                screen (default 10).
#' @param threshold_low_variance  Minimum covariate standard deviation kept
#'                                by `filter_low_variance()` (default 0.01).
#' @param threshold_cor           Pairwise correlation threshold used by
#'                                `build_covariate_combinations()` when
#'                                `pca = FALSE` (default 0.6).
#' @param min_cor                 Minimum |correlation| with response kept by
#'                                `filter_by_correlation()` (default 0.1).
#' @param max_size_covariates     Maximum covariates per formula when
#'                                `pca = FALSE` (default 3).
#' @param levels                  Prediction-interval coverage levels passed
#'                                through to `run_grid_search()` (default
#'                                c(50, 80, 90, 95)).
#' @param max_order               Named vector of upper bounds for the
#'                                (p,d,q)(P,D,Q) grid search (default
#'                                c(p=2, d=1, q=2, P=1, D=0, Q=1)).
#' @param n_cores                 Parallel workers for `run_grid_search()`
#'                                (default detectCores() - 1).
#' @param method                  Estimation method passed to `Arima()`
#'                                (default "CSS-ML").
#' @param lambda                  Box-Cox lambda; NULL (default) uses
#'                                log1p/expm1 transform instead.
#' @param pca                     Logical: reduce covariates via PCA instead
#'                                of exhaustive combination search (default
#'                                TRUE).
#' @param pca_var_threshold       Cumulative variance threshold for PCA
#'                                component selection (default 0.90).
#' @param k                       Maximum number of principal components to
#'                                build formulas from, when `pca = TRUE`
#'                                (default 5).
#' @param index_cor_threshold     Minimum partial correlation for ENSO/IOD/PDO
#'                                to be kept alongside PCs (default 0.3; see
#'                                `filter_redundant_indices()`).
#' @param fixed_stat_par          Logical: lock differencing orders to a
#'                                single KPSS/OCSB-recommended value instead
#'                                of searching d/D via CV WIS (default FALSE).
#' @param bootstrap               Logical: use bootstrapped simulation paths
#'                                for prediction intervals (default TRUE).
#' @param npaths                  Number of bootstrap simulation paths
#'                                (default 1000).
#' @param concluded_states_path   Output path for the updated concluded-states
#'                                checkpoint file.
#' @param concluded_cities_path   Output path for the updated concluded-cities
#'                                checkpoint file.
#'
#' @return A list with `best_order`, `best_formula`, `best_pred`,
#'         `final_metrics` (per train/target split), and `all_metrics`
#'         (full formula/order leaderboard). Also writes CSV results to
#'         `sarimax/results/metrics/` and updates the concluded-states/cities
#'         checkpoint file as a side effect.
run_model_selection <- function(
  state = NULL,
  concluded_states = NULL,
  city = NULL,
  concluded_cities = NULL,
  disease = "dengue",
  metric = "wis",
  sample_size = 10,
  threshold_low_variance = 0.01,
  threshold_cor = 0.6,
  min_cor = 0.1,
  max_size_covariates = 3,
  levels = c(50, 80, 90, 95),
  max_order = c(p = 2, d = 1, q = 2, P = 1, D = 0, Q = 1),
  n_cores = parallel::detectCores() - 1,
  method = "CSS-ML",
  lambda = NULL,
  pca = T,
  pca_var_threshold = 0.90,
  k = 5,
  index_cor_threshold = 0.3,  # min |partial cor| for enso/iod/pdo to be kept (see filter_redundant_indices)
  fixed_stat_par = F,         # search d/D via CV WIS instead of locking to a KPSS/OCSB pick
  bootstrap = TRUE,           # simulate forecast paths for better-calibrated intervals
  npaths = 1000,
  concluded_states_path = "sarimax/results/concluded_states_dengue.csv",  # overridable so test
                                                                   # runs don't touch the
                                                                   # real checkpoint file
  concluded_cities_path = "sarimax/results/concluded_cities_dengue.csv"  # overridable so test
                                                                   # runs don't touch the
                                                                   # real checkpoint file
) {
  if (!is.null(state) | !is.null(city)) {
    if (!is.null(state)) {
      if (state %in% concluded_states$state) {
        message("State ", state, " already processed. Skipping.")
        return(NULL)
      }
      message("Processing state: ", state)
      level <- "state"; id <- state
    } else if (!is.null(city)) {
      if (city %in% concluded_cities$city) {
        message("City ", city, " already processed. Skipping.")
        return(NULL)
      }
      message("Processing city: ", city)
      level <- "city"; id <- city
    }
  }

  # load_unit_data() (sarimax/src/utils.r) normalizes the state-vs-city
  # column-name differences (municipality files use "casos" and lack the
  # "_mean" suffix on weather variables) so everything below -- get_candidates(),
  # pca_all(), fit_sarimax() via run_grid_search(), reformulate(response="cases") --
  # sees one consistent schema regardless of level.
  dengue_state <- load_unit_data(disease = disease, level = level, id = id)

  train_ids  <- paste0("train_", 1:4)
  
  candidates <- get_candidates(dengue_state)
  # filter_low_variance()/filter_by_correlation() each message() every column
  # they drop, which -- with up to 9 lag/rolling columns per base variable --
  # is a lot of noise that doesn't say what actually matters: which
  # covariates end up driving the model. Silenced here; what's actually used
  # for PCA (or, when pca = FALSE, for the covariate-combination search) is
  # printed explicitly below instead, once it's known.
  candidates <- suppressMessages(filter_low_variance(
    dengue_state[
      dengue_state[[train_ids[1]]] == 1,
    ],
    candidates,
    threshold = threshold_low_variance
  ))
  candidates <- suppressMessages(filter_by_correlation(
    dengue_state[
      dengue_state[[train_ids[1]]] == 1,
    ],
    candidates,
    min_cor = min_cor
  ))

  if (pca) {
    non_index_candidates <- candidates[!grepl("enso|iod|pdo", candidates)]
    climate_indices <- intersect(c("enso", "iod", "pdo"), candidates)

    # Collapse each base local-weather variable's lag/rolling family down to
    # its single best-correlated-with-response variant BEFORE handing
    # anything to PCA -- exactly the same collapse the `pca = FALSE` branch
    # below already applies before build_covariate_combinations(). Without
    # this, a variable that still has several lag/rolling columns surviving
    # filtering (all fixed transforms of one series, not genuinely distinct
    # series -- see src/forecast_covariates.r) would contribute that many
    # near-duplicate, mutually correlated columns to the PCA input matrix,
    # letting whichever variable happened to retain the most transforms
    # dominate the leading components rather than each distinct physical
    # driver contributing roughly one dimension.
    if (length(non_index_candidates) > 0) {
      non_index_candidates <- select_best_per_variable(
        dengue_state[
          dengue_state[[train_ids[1]]] == 1,
        ],
        non_index_candidates
      )
    }

    if (length(non_index_candidates) == 0) {
      # Variance/correlation filtering dropped every local-weather covariate
      # for this state/city, leaving a zero-column matrix. prcomp()/svd()
      # can't factorize that ("a dimension is zero"), so skip PCA entirely
      # and fall back to climate indices alone (if any survived the
      # correlation filter) or a plain covariate-free SARIMA.
      message(
        "No local-weather covariates survived filtering for ",
        if (!is.null(state)) state else city,
        " — skipping PCA; fitting with climate indices only (if any) or no covariates."
      )
      pcs <- character(0)
      max_k <- 0
      # No PCs to test redundancy against, so keep whatever indices already
      # passed filter_by_correlation() above as-is.
      kept_indices <- climate_indices
      if (length(kept_indices) == 0) {
        message(
          "No climate indices survived filtering either for ",
          if (!is.null(state)) state else city,
          " — fitting a covariate-free SARIMA."
        )
      }
      formulas <- list(safe_reformulate(kept_indices))
      candidates <- kept_indices
    } else {
      message(sprintf(
        "Applying PCA to %d covariate(s): %s",
        length(non_index_candidates), paste(non_index_candidates, collapse = ", ")
      ))
      pca_result <- pca_all(
        data = dengue_state,
        candidates = non_index_candidates,
        var_threshold = pca_var_threshold
      )
      dengue_state <- pca_result$data
      pcs <- pca_result$variables
      # Only go up to the components selected by var_threshold
      max_k <- min(k, length(pcs))

      # ENSO/IOD/PDO were deliberately excluded from the PCA matrix above so they
      # aren't blended away by an unsupervised projection. Test whether they carry
      # signal beyond what the retained local-weather PCs already explain, and add
      # back any that do, instead of dropping them outright.
      kept_indices <- character(0)
      if (length(climate_indices) > 0) {
        retained <- filter_redundant_indices(
          data       = dengue_state,
          covariates = c(pcs, climate_indices),
          indices    = climate_indices,
          threshold  = index_cor_threshold
        )
        kept_indices <- intersect(climate_indices, retained)
      }

      formulas <- lapply(seq_len(max_k), function(i) {
        safe_reformulate(pcs[1:i])
      })
      if (length(kept_indices) > 0) {
        formulas <- c(formulas, lapply(seq_len(max_k), function(i) {
          safe_reformulate(c(pcs[1:i], kept_indices))
        }))
      }
      candidates <- c(pcs, kept_indices)
    }
  } else {
    candidates <- select_best_per_variable(
      dengue_state[
        dengue_state[[train_ids[1]]] == 1,
      ],
      candidates
    )
    if (length(candidates) == 0) {
      # build_covariate_combinations() can't factor a zero-column correlation
      # matrix (errors "'x' is empty" out of cor()) -- same "everything got
      # filtered out" situation the pca = TRUE branch handles above, just
      # reached via the exhaustive-combination path instead of PCA.
      message(
        "No covariates survived filtering for ",
        if (!is.null(state)) state else city,
        " — fitting a covariate-free SARIMA."
      )
      formulas <- list(safe_reformulate(character(0)))
    } else {
      message(sprintf(
        "Building covariate combinations from %d candidate(s): %s",
        length(candidates), paste(candidates, collapse = ", ")
      ))
      combos <- build_covariate_combinations(
        data       = dengue_state,
        covariates = candidates,
        threshold  = threshold_cor,
        max_size   = max_size_covariates
      )
      formulas   <- lapply(combos, \(vars) safe_reformulate(vars))
    }
  }
  
  if (pca) {
    sample_formulas <- formulas
  } else {
    sample_formulas <- formulas[sample(length(formulas), min(sample_size, length(formulas)))]
  }
  
  order_screen <- run_grid_search(
    data      = dengue_state,
    formulas  = sample_formulas,
    levels    = levels,
    max_order = max_order,
    n_cores   = n_cores,
    method    = method,
    lambda    = lambda,
    fixed_stat_par = fixed_stat_par,
    bootstrap = bootstrap,
    npaths    = npaths,
    # Share this unit's covariate forecasts with fit.r/evaluate_covariate_forecasts.r
    # (and a later re-run of this same function) via the persistent table --
    # see covariate_forecast_table_path().
    table_path = covariate_forecast_table_path(disease, level),
    id         = id
  )

  # If EVERY formula/order combination failed to fit, group_by()...
  # filter(all(!failed)) below removes every group, and the metric column
  # (e.g. "wis") never got created at all since compute_metrics() was never
  # reached for any of them. That used to surface three calls later as a
  # cryptic "subscript out of bounds" from parse_order() on an empty
  # best_order -- worse, `get(metric)` can silently resolve to an unrelated
  # same-named function from an attached package (scoringutils exports a
  # `wis()` function) instead of erroring "object not found", so the
  # underlying all-failed condition was easy to miss entirely. Catch it here
  # instead, right where the actual cause -- and run_grid_search()'s
  # "Failure reasons:" summary above -- are still in view.
  if (nrow(order_screen$metrics) == 0 || all(order_screen$metrics$failed)) {
    stop(sprintf(
      "run_model_selection(): every formula/order combination failed to fit for %s -- see the 'Failure reasons:' summary printed above by run_grid_search().",
      if (!is.null(state)) state else city
    ))
  }

  metric_name <- paste0("mean_", metric)
  best_par <- order_screen$metrics |>
    group_by(formula_id, order) |>
    filter(all(!failed)) |>
    summarise(mean_metric = mean(.data[[metric]]), .groups = "drop") |>
    arrange(mean_metric) |>
    slice(1)
  best_order <- best_par$order
  ord <- parse_order(best_order)

  if (pca) {
    best_formula <- best_par$formula_id
    metrics <- order_screen$metrics |>
      group_by(formula_id, order) |>
      filter(all(!failed)) |>
      summarise(mean_metric = mean(.data[[metric]]), .groups = "drop") |>
      arrange(mean_metric)
    best_pred <- order_screen$predictions |>
      filter(formula_id == best_formula, order == best_order)
  } else {
    formula_results <- run_grid_search(
      data      = dengue_state,
      formulas  = formulas,
      levels    = levels,
      fixed_order = c(p = ord$order[1], d = ord$order[2], q = ord$order[3],
                      P = ord$seasonal_order[1], D = ord$seasonal_order[2], Q = ord$seasonal_order[3]),
      n_cores   = n_cores,
      method    = method,
      lambda    = lambda,
      bootstrap = bootstrap,
      npaths    = npaths,
      table_path = covariate_forecast_table_path(disease, level),
      id         = id
    )
    if (nrow(formula_results$metrics) == 0 || all(formula_results$metrics$failed)) {
      stop(sprintf(
        "run_model_selection(): every formula failed to fit at the chosen order for %s -- see the 'Failure reasons:' summary printed above by run_grid_search().",
        if (!is.null(state)) state else city
      ))
    }
    metrics <- formula_results$metrics |>
      group_by(formula_id, order) |>
      filter(all(!failed)) |>        # keep only groups where ALL splits succeeded
      summarise(mean_metric = mean(.data[[metric]]), .groups = "drop") |>
      arrange(mean_metric)

    best_formula <- metrics$formula_id[1]
    best_pred <- formula_results$predictions |>
      filter(formula_id == best_formula, order == best_order)
  }

  train_ids <- paste0("train_", 1:4)
  target_ids <- paste0("target_", 1:4)

  final_metrics <- lapply(seq_along(train_ids), function(i) {
    actual    <- dengue_state$cases[dengue_state[[target_ids[i]]] == 1]
    pred     <- best_pred |> filter(train_id == train_ids[i])
    pred_df  <- pred |> select(date, pred, starts_with("lower_"), starts_with("upper_"))
    compute_metrics(pred, actual)
  }) |> bind_rows() |> mutate(target_id = target_ids)
  # save files — use whichever identifier (state or city) this run was for,
  # so per-city runs don't all collide on the same output file.
  id <- if (!is.null(state)) state else city
  write_csv(final_metrics, file.path("sarimax/results/metrics/", paste0(metric_name, "_best_model_", disease, "_", id, ".csv")))
  write_csv(metrics, file.path("sarimax/results/metrics/", paste0("metrics_all_formulas_", disease, "_", id, ".csv")))

  # update concluded states
  if (!is.null(state)) {
    concluded_states <- rbind(concluded_states, data.frame(state = state, stringsAsFactors = FALSE))
    write_csv(concluded_states, concluded_states_path)
  } else if (!is.null(city)) {
    concluded_cities <- rbind(concluded_cities, data.frame(city = city, stringsAsFactors = FALSE))
    write_csv(concluded_cities, concluded_cities_path)
  }
  return(list(
    best_order = best_order,
    best_formula = best_formula,
    best_pred = best_pred,
    final_metrics = final_metrics,
    all_metrics = metrics
  ))
}


states <- c("AC", "AL", "AM", "AP", "BA", "CE", "DF", "GO", "MA",
            "MG", "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ",
            "RN", "RO", "RR", "RS", "SC", "SE", "SP", "TO")


#-----------------------------------------------------------------------------------
# Run model selection for each state and save results
#-----------------------------------------------------------------------------------

results_state_dengue <- lapply(states, function(st) {
  concluded_states <- tryCatch({
    read_csv("sarimax/results/concluded_states_dengue.csv", show_col_types = FALSE)
  }, error = function(e) {
    data.frame(state = character(0), stringsAsFactors = FALSE)
  })
  run_model_selection(
    state = st,
    concluded_states = concluded_states
  )
}
)

# get summary
best_wis_df_dengue <- lapply(states, function(st) {
  file_name <- paste0("sarimax/results/metrics/metrics_all_formulas_dengue_", st, ".csv")
  data <- read_csv(file_name, show_col_types = FALSE)
  data[1, ] |> mutate(state = st, .before = 1)
}) |> bind_rows()
write_csv(best_wis_df_dengue, "sarimax/results/metrics/best_wis_dengue_all_states.csv")

results_state_chikungunya <- lapply(states, function(st) {
  concluded_states <- tryCatch({
    read_csv("sarimax/results/concluded_states_chikungunya.csv", show_col_types = FALSE)
  }, error = function(e) {
    data.frame(state = character(0), stringsAsFactors = FALSE)
  })
  run_model_selection(
    state = st,
    concluded_states = concluded_states,
    disease = "chikungunya",
    concluded_states_path = "sarimax/results/concluded_states_chikungunya.csv"
  )
}
)

best_wis_df_chikungunya <- lapply(states, function(st) {
  file_name <- paste0("sarimax/results/metrics/metrics_all_formulas_chikungunya_", st, ".csv")
  data <- read_csv(file_name, show_col_types = FALSE)
  data[1, ] |> mutate(state = st, .before = 1)
}) |> bind_rows()
write_csv(best_wis_df_chikungunya, "sarimax/results/metrics/best_wis_chikungunya_all_states.csv")

#-----------------------------------------------------------------------------------
# Run model selection for each city and save results
#-----------------------------------------------------------------------------------

cities_dengue <- c(
  2931350,
  2933307,
  2302503,
  3119401,
  3549805,
  3541406,
  1200401,
  1200203,
  1716109,
  4113700,
  4103701,
  4104808,
  5201405,
  5102637,
  5215231
)

cities_chikungunya <- c(
  2211001,
  2931350,
  3143302,
  3119401,
  1721000,
  1716109,
  4104808,
  4219507,
  5103403,
  5102637
)

results_city_dengue <- lapply(cities_dengue, function(city) {
  concluded_cities <- tryCatch({
    read_csv("sarimax/results/concluded_cities_dengue.csv", show_col_types = FALSE)
  }, error = function(e) {
    data.frame(city = character(0), stringsAsFactors = FALSE)
  })
  run_model_selection(
    city = city,
    concluded_cities = concluded_cities
  )
}
)

best_wis_df_dengue_cities <- lapply(cities_dengue, function(city) {
  file_name <- paste0("sarimax/results/metrics/metrics_all_formulas_dengue_", city, ".csv")
  data <- read_csv(file_name, show_col_types = FALSE)
  data[1, ] |> mutate(city = city, .before = 1)
}) |> bind_rows()
write_csv(best_wis_df_dengue_cities, "sarimax/results/metrics/best_wis_dengue_all_cities.csv")

results_city_chikungunya <- lapply(cities_chikungunya, function(city) {
  concluded_cities <- tryCatch({
    read_csv("sarimax/results/concluded_cities_chikungunya.csv", show_col_types = FALSE)
  }, error = function(e) {
    data.frame(city = character(0), stringsAsFactors = FALSE)
  })
  run_model_selection(
    city = city,
    concluded_cities = concluded_cities,
    disease = "chikungunya",
    concluded_cities_path = "sarimax/results/concluded_cities_chikungunya.csv"
  )
}
)

best_wis_df_chikungunya_cities <- lapply(cities_chikungunya, function(city) {
  file_name <- paste0("sarimax/results/metrics/metrics_all_formulas_chikungunya_", city, ".csv")
  data <- read_csv(file_name, show_col_types = FALSE)
  data[1, ] |> mutate(city = city, .before = 1)
}) |> bind_rows()
write_csv(best_wis_df_chikungunya_cities, "sarimax/results/metrics/best_wis_chikungunya_all_cities.csv")