source("sarimax/src/utils.r")

# Load support tables once; they are reused by every state/city run.
forecasting_climate_all <- read_csv(
  "processed_data/climate/forecasting_climate.csv.gz",
  show_col_types = FALSE
)

geocode_map_all <- read_csv(
  "raw_data/dengue.csv.gz",
  show_col_types = FALSE
) |>
  distinct(geocode, uf, uf_code)
  
#' Run full covariate/order model selection for one state or city
#'
#' End-to-end pipeline for a single state (or city) and disease: loads the
#' aggregated covariate file, filters candidate covariates (low-variance,
#' weak-correlation), optionally reduces them via fold-consistent PCA
#' (`pca_all()`) and tests whether ENSO/IOD/PDO add signal beyond the
#' retained PCs (`filter_redundant_indices()`), builds a set of candidate
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
  n_cores = max(1L, parallel::detectCores() - 1L),
  method = "CSS-ML",
  lambda = NULL,
  pca = T,
  pca_var_threshold = 0.90,
  k = 5,
  index_cor_threshold = 0.3,  # min |partial cor| for enso/iod/pdo to be kept (see filter_redundant_indices)
  fixed_stat_par = F,         # search d/D via CV WIS instead of locking to a KPSS/OCSB pick
  bootstrap = TRUE,           # simulate forecast paths for better-calibrated intervals
  npaths = 1000,
  climatology_years = 5,
  smooth_weeks = 2,
  concluded_states_path = "sarimax/results/concluded_states_dengue.csv",  # overridable so test
                                                                   # runs don't touch the
                                                                   # real checkpoint file
  concluded_cities_path = "sarimax/results/concluded_cities_dengue.csv"  # overridable so test
                                                                   # runs don't touch the
                                                                   # real checkpoint file
) {
  if (is.null(state) && is.null(city))
    stop("Provide either `state` or `city`.")
  if (!is.null(state) && !is.null(city))
    stop("Provide only one of `state` or `city`.")
  
  if (!is.null(state)) {
    concluded_states <- if (is.null(concluded_states)) {
      data.frame(state = character(0), stringsAsFactors = FALSE)
    } else concluded_states

    if (state %in% concluded_states$state) {
      message("State ", state, " already processed. Skipping.")
      return(NULL)
    }

    message("Processing state: ", state)
    file_name <- paste0(
      "processed_data/", disease, "/",
      disease, "_", state, "_agg.csv.gz"
    )
    unit_type <- "state"
    unit_id <- state
  } else {
    concluded_cities <- if (is.null(concluded_cities)) {
      data.frame(city = character(0), stringsAsFactors = FALSE)
    } else concluded_cities

    if (city %in% concluded_cities$city) {
      message("City ", city, " already processed. Skipping.")
      return(NULL)
    }

    message("Processing city: ", city)
    file_name <- paste0(
      "processed_data/", disease, "/sel_cities/",
      disease, "_", city, ".csv.gz"
    )
    unit_type <- "city"
    unit_id <- city
  }

    master_data <- read_csv(file_name, show_col_types = FALSE)

  model_data <- prepare_imdc_model_data(
    data = master_data,
    forecasting_climate = forecasting_climate_all,
    geocode_map = geocode_map_all,
    unit_type = unit_type,
    unit_id = unit_id,
    pca = pca,
    pca_var_threshold = pca_var_threshold,
    k = k,
    threshold_low_variance = threshold_low_variance,
    threshold_cor = threshold_cor,
    min_cor = min_cor,
    max_size_covariates = max_size_covariates,
    index_cor_threshold = index_cor_threshold,
    climatology_years = climatology_years,
    smooth_weeks = smooth_weeks
  )

  formulas <- model_data$formulas

  if (length(formulas) == 0)
    formulas <- list(reformulate(character(0), response = "cases"))

  if (pca) {
    sample_formulas <- formulas
  } else {
    set.seed(123)
    sample_formulas <- formulas[
      sample(length(formulas), min(sample_size, length(formulas)))
    ]
  }

  # Screen ARIMA orders using the same leakage-free split construction used
  # later for the final model comparison.
  order_screen <- run_grid_search(
    splits   = model_data$splits,
    formulas = sample_formulas,
    levels   = levels,
    max_order = max_order,
    n_cores   = n_cores,
    method    = method,
    lambda    = lambda,
    fixed_stat_par = fixed_stat_par,
    bootstrap = bootstrap,
    npaths    = npaths
  )

  summarise_grid <- function(x) {
    x |>
      group_by(formula_id, order) |>
      filter(all(!failed)) |>
      summarise(
        mean_metric = mean(.data[[metric]], na.rm = TRUE),
        .groups = "drop"
      ) |>
      arrange(mean_metric)
  }

  order_metrics <- summarise_grid(order_screen$metrics)
  if (nrow(order_metrics) == 0)
    stop("No formula/order combination succeeded in all four backtests.")
    
  best_order <- order_metrics$order[1]
  ord <- parse_order(best_order)

  if (pca) {
    metrics <- order_metrics
    best_formula <- metrics$formula_id[1]
    selected_results <- order_screen
  } else {
    formula_results <- run_grid_search(
      splits = model_data$splits,
      formulas = formulas,
      levels = levels,
      fixed_order = c(
        p = ord$order[1], d = ord$order[2], q = ord$order[3],
        P = ord$seasonal_order[1],
        D = ord$seasonal_order[2],
        Q = ord$seasonal_order[3]
      ),
      n_cores = n_cores,
      method = method,
      lambda = lambda,
      bootstrap = bootstrap,
      npaths = npaths
    )

    metrics <- summarise_grid(formula_results$metrics)
    if (nrow(metrics) == 0)
      stop("No covariate formula succeeded in all four backtests.")

    best_formula <- metrics$formula_id[1]
    selected_results <- formula_results
  }

  # Persist the preprocessing/model settings together with the leaderboard so
  # fit.r can reproduce exactly the recipe that was validated by WIS.
  metrics <- metrics |>
    mutate(
      pca = pca,
      pca_var_threshold = pca_var_threshold,
      k = k,
      threshold_low_variance = threshold_low_variance,
      threshold_cor = threshold_cor,
      min_cor = min_cor,
      max_size_covariates = max_size_covariates,
      index_cor_threshold = index_cor_threshold,
      climatology_years = climatology_years,
      smooth_weeks = smooth_weeks,
      method = method,
      lambda = if (is.null(lambda)) NA_real_ else as.numeric(lambda),
      bootstrap = bootstrap,
      npaths = npaths
    )

  best_pred <- selected_results$predictions |>
    filter(
      formula_id == best_formula,
      order == best_order
    )

  final_metrics <- selected_results$metrics |>
    filter(
      formula_id == best_formula,
      order == best_order
    ) |>
    arrange(target_id) |>
    mutate(
      pca = pca,
      pca_var_threshold = pca_var_threshold,
      k = k,
      climatology_years = climatology_years,
      smooth_weeks = smooth_weeks
    )

  if (!dir.exists("sarimax/results/metrics"))
    dir.create("sarimax/results/metrics", recursive = TRUE)

  # save files — use whichever identifier (state or city) this run was for,
  # so per-city runs don't all collide on the same output file.
  id <- if (!is.null(state)) state else city
  metric_name <- paste0("mean_", metric)

  write_csv(final_metrics, file.path("sarimax/results/metrics/", paste0(metric_name, "_best_model_", disease, "_", id, ".csv")))
  write_csv(metrics, file.path("sarimax/results/metrics/", paste0("metrics_all_formulas_", disease, "_", id, ".csv")))

  # update concluded states
  if (!is.null(state)) {
    concluded_states <- rbind(concluded_states, data.frame(state = state, stringsAsFactors = FALSE))
    write_csv(concluded_states, concluded_states_path)
  } else {
    concluded_cities <- rbind(concluded_cities, data.frame(city = city, stringsAsFactors = FALSE))
    write_csv(concluded_cities, concluded_cities_path)
  }

  return(list(
    best_order = best_order,
    best_formula = best_formula,
    best_pred = best_pred,
    final_metrics = final_metrics,
    all_metrics = metrics,
    recipe = model_data
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