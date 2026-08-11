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

if (!dir.exists("sarimax/results/preds")) {
  dir.create("sarimax/results/preds", recursive = TRUE)
}

#-----------------------------------------------------------------------------------
# State-level dengue: four validation forecasts + final 2026-2027 forecast
#-----------------------------------------------------------------------------------

best_wis_df_dengue_state <- read_csv("sarimax/results/metrics/best_wis_dengue_all_states.csv", show_col_types = FALSE)

preds_dengue_state <- lapply(best_wis_df_dengue_state$state, function(st) {
  file_name <- paste0("processed_data/dengue/dengue_", st, "_agg.csv.gz")
  d <- read_csv(file_name, show_col_types = FALSE)

  best_row <- best_wis_df_dengue_state |>
    filter(state == st) |>
    slice(1)

  model_data <- prepare_imdc_model_data(
    data = d,
    forecasting_climate = forecasting_climate_all,
    geocode_map = geocode_map_all,
    unit_type = "state",
    unit_id = st,
    pca = best_row$pca[1],
    pca_var_threshold = best_row$pca_var_threshold[1],
    k = best_row$k[1],
    threshold_low_variance = best_row$threshold_low_variance[1],
    threshold_cor = best_row$threshold_cor[1],
    min_cor = best_row$min_cor[1],
    max_size_covariates = best_row$max_size_covariates[1],
    index_cor_threshold = best_row$index_cor_threshold[1],
    climatology_years = best_row$climatology_years[1],
    smooth_weeks = best_row$smooth_weeks[1]
  )

  best_order <- best_row |> pull(order)
  ord <- parse_order(best_order)
  best_formula <- best_row |> pull(formula_id)
  if (is.na(best_formula) || identical(best_formula, "")) {
    best_formula <- character(0)
  }

  lambda_value <- if (is.na(best_row$lambda[1])) NULL else best_row$lambda[1]

  train_ids  <- paste0("train_", 1:4)
  target_ids <- paste0("target_", 1:4)

  pred_target <- lapply(seq_along(train_ids), function(i) {
    train_id <- train_ids[i]
    target_id <- target_ids[i]
    split <- model_data$splits[[i]]

    fit <- fit_sarimax(
      data = split$data,
      formula = best_formula,
      order = c(ord$order[1], ord$order[2], ord$order[3]),
      seasonal = list(order = c(ord$seasonal_order[1], ord$seasonal_order[2], ord$seasonal_order[3]), period = 52),
      train_id = NULL,
      method = best_row$method[1],
      lambda = lambda_value,
      bootstrap = best_row$bootstrap[1],
      npaths = best_row$npaths[1],
      seed = 123 + i
    )

    write_csv(fit, file.path("sarimax/results/preds/", paste0("pred_dengue_", st, "_", target_id, ".csv")))
    fit
  })

  names(pred_target) <- target_ids

  # Final challenge forecast: cutoff EW25/2026, target EW41/2026-EW40/2027.
  final_split <- prepare_imdc_final_data(
    data = d,
    recipe = model_data,
    forecasting_climate = forecasting_climate_all,
    geocode_map = geocode_map_all,
    unit_type = "state",
    unit_id = st
  )

  final_fit <- fit_sarimax(
    data = final_split$data,
    formula = best_formula,
    order = c(ord$order[1], ord$order[2], ord$order[3]),
    seasonal = list(order = c(ord$seasonal_order[1], ord$seasonal_order[2], ord$seasonal_order[3]), period = 52),
    train_id = NULL,
    method = best_row$method[1],
    lambda = lambda_value,
    bootstrap = best_row$bootstrap[1],
    npaths = best_row$npaths[1],
    seed = 128
  )

  write_csv(final_fit, file.path("sarimax/results/preds/", paste0("pred_dengue_", st, "_final.csv")))
  pred_target[["final"]] <- final_fit

  bind_rows(pred_target, .id = "target_id") |>
    mutate(state = st)
})

#-----------------------------------------------------------------------------------
# State-level chikungunya: four validation forecasts
#-----------------------------------------------------------------------------------

best_wis_df_chikungunya_state <- read_csv("sarimax/results/metrics/best_wis_chikungunya_all_states.csv", show_col_types = FALSE)
state_warnings_chikungunya <- list()

preds_chikungunya_state <- lapply(best_wis_df_chikungunya_state$state, function(st) {
  file_name <- paste0("processed_data/chikungunya/chikungunya_", st, "_agg.csv.gz")
  d <- read_csv(file_name, show_col_types = FALSE)

  best_row <- best_wis_df_chikungunya_state |>
    filter(state == st) |>
    slice(1)

  model_data <- prepare_imdc_model_data(
    data = d,
    forecasting_climate = forecasting_climate_all,
    geocode_map = geocode_map_all,
    unit_type = "state",
    unit_id = st,
    pca = best_row$pca[1],
    pca_var_threshold = best_row$pca_var_threshold[1],
    k = best_row$k[1],
    threshold_low_variance = best_row$threshold_low_variance[1],
    threshold_cor = best_row$threshold_cor[1],
    min_cor = best_row$min_cor[1],
    max_size_covariates = best_row$max_size_covariates[1],
    index_cor_threshold = best_row$index_cor_threshold[1],
    climatology_years = best_row$climatology_years[1],
    smooth_weeks = best_row$smooth_weeks[1]
  )

  best_order <- best_row |> pull(order)
  ord <- parse_order(best_order)
  best_formula <- best_row |> pull(formula_id)
  if (is.na(best_formula) || identical(best_formula, "")) {
    best_formula <- character(0)
  }

  lambda_value <- if (is.na(best_row$lambda[1])) NULL else best_row$lambda[1]

  train_ids  <- paste0("train_", 1:4)
  target_ids <- paste0("target_", 1:4)

  pred_target <- lapply(seq_along(train_ids), function(i) {
    train_id <- train_ids[i]
    target_id <- target_ids[i]
    split <- model_data$splits[[i]]

    fit <- fit_sarimax(
      data = split$data,
      formula = best_formula,
      order = c(ord$order[1], ord$order[2], ord$order[3]),
      seasonal = list(order = c(ord$seasonal_order[1], ord$seasonal_order[2], ord$seasonal_order[3]), period = 52),
      train_id = NULL,
      method = best_row$method[1],
      lambda = lambda_value,
      bootstrap = best_row$bootstrap[1],
      npaths = best_row$npaths[1],
      seed = 123 + i
    )

    # ── Report warnings with state/target context ──────────────────────────
    w <- attr(fit, "warnings")
    if (!is.null(w)) {
      state_warnings_chikungunya[[st]] <<- c(state_warnings_chikungunya[[st]],
        setNames(w, rep(target_id, length(w))))
      message(sprintf("[%s | %s] %d warning(s):\n%s",
                      st, target_id, length(w),
                      paste0("  - ", w, collapse = "\n")))
    }
    write_csv(fit, file.path("sarimax/results/preds/", paste0("pred_chikungunya_", st, "_", target_id, ".csv")))
    fit
  })
  names(pred_target) <- target_ids
  bind_rows(pred_target, .id = "target_id") |>
    mutate(state = st)
})

#-----------------------------------------------------------------------------------
# Retry chikungunya states that produced actionable warnings
#-----------------------------------------------------------------------------------

states_to_retry_chikungunya <- names(state_warnings_chikungunya)
# states_to_retry_chikungunya <- c("DF", "MS", "RO", "SP")
state_warnings_chikungunya <- list()
preds_chikungunya_retry    <- list()
resolved_formulas_chikungunya <- tibble(state = character(), formula_id = character(), order = character(), row = integer())

for (st in states_to_retry_chikungunya) {
  file_name <- paste0("processed_data/chikungunya/chikungunya_", st, "_agg.csv.gz")
  d <- read_csv(file_name, show_col_types = FALSE)

  metrics_file <- paste0("sarimax/results/metrics/metrics_all_formulas_chikungunya_", st, ".csv")
  metrics_df   <- read_csv(metrics_file, show_col_types = FALSE)

  i          <- 2
  resolved   <- FALSE

  while (i <= nrow(metrics_df) && !resolved) {
    message(sprintf("[%s] Trying formula row %d of %d: %s",
                    st, i, nrow(metrics_df), metrics_df$formula_id[i]))

    best_row <- metrics_df[i, ]

    model_data <- prepare_imdc_model_data(
      data = d,
      forecasting_climate = forecasting_climate_all,
      geocode_map = geocode_map_all,
      unit_type = "state",
      unit_id = st,
      pca = best_row$pca[1],
      pca_var_threshold = best_row$pca_var_threshold[1],
      k = best_row$k[1],
      threshold_low_variance = best_row$threshold_low_variance[1],
      threshold_cor = best_row$threshold_cor[1],
      min_cor = best_row$min_cor[1],
      max_size_covariates = best_row$max_size_covariates[1],
      index_cor_threshold = best_row$index_cor_threshold[1],
      climatology_years = best_row$climatology_years[1],
      smooth_weeks = best_row$smooth_weeks[1]
    )

    best_order   <- best_row$order[1]
    ord          <- parse_order(best_order)
    best_formula <- best_row$formula_id[1]
    if (is.na(best_formula) || identical(best_formula, "")) {
      best_formula <- character(0)
    }

    lambda_value <- if (is.na(best_row$lambda[1])) NULL else best_row$lambda[1]

    train_ids    <- paste0("train_", 1:4)
    target_ids   <- paste0("target_", 1:4)

    # Reset warnings for this attempt
    attempt_warnings <- list()

    pred_target <- lapply(seq_along(train_ids), function(j) {
      train_id  <- train_ids[j]
      target_id <- target_ids[j]
      split <- model_data$splits[[j]]

      fit <- fit_sarimax(
        data = split$data,
        formula = best_formula,
        order = c(ord$order[1], ord$order[2], ord$order[3]),
        seasonal = list(order = c(ord$seasonal_order[1], ord$seasonal_order[2], ord$seasonal_order[3]), period = 52),
        train_id = NULL,
        method = best_row$method[1],
        lambda = lambda_value,
        bootstrap = best_row$bootstrap[1],
        npaths = best_row$npaths[1],
        seed = 123 + j
      )

      w <- attr(fit, "warnings")
      if (!is.null(w)) {
        attempt_warnings[[target_id]] <<- w
        message(sprintf("[%s | %s] %d warning(s):\n%s",
                        st, target_id, length(w),
                        paste0("  - ", w, collapse = "\n")))
      }

      fit
    })
    names(pred_target) <- target_ids

    # Check if this attempt is clean (no actionable warnings)
    actionable <- unlist(lapply(attempt_warnings, function(w) {
      any(grepl("auto.arima|unreliable|failed", w, ignore.case = TRUE))
    }))

    if (length(actionable) == 0 || !any(actionable)) {
      # Clean fit — save results and move on
      resolved <- TRUE
      state_warnings_chikungunya[[st]] <- attempt_warnings

      for (j in seq_along(train_ids)) {
        write_csv(
          pred_target[[j]],
          file.path("sarimax/results/preds/",
                    paste0("pred_chikungunya_", st, "_", target_ids[j], ".csv"))
        )
      }

      preds_chikungunya_retry[[st]] <- bind_rows(pred_target, .id = "target_id") |>
        mutate(state = st)

      message(sprintf("[%s] Resolved with formula row %d: %s", st, i, best_formula))
      resolved_formulas_chikungunya <- resolved_formulas_chikungunya |>
        bind_rows(tibble(state = st, formula_id = best_formula, order = best_order, row = i))
    } else {
      message(sprintf("[%s] Formula row %d still has warnings; trying next.", st, i))
    }

    i <- i + 1
  }

  if (!resolved) {
    message(sprintf("[%s] All %d formulas exhausted; could not resolve warnings.", st, nrow(metrics_df)))
    state_warnings_chikungunya[[st]] <- attempt_warnings
  }
}

write_csv(resolved_formulas_chikungunya, "sarimax/results/metrics/resolved_formulas_chikungunya.csv")

#-----------------------------------------------------------------------------------
# City-level dengue: four validation forecasts + final 2026-2027 forecast
#-----------------------------------------------------------------------------------

best_wis_df_dengue_cities <- read_csv("sarimax/results/metrics/best_wis_dengue_all_cities.csv", show_col_types = FALSE)

preds_dengue_cities <- lapply(best_wis_df_dengue_cities$city, function(city) {
  file_name <- paste0("processed_data/dengue/sel_cities/dengue_", city, ".csv.gz")
  d <- read_csv(file_name, show_col_types = FALSE)

  best_row <- best_wis_df_dengue_cities |>
    filter(.data$city == !!city) |>
    slice(1)

  model_data <- prepare_imdc_model_data(
    data = d,
    forecasting_climate = forecasting_climate_all,
    geocode_map = geocode_map_all,
    unit_type = "city",
    unit_id = city,
    pca = best_row$pca[1],
    pca_var_threshold = best_row$pca_var_threshold[1],
    k = best_row$k[1],
    threshold_low_variance = best_row$threshold_low_variance[1],
    threshold_cor = best_row$threshold_cor[1],
    min_cor = best_row$min_cor[1],
    max_size_covariates = best_row$max_size_covariates[1],
    index_cor_threshold = best_row$index_cor_threshold[1],
    climatology_years = best_row$climatology_years[1],
    smooth_weeks = best_row$smooth_weeks[1]
  )

  best_order <- best_row |> pull(order)
  ord <- parse_order(best_order)
  best_formula <- best_row |> pull(formula_id)
  if (is.na(best_formula) || identical(best_formula, "")) {
    best_formula <- character(0)
  }

  lambda_value <- if (is.na(best_row$lambda[1])) NULL else best_row$lambda[1]

  train_ids  <- paste0("train_", 1:4)
  target_ids <- paste0("target_", 1:4)

  pred_target <- lapply(seq_along(train_ids), function(i) {
    train_id <- train_ids[i]
    target_id <- target_ids[i]
    split <- model_data$splits[[i]]

    fit <- fit_sarimax(
      data = split$data,
      formula = best_formula,
      order = c(ord$order[1], ord$order[2], ord$order[3]),
      seasonal = list(order = c(ord$seasonal_order[1], ord$seasonal_order[2], ord$seasonal_order[3]), period = 52),
      train_id = NULL,
      method = best_row$method[1],
      lambda = lambda_value,
      bootstrap = best_row$bootstrap[1],
      npaths = best_row$npaths[1],
      seed = 123 + i
    )

    write_csv(fit, file.path("sarimax/results/preds/", paste0("pred_dengue_", city, "_", target_id, ".csv")))
    fit
  })

  names(pred_target) <- target_ids

  # Final challenge forecast: cutoff EW25/2026, target EW41/2026-EW40/2027.
  final_split <- prepare_imdc_final_data(
    data = d,
    recipe = model_data,
    forecasting_climate = forecasting_climate_all,
    geocode_map = geocode_map_all,
    unit_type = "city",
    unit_id = city
  )

  final_fit <- fit_sarimax(
    data = final_split$data,
    formula = best_formula,
    order = c(ord$order[1], ord$order[2], ord$order[3]),
    seasonal = list(order = c(ord$seasonal_order[1], ord$seasonal_order[2], ord$seasonal_order[3]), period = 52),
    train_id = NULL,
    method = best_row$method[1],
    lambda = lambda_value,
    bootstrap = best_row$bootstrap[1],
    npaths = best_row$npaths[1],
    seed = 128
  )

  write_csv(final_fit, file.path("sarimax/results/preds/", paste0("pred_dengue_", city, "_final.csv")))
  pred_target[["final"]] <- final_fit

  bind_rows(pred_target, .id = "target_id") |>
    mutate(city = city)
})


#-----------------------------------------------------------------------------------
# City-level chikungunya: four validation forecasts
#-----------------------------------------------------------------------------------

best_wis_df_chikungunya_cities <- read_csv("sarimax/results/metrics/best_wis_chikungunya_all_cities.csv", show_col_types = FALSE)

preds_chikungunya_cities <- lapply(best_wis_df_chikungunya_cities$city, function(city) {
  file_name <- paste0("processed_data/chikungunya/sel_cities/chikungunya_", city, ".csv.gz")
  d <- read_csv(file_name, show_col_types = FALSE)

  best_row <- best_wis_df_chikungunya_cities |>
    filter(.data$city == !!city) |>
    slice(1)

  model_data <- prepare_imdc_model_data(
    data = d,
    forecasting_climate = forecasting_climate_all,
    geocode_map = geocode_map_all,
    unit_type = "city",
    unit_id = city,
    pca = best_row$pca[1],
    pca_var_threshold = best_row$pca_var_threshold[1],
    k = best_row$k[1],
    threshold_low_variance = best_row$threshold_low_variance[1],
    threshold_cor = best_row$threshold_cor[1],
    min_cor = best_row$min_cor[1],
    max_size_covariates = best_row$max_size_covariates[1],
    index_cor_threshold = best_row$index_cor_threshold[1],
    climatology_years = best_row$climatology_years[1],
    smooth_weeks = best_row$smooth_weeks[1]
  )

  best_order <- best_row |> pull(order)
  ord <- parse_order(best_order)
  best_formula <- best_row |> pull(formula_id)
  if (is.na(best_formula) || identical(best_formula, "")) {
    best_formula <- character(0)
  }

  lambda_value <- if (is.na(best_row$lambda[1])) NULL else best_row$lambda[1]

  train_ids  <- paste0("train_", 1:4)
  target_ids <- paste0("target_", 1:4)

  pred_target <- lapply(seq_along(train_ids), function(i) {
    train_id <- train_ids[i]
    target_id <- target_ids[i]
    split <- model_data$splits[[i]]

    fit <- fit_sarimax(
      data = split$data,
      formula = best_formula,
      order = c(ord$order[1], ord$order[2], ord$order[3]),
      seasonal = list(order = c(ord$seasonal_order[1], ord$seasonal_order[2], ord$seasonal_order[3]), period = 52),
      train_id = NULL,
      method = best_row$method[1],
      lambda = lambda_value,
      bootstrap = best_row$bootstrap[1],
      npaths = best_row$npaths[1],
      seed = 123 + i
    )

    write_csv(fit, file.path("sarimax/results/preds/", paste0("pred_chikungunya_", city, "_", target_id, ".csv")))
    fit
  })

  names(pred_target) <- target_ids
  bind_rows(pred_target, .id = "target_id") |>
    mutate(city = city)
})
