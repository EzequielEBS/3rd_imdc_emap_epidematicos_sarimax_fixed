# agg_data_uf.r — Step 4 of the data-prep pipeline: city → state aggregation
# and feature engineering.
#
# This is the script that actually produces the modeling-ready tables
# consumed by sarimax/src/*.r. For each disease it: (1) sums case counts and
# averages weather covariates from city/geocode level up to state (uf)
# level; (2) creates lagged (4/8/12/16-week) and rolling-mean
# (3/6/9/12-month) versions of every weather covariate, to let the SARIMAX
# models use past climate as a leading indicator of future transmission;
# and (3) writes one CSV per state to processed_data/<disease>/<disease>_<uf>_agg.csv.gz
# (the files read by model_sel.r and fit.r). Espírito Santo (ES) is
# excluded from the state-level outputs.
#
# Run after data_prep/handle_na.r. This is the last data-prep step before
# modeling.
library(tidyverse)
library(runner)

dengue_merged <- read_csv("processed_data/dengue/dengue_merged.csv.gz")
chikungunya_merged <- read_csv("processed_data/chikungunya/chikungunya_merged.csv.gz")

#' Most frequent value in a vector (statistical mode)
#'
#' Used to aggregate the categorical Köppen climate classification and biome
#' columns from city level up to a single representative state-level value
#' (the most common category among that state's cities).
#'
#' @param x A vector (typically character/factor).
#' @return The single most frequent value in `x`.
mode_stat <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# ── Aggregate dengue from city to state level ──────────────────────────────
# Cases are summed across cities; weather covariates are averaged; the
# categorical Köppen/biome columns use mode_stat(); ocean-climate indices,
# population, year and the train/target indicators are constant within a
# state/epiweek so first() simply picks them up unchanged.
dengue_uf <- dengue_merged %>%
  group_by(uf, uf_code, date, epiweek) %>%
  summarise(
    cases = sum(casos, na.rm = TRUE),
    temp_min_mean = mean(temp_min, na.rm = TRUE),
    temp_med_mean = mean(temp_med, na.rm = TRUE),
    temp_max_mean = mean(temp_max, na.rm = TRUE),
    precip_min_mean = mean(precip_min, na.rm = TRUE),
    precip_med_mean = mean(precip_med, na.rm = TRUE),
    precip_max_mean = mean(precip_max, na.rm = TRUE),
    pressure_min_mean = mean(pressure_min, na.rm = TRUE),
    pressure_med_mean = mean(pressure_med, na.rm = TRUE),
    pressure_max_mean = mean(pressure_max, na.rm = TRUE),
    rel_humid_min_mean = mean(rel_humid_min, na.rm = TRUE),
    rel_humid_med_mean = mean(rel_humid_med, na.rm = TRUE),
    rel_humid_max_mean = mean(rel_humid_max, na.rm = TRUE),
    thermal_range_mean = mean(thermal_range, na.rm = TRUE),
    rainy_days_mean = mean(rainy_days, na.rm = TRUE),
    koppen = mode_stat(koppen),
    biome = mode_stat(biome),
    enso = first(enso),
    iod = first(iod),
    pdo = first(pdo),
    pop = sum(population, na.rm = TRUE),
    year = first(year),
    train_1 = first(train_1),
    train_2 = first(train_2),
    train_3 = first(train_3),
    train_4 = first(train_4),
    target_1 = first(target_1),
    target_2 = first(target_2),
    target_3 = first(target_3),
    target_4 = first(target_4)
  ) %>%
  ungroup()

# ── Aggregate chikungunya from city to state level (mirrors dengue above) ──
chikungunya_uf <- chikungunya_merged %>%
  group_by(uf, uf_code, date, epiweek) %>%
  summarise(
    cases = sum(casos, na.rm = TRUE),
    temp_min_mean = mean(temp_min, na.rm = TRUE),
    temp_med_mean = mean(temp_med, na.rm = TRUE),
    temp_max_mean = mean(temp_max, na.rm = TRUE),
    precip_min_mean = mean(precip_min, na.rm = TRUE),
    precip_med_mean = mean(precip_med, na.rm = TRUE),
    precip_max_mean = mean(precip_max, na.rm = TRUE),
    pressure_min_mean = mean(pressure_min, na.rm = TRUE),
    pressure_med_mean = mean(pressure_med, na.rm = TRUE),
    pressure_max_mean = mean(pressure_max, na.rm = TRUE),
    rel_humid_min_mean = mean(rel_humid_min, na.rm = TRUE),
    rel_humid_med_mean = mean(rel_humid_med, na.rm = TRUE),
    rel_humid_max_mean = mean(rel_humid_max, na.rm = TRUE),
    thermal_range_mean = mean(thermal_range, na.rm = TRUE),
    rainy_days_mean = mean(rainy_days, na.rm = TRUE),
    koppen = mode_stat(koppen),
    biome = mode_stat(biome),
    enso = first(enso),
    iod = first(iod),
    pdo = first(pdo),
    pop = sum(population, na.rm = TRUE),
    year = first(year),
    train_1 = first(train_1),
    train_2 = first(train_2),
    train_3 = first(train_3),
    train_4 = first(train_4),
    target_1 = first(target_1),
    target_2 = first(target_2),
    target_3 = first(target_3),
    target_4 = first(target_4)
  ) %>%
  ungroup()

# create lagged covariates
dengue_uf <- dengue_uf %>%
  arrange(epiweek) %>%
  group_by(uf) %>%
  mutate(
    across(
      c(temp_min_mean, temp_med_mean, temp_max_mean, precip_min_mean, precip_med_mean, precip_max_mean,
        pressure_min_mean, pressure_med_mean, pressure_max_mean, rel_humid_min_mean,
        rel_humid_med_mean, rel_humid_max_mean, thermal_range_mean, rainy_days_mean),
      list(
        lag4 = ~ lag(., 4),
        lag8 = ~ lag(., 8),
        lag12 = ~ lag(., 12),
        lag16 = ~ lag(., 16)
      ),
      .names = "{col}_{fn}"
    )
  ) %>%
  ungroup()

chikungunya_uf <- chikungunya_uf %>%
  arrange(epiweek) %>%
  group_by(uf) %>%
  mutate(
    across(
      c(temp_min_mean, temp_med_mean, temp_max_mean, precip_min_mean, precip_med_mean, precip_max_mean,
        pressure_min_mean, pressure_med_mean, pressure_max_mean, rel_humid_min_mean,
        rel_humid_med_mean, rel_humid_max_mean, thermal_range_mean, rainy_days_mean),
      list(
        lag4 = ~ lag(., 4),
        lag8 = ~ lag(., 8),
        lag12 = ~ lag(., 12),
        lag16 = ~ lag(., 16)
      ),
      .names = "{col}_{fn}"
    )
  ) %>%
  ungroup()

# create summary covariates
dengue_uf <- dengue_uf %>%
  arrange(epiweek) %>%
  group_by(uf) %>%
  mutate(
    across(
      c(temp_min_mean, temp_med_mean, temp_max_mean, precip_min_mean, precip_med_mean, precip_max_mean,
        pressure_min_mean, pressure_med_mean, pressure_max_mean, rel_humid_min_mean,
        rel_humid_med_mean, rel_humid_max_mean, thermal_range_mean, rainy_days_mean),
      list(
        # Calculates the mean of the previous 12 weeks
        mean_3mo = ~ dplyr::lag(mean_run(., k = 12, na_rm = TRUE)),
        mean_6mo = ~ dplyr::lag(mean_run(., k = 24, na_rm = TRUE)),
        mean_9mo = ~ dplyr::lag(mean_run(., k = 36, na_rm = TRUE)),
        mean_12mo = ~ dplyr::lag(mean_run(., k = 48, na_rm = TRUE))
      ),
      .names = "{col}_{fn}"
    )
  ) %>%
  ungroup()

# ── Summary (rolling-mean) covariates for chikungunya (mirrors dengue above) ──
chikungunya_uf <- chikungunya_uf %>%
  arrange(epiweek) %>%
  group_by(uf) %>%
  mutate(
    across(
      c(temp_min_mean, temp_med_mean, temp_max_mean, precip_min_mean, precip_med_mean, precip_max_mean,
        pressure_min_mean, pressure_med_mean, pressure_max_mean, rel_humid_min_mean,
        rel_humid_med_mean, rel_humid_max_mean, thermal_range_mean, rainy_days_mean),
      list(
        # Calculates the mean of the previous 12 weeks
        mean_3mo = ~ dplyr::lag(mean_run(., k = 12, na_rm = TRUE)),
        mean_6mo = ~ dplyr::lag(mean_run(., k = 24, na_rm = TRUE)),
        mean_9mo = ~ dplyr::lag(mean_run(., k = 36, na_rm = TRUE)),
        mean_12mo = ~ dplyr::lag(mean_run(., k = 48, na_rm = TRUE))
      ),
      .names = "{col}_{fn}"
    )
  ) %>%
  ungroup()


# count NAs by column
na_counts_dengue <- sapply(dengue_uf, function(x) sum(is.na(x)))
na_counts_chikungunya <- sapply(chikungunya_uf, function(x) sum(is.na(x)))

lapply(unique(dengue_uf$uf), function(uf) {
  if (uf == "ES") {
    return()  # Skip ES
  }
  dengue_uf %>%
    filter(uf == !!uf) %>%
    write_csv(paste0("processed_data/dengue/dengue_", uf, "_agg.csv.gz"))
})

lapply(unique(chikungunya_uf$uf), function(uf) {
  if (uf == "ES") {
    return()  # Skip ES
  }
  chikungunya_uf %>%
    filter(uf == !!uf) %>%
    write_csv(paste0("processed_data/chikungunya/chikungunya_", uf, "_agg.csv.gz"))
})
