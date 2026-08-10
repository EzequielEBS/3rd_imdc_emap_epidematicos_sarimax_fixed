# sel_cities.r — Optional data-prep step: extract a focal set of cities.
#
# Carves out the selected geocode-level series (one CSV per city) from the
# already-merged dengue/chikungunya tables. The same climate aliases and
# historical lag/rolling columns used by the state-level SARIMAX inputs are
# added here so the city files expose the same candidate-covariate names.
#
# The cutoff-specific model step must still rebuild future covariates for each
# backtest; these are the historical master features only.

library(tidyverse)
library(runner)

# IBGE geocodes of state-capital/large cities tracked for dengue diagnostics.
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

# IBGE geocodes of state-capital/large cities tracked for chikungunya diagnostics.
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

dengue_merged <- read_csv("processed_data/dengue/dengue_merged.csv.gz")
chikungunya_merged <- read_csv("processed_data/chikungunya/chikungunya_merged.csv.gz")

# Verify that the hard-coded challenge lists still match the official
# target_city flags carried by the epidemiological datasets.
stopifnot(setequal(
  cities_dengue,
  dengue_merged %>% filter(target_city) %>% distinct(geocode) %>% pull(geocode)
))

stopifnot(setequal(
  cities_chikungunya,
  chikungunya_merged %>% filter(target_city) %>% distinct(geocode) %>% pull(geocode)
))

add_city_features <- function(data) {
  data <- data %>%
    arrange(geocode, epiweek) %>%
    group_by(geocode) %>%
    mutate(
      cases = casos,
      temp_min_mean = temp_min,
      temp_med_mean = temp_med,
      temp_max_mean = temp_max,
      precip_min_mean = precip_min,
      precip_med_mean = precip_med,
      precip_max_mean = precip_max,
      pressure_min_mean = pressure_min,
      pressure_med_mean = pressure_med,
      pressure_max_mean = pressure_max,
      rel_humid_min_mean = rel_humid_min,
      rel_humid_med_mean = rel_humid_med,
      rel_humid_max_mean = rel_humid_max,
      thermal_range_mean = thermal_range,
      rainy_days_mean = rainy_days,
      pop = population
    ) %>%
    ungroup()

  data <- data %>%
    group_by(geocode) %>%
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

  data <- data %>%
    group_by(geocode) %>%
    mutate(
      across(
        c(temp_min_mean, temp_med_mean, temp_max_mean, precip_min_mean, precip_med_mean, precip_max_mean,
          pressure_min_mean, pressure_med_mean, pressure_max_mean, rel_humid_min_mean,
          rel_humid_med_mean, rel_humid_max_mean, thermal_range_mean, rainy_days_mean),
        list(
          mean_3mo = ~ dplyr::lag(mean_run(., k = 12, na_rm = TRUE)),
          mean_6mo = ~ dplyr::lag(mean_run(., k = 24, na_rm = TRUE)),
          mean_9mo = ~ dplyr::lag(mean_run(., k = 36, na_rm = TRUE)),
          mean_12mo = ~ dplyr::lag(mean_run(., k = 48, na_rm = TRUE))
        ),
        .names = "{col}_{fn}"
      )
    ) %>%
    ungroup()

  data
}

dengue_merged <- add_city_features(dengue_merged)
chikungunya_merged <- add_city_features(chikungunya_merged)

lapply(cities_dengue, function(city) {
  if (!dir.exists("processed_data/dengue/sel_cities")) {
    dir.create("processed_data/dengue/sel_cities")
  }
  dengue_merged %>%
    filter(geocode == city) %>%
    write_csv(paste0("processed_data/dengue/sel_cities/dengue_", city, ".csv.gz"))
})

lapply(cities_chikungunya, function(city) {
  if (!dir.exists("processed_data/chikungunya/sel_cities")) {
    dir.create("processed_data/chikungunya/sel_cities")
  }
  chikungunya_merged %>%
    filter(geocode == city) %>%
    write_csv(paste0("processed_data/chikungunya/sel_cities/chikungunya_", city, ".csv.gz"))
})
