# merge_dengue.r — Step 2 of the data-prep pipeline (dengue branch).
#
# Joins the raw dengue case series (geocode x epiweek, with train/target
# split indicators already provided by Mosqlimate) against climate,
# environmental (Köppen climate classification, biome), ocean-climate index
# (ENSO/IOD/PDO), and population covariates. Population is extrapolated one
# year forward (to 2026) by repeating each geocode's 2025 value, since
# official population estimates lag the current year. Writes
# processed_data/dengue/dengue_merged.csv.gz.
#
# Run after prep_climate.r; followed by data_prep/handle_na.r.
library(tidyverse)
library(lubridate)
library(aweek)

dengue <- read_csv("raw_data/dengue.csv.gz")

#If the official update file is added to raw_data, use it with the same update rule
# applied to the other 2026 update files.
if (file.exists("raw_data/dengue_update_2026.csv.gz")) {
  dengue_update <- read_csv("raw_data/dengue_update_2026.csv.gz")
  dengue <- rows_upsert(
    dengue,
    dengue_update,
    by = c("geocode", "epiweek")
  )
}

climate <- read_csv("processed_data/climate/climate.csv.gz")
env_vars <- read_csv("raw_data/environ_vars.csv.gz")
ocean <- read_csv("raw_data/ocean_climate_oscillations.csv.gz")
pop <- read_csv("raw_data/datasus_population_2001_2025.csv.gz")
map_regional_health <- read_csv("raw_data/map_regional_health.csv")

# Prepare data to merge: derive join keys and drop columns that would
# otherwise collide across tables (each table keeps only one `date`).
dengue <- dengue %>%
  mutate(year = epiweek %/% 100)
climate <- climate |> dplyr::select(-date)
ocean <- ocean |> dplyr::select(-date)

# Extend population two years past the last available DATASUS estimate (2025)
# by carrying the 2025 value forward to 2026 and 2027.
pop_2026 <- pop %>%
  filter(year == 2025) %>%
  mutate(year = 2026)

pop_2027 <- pop %>%
  filter(year == 2025) %>%
  mutate(year = 2027)

pop <- bind_rows(pop, pop_2026, pop_2027)

# Merge data: left-join everything onto the dengue case series so every
# case row is preserved even if a covariate table is missing that key.
dengue_merged <- dengue %>%
  left_join(climate, by = c("epiweek", "geocode")) %>%
  left_join(env_vars, by = c("geocode", "uf_code")) %>%
  left_join(ocean, by = "epiweek") %>%
  left_join(pop, by = c("geocode", "year"))

# Count NAs by column
na_counts <- sapply(dengue_merged, function(x) sum(is.na(x)))

if (!dir.exists("processed_data/dengue")) {
  dir.create("processed_data/dengue")
}

write_csv(dengue_merged, "processed_data/dengue/dengue_merged.csv.gz")
