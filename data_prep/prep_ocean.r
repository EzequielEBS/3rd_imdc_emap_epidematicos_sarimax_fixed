# prep_ocean.r — Prepare ocean-climate indices before the disease merges.
#
# Consolidates the base and 2026 update files for ENSO/IOD/PDO. When the same
# epiweek exists in both files, the update version is kept. Any remaining
# internal gaps are filled only from the most recent previous observation
# (last observation carried forward), so no future value is used.
#
# Run before merge_dengue.r / merge_chikungunya.r.
library(tidyverse)

ocean <- read_csv("raw_data/ocean_climate_oscillations.csv.gz")
ocean_update <- read_csv("raw_data/ocean_climate_oscillations_update_2026.csv.gz")

ocean <- rows_upsert(
  ocean,
  ocean_update,
  by = "epiweek"
) %>%
  arrange(epiweek) %>%
  fill(enso, iod, pdo, .direction = "down")

if (!dir.exists("processed_data/ocean")) {
  dir.create("processed_data/ocean")
}

write_csv(ocean, "processed_data/ocean/ocean_climate_oscillations.csv.gz")