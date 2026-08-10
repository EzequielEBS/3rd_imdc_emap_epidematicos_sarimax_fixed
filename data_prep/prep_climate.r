# prep_climate.r — Step 1 of the data-prep pipeline.
#
# Gap-fills climate data for island municipalities that have no weather
# station / climate-grid coverage of their own, by copying the climate
# series of their nearest mainland municipality (with data) under the
# island's own geocode. Reads raw_data/climate.csv.gz and writes the
# completed table to processed_data/climate/climate.csv.gz, which is then
# consumed by data_prep/merge_dengue.r and data_prep/merge_chikungunya.r.
#
# Run before merge_dengue.r / merge_chikungunya.r.
library(tidyverse)
library(sf)
library(geobr)

climate <- read_csv("raw_data/climate.csv.gz")

# Update rule: keep base-only rows, add update-only rows, and prefer the
# update version when the same geocode/epiweek exists in both files.
climate <- rows_upsert(
  climate,
  climate_update,
  by = c("geocode", "epiweek")
)

forecasting_climate <- read_csv("raw_data/forecasting_climate.csv.gz")
forecasting_climate_update <- read_csv("raw_data/forecasting_climate_update_2026.csv.gz")

forecasting_climate <- rows_upsert(
  forecasting_climate,
  forecasting_climate_update,
  by = c("geocode", "reference_month", "forecast_months_ahead")
)

# Download municipality boundaries (2020 IBGE mesh) to get coordinates for
# every geocode, since the raw climate data has no spatial reference of its
# own beyond the geocode itself.
munis <- read_municipality(year = 2020) %>%
  mutate(geocode = as.integer(code_muni))

# Get centroids
munis <- munis %>%
  mutate(centroid = st_centroid(geom),
         lon = st_coordinates(centroid)[,1],
         lat = st_coordinates(centroid)[,2])

# Known island municipalities without their own climate series: Fernando de
# Noronha (PE), Ilhabela (SP), and Florianópolis-area Ilha de
# Santa Catarina/Ilha do Mel-type codes. Hard-coded because there is no
# generic "is an island" flag in the IBGE mesh.
island_codes <- c(2916104, 2605459, 2919926)

islands   <- munis %>% filter(geocode %in% island_codes)
mainland  <- munis %>% filter(!geocode %in% island_codes)

# Further restrict mainland candidates to geocodes that actually have
# climate data, so nearest-neighbor matching never picks a donor with no
# series of its own.
mainland  <- mainland %>% filter(geocode %in% unique(climate$geocode))

# For each island, find its closest mainland municipality (by centroid).
nearest_idx <- st_nearest_feature(islands, mainland)

island_to_mainland <- tibble(
  island_geocode   = islands$geocode,
  mainland_geocode = mainland$geocode[nearest_idx]
)

# Copy the matched mainland municipality's climate series onto the island's
# own geocode, then append it to the original table.
climate_islands <- climate %>%
  filter(geocode %in% island_to_mainland$mainland_geocode) %>%
  left_join(island_to_mainland, by = c("geocode" = "mainland_geocode")) %>%
  mutate(geocode = island_geocode) %>%
  select(-island_geocode) %>%
  anti_join(climate %>% select(geocode, epiweek), by = c("geocode", "epiweek"))

climate_full <- bind_rows(climate, climate_islands)

if (!dir.exists("processed_data/climate")) {
  dir.create("processed_data/climate")
}

write_csv(climate_full, "processed_data/climate/climate.csv.gz")
write_csv(forecasting_climate, "processed_data/climate/forecasting_climate.csv.gz")