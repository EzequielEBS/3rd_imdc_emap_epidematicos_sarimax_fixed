# prep_climate.r — Step 1 of the data-prep pipeline.
#
# Gap-fills climate data for island municipalities that have no weather
# station / climate-grid coverage of their own, by copying the climate
# series of their nearest mainland municipality (with data) under the
# island's own geocode. Reads raw_data/climate.csv.gz -- combined with
# raw_data/climate_update_2026.csv.gz when that incremental file is present
# (see merge_raw_update(), below) -- and writes the completed table to
# processed_data/climate/climate.csv.gz, which is then consumed by
# data_prep/merge_dengue.r and data_prep/merge_chikungunya.r.
#
# Run before merge_dengue.r / merge_chikungunya.r.
library(tidyverse)
library(sf)
library(geobr)

#' Combine a base raw-data extract with an incremental "_update" file
#'
#' Update rows take precedence over base rows sharing the same key -- the
#' update is treated as Infodengue/Mosqlimate's revised/consolidated data
#' for the weeks it covers -- and any (key) combinations present only in
#' the update file are appended. If no update file exists at `update_path`,
#' `base` is returned completely unchanged, so this script keeps working
#' exactly as before once an update is no longer being supplied.
#'
#' @param base Data frame already read from the base raw file.
#' @param update_path Path to the corresponding "_update_2026" file.
#' @param key_cols Character vector of columns identifying a unique row
#'   (e.g. c("geocode", "epiweek")).
#' @return `base` with any updated/new rows from `update_path` merged in.
merge_raw_update <- function(base, update_path, key_cols) {
  if (!file.exists(update_path)) {
    return(base)
  }
  update <- read_csv(update_path)
  # Some raw update extracts carry a stray unnamed index column (e.g. from
  # write.csv(row.names = TRUE) upstream); drop it so bind_rows() doesn't
  # pick up an extra, all-but-useless column.
  update <- update %>% select(-any_of(c("...1", "")))
  base_kept <- base %>% anti_join(update, by = key_cols)
  bind_rows(base_kept, update)
}

climate <- read_csv("raw_data/climate.csv.gz")
climate <- merge_raw_update(climate, "raw_data/climate_update_2026.csv.gz", c("geocode", "epiweek"))

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
  select(-island_geocode)

climate_full <- bind_rows(climate, climate_islands)

if (!dir.exists("processed_data/climate")) {
  dir.create("processed_data/climate")
}

write_csv(climate_full, "processed_data/climate/climate.csv.gz")