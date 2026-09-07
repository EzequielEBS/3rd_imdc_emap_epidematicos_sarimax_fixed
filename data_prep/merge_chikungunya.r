# merge_chikungunya.r — Step 2 of the data-prep pipeline (chikungunya branch).
#
# Mirrors merge_dengue.r exactly, but for the chikungunya case series: joins
# climate, environmental, ocean-climate index, and population covariates
# onto the raw chikungunya table and writes
# processed_data/chikungunya/chikungunya_merged.csv.gz. chikungunya.csv.gz
# and ocean_climate_oscillations.csv.gz are each combined with their
# "_update_2026" counterpart when present (see merge_raw_update(), below);
# climate.csv.gz was already combined with its own update by
# prep_climate.r, so it needs no merging here.
#
# Run after prep_climate.r; followed by data_prep/handle_na.r.
library(tidyverse)
library(lubridate)
library(aweek)

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

chikungunya <- read_csv("raw_data/chikungunya.csv.gz")
chikungunya <- merge_raw_update(chikungunya, "raw_data/chikungunya_update_2026.csv.gz", c("geocode", "epiweek"))
climate <- read_csv("processed_data/climate/climate.csv.gz")
env_vars <- read_csv("raw_data/environ_vars.csv.gz")
ocean <- read_csv("raw_data/ocean_climate_oscillations.csv.gz")
# Merge the ocean-oscillation update on `date` (not `epiweek`) because the
# epiweek this script actually uses is recomputed from `date` just below --
# joining on the update file's own epiweek column first could disagree with
# that recomputation.
ocean <- merge_raw_update(ocean, "raw_data/ocean_climate_oscillations_update_2026.csv.gz", c("date"))
pop <- read_csv("raw_data/datasus_population_2001_2025.csv.gz")
map_regional_health <- read_csv("raw_data/map_regional_health.csv")

# Prepare data to merge: derive join keys and drop columns that would
# otherwise collide across tables (each table keeps only one `date`).
chikungunya <- chikungunya %>%
  mutate(year = year(date))
climate <- climate |> dplyr::select(-date)
ocean <- ocean %>%
  mutate(
    epiweek = as.integer(
      paste0(epiyear(date), sprintf("%02d", epiweek(date)))
    )
  ) |>
  dplyr::select(-date)
# Extend population one year past the last available DATASUS estimate (2025)
# by carrying the 2025 value forward to 2026.
pop <- rbind(
  pop,
  lapply(unique(pop$geocode), function(code) {
    data.frame(
      geocode = code,
      year = 2026,
      population = pop$population[pop$geocode == code & pop$year == 2025]
    )
  }) %>% bind_rows()
)

# Merge data: left-join everything onto the chikungunya case series so every
# case row is preserved even if a covariate table is missing that key.
chikungunya_merged <- chikungunya %>%
  left_join(climate, by = c("epiweek", "geocode")) %>%
  left_join(env_vars, by = c("geocode", "uf_code")) %>%
  left_join(ocean, by = "epiweek") %>%
  left_join(pop, by = c("geocode", "year"))

# Count NAs by column
na_counts <- sapply(chikungunya_merged, function(x) sum(is.na(x)))

if (!dir.exists("processed_data/chikungunya")) {
  dir.create("processed_data/chikungunya")
}
write_csv(chikungunya_merged, "processed_data/chikungunya/chikungunya_merged.csv.gz")
