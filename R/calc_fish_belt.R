#' Calculate parrotfish erosion rates from belt data
#'
#'@author Rebecca Weible
#'
#'@param data fish survey data.
#'@param rates_dbase Erosion rates database to use.
#'
#'@import Rmisc
#'@import tidyverse
#'@import tibble
#'@import dplyr
#'
#'@export calc_fish_belt
#'
#'@examples
#'fish_data <- read.csv("CB_FishBelt_alldata.csv", na = "", check.names = FALSE)
#'
#'fish_belt <- calc_fish_belt(data = fish_data, rates_dbase = rates_dbase,
#'full_summary = TRUE)

calc_fish_belt <- function(data, 
                           rates_dbase, 
                           fixed_metadata = sites_metadata,
                           full_summary = TRUE) {
  
  # Calculate erosion rates per fish ----------------------------------------
  
  data <- prep_spc(data, fixed_metadata)
  
  calc_eros_fish_output <- calc_eros_fish(data, rates_dbase)
  
  # Calculate metrics ---------------------

  metrics_map <- list(
    density    = "FISH_DENSITY_ABUNDANCE_HA",
    biomass    = "FISH_BIOMASS_KG_HA",
    bioerosion = "FISH_EROSION_KG_M2_YR"
  )
  
  # Iterate over the map and bind results into one dataframe
  species_table <- purrr::imap_dfr(metrics_map, function(label, metric_key) {
    suppressWarnings(
      summarize_fish_metrics(
        data = calc_eros_fish_output,
        metric = metric_key,
        level = "transect",
        summarize_by = "species"
      )
    ) %>%
      dplyr::mutate(METRIC = label)
  })
  
  # Final Summary -----------------------------------------------------------
  summary_belt_erosion <- summarize_fish_erosion(species_table, full_summary)
  
  # Return Logic ( --------------------------
  if (full_summary) {
    return(list(
      fish_erosion_transect = summary_belt_erosion$fish_erosion_transect,
      fish_erosion_site     = summary_belt_erosion$fish_erosion_site
    ))
  } else {
    return(summary_belt_erosion$fish_erosion_site)
  }
}
