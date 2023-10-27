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
#'@export process_fish
#'
#'@examples
#'fish_data <- read.csv("CB_FishBelt_alldata.csv", na = "", check.names = FALSE)
#'
#'fish_belt <- calc_fish_belt(data = fish_data, dbase_type = "Kindinger",
#'full_summary = TRUE)

calc_fish_belt <- function(data,
                           dbasetype = rates_dbase,
                           full_summary = TRUE) {
  
    # FOR BELT DATA ----------------------------------------------------------------
    
    # Calculate erosion rates per fish -------------------------------------------
    
    calc_eros_fish_output <- calc_eros_fish(data, rates_dbase)
    
    # Calculate bioerosion metrics per grazing type per site ---------------------
    
    density_average <- suppressWarnings(
      summarize_fish_metrics(
        data = calc_eros_fish_output,
        metric = "density",
        level = "transect",
        summarize_by = "species"
      )
    ) %>%
      add_column(METRIC = "FISH_DENSITY_ABUNDANCE_HA")
    
    biomass_average <- suppressWarnings(
      summarize_fish_metrics(
        data = calc_eros_fish_output,
        metric = "biomass",
        level = "transect",
        summarize_by = "species"
      )
    ) %>%
      tibble::add_column(METRIC = "FISH_BIOMASS_KG_HA")
    
    bioerosion_average <- suppressWarnings(
      summarize_fish_metrics(
        data = calc_eros_fish_output,
        metric = "bioerosion",
        level = "transect",
        summarize_by = "species"
      )
    ) %>%
      tibble::add_column(METRIC = "FISH_EROSION_KG_M2_YR")
    
    
    species_table <-
      rbind(density_average, biomass_average, bioerosion_average)
    
    
    summary_belt_erosion <- summarize_fish_erosion(species_table, full_summary)
    
    
    return(summary_belt_erosion)
    
    if (full_summary == TRUE) {
      return(list(
        fish_erosion_transect = summary_belt_erosion$fish_erosion_transect,
        fish_erosion_site = summary_belt_erosion$fish_erosion_site)
      )
    }
    
    if (full_summary == FALSE) {
      return(summary_belt_erosion$fish_erosion_site)
    }
  
  
}
