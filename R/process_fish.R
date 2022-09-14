#' Calculate parrotfish erosion rates from belt data
#'
#'@author Rebecca Weible
#'
#'@param data Parrotfish belt data.
#'@param rates_dbase Erosion rates database to use. Choose either Indo-Pacific
#'ReefBudget ("rates_dbase = "IPRB") or U.S. Pacific Islands rates developed
#'by Tye Kindinger, NOAA PIFSC ("rates_dbase = "Kindinger").
#'
#'@import Rmisc
#'@import tidyverse
#'@import dplyr
#'
#'@export process_fish

process_fish <- function(data,
                         rates_dbase = "Kindinger",
                         full_summary = FALSE) {

  # Calculate erosion rates per fish -------------------------------------------

  calc_eros_fish_output <- calc_eros_fish(data,
                                          rates_dbase = rates_dbase)

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
    add_column(METRIC = "FISH_BIOMASS_KG_HA")

  bioerosion_average <- suppressWarnings(
    summarize_fish_metrics(
      data = calc_eros_fish_output,
      metric = "bioerosion",
      level = "transect",
      summarize_by = "species"
    )
  ) %>%
    add_column(METRIC = "FISH_EROSION_KG_M2_YR")


  species_table <-
    rbind(density_average, biomass_average, bioerosion_average)


  summary_fish_erosion <- summarize_fish_erosion(species_table, full_summary)

  return(summary_fish_erosion)

  if (full_summary == TRUE) {
    return(list(
      fish_erosion_transect = fish_erosion_transect,
      fish_erosion_site = fish_erosion_site)
    )
  }

  if (full_summary == FALSE) {
    return(summary_fish_erosion)
  }


}
