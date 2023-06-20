#' Process parrotfish erosion rates from belt and spc data
#'
#'@author Rebecca Weible
#'
#'@param data_belt fish belt survey data.
#'#'@param data_spc fish spc survey data.
#'@param dbase_type Erosion rates database to use. Choose either Indo-Pacific
#'ReefBudget ("dbase_type = "IPRB") or U.S. Pacific Islands rates developed
#'by Tye Kindinger, NOAA PIFSC ("dbase_type = "Kindinger").
#'#'@param sites_associated Location data was collected. Choose either Oahu ("sites_associated = "OAH"),
#'or Mariana Islands ("sites_associated = "MARIAN").
#'
#'@import Rmisc
#'@import tidyverse
#'@import dplyr
#'
#'@export process_fish
#'
#'@examples
#'fish_data <- read.csv("CB_FishBelt_alldata.csv", na = "", check.names = FALSE)
#'
#'fish_belt <- process_fish(data_belt = fish_data_belt, data_spc = fish_data_spc, 
#'dbase_type = "Kindinger", sites_associated = "OAH")

process_fish <- function(data_belt,
                         data_spc,
                         dbase_type = c("IPRB", "Kindinger"),
                         sites_associated = c("OAH", "MARIAN")) {
  

  
  fish_belt <- calc_fish_belt(
                data = data_belt, 
                dbase_type = dbase_type, 
                full_summary = TRUE)
  
  
  fish_fixed_spc <- calc_fish_fixed_spc(
                      data = data_spc, 
                      dbase_type = dbase_type, 
                      sites_associated = sites_associated)
  
  
  fish_strs_spc <- calc_fish_strs_spc(
                    data = data_spc, 
                    dbase_type = dbase_type, 
                    sites_associated = sites_associated)
  
  
  
  if (sites_associated == "OAH") {
    return(rbind(fish_belt$fish_erosion_site, 
                 fish_fixed_spc %>%
                   # If Excavator functional group is missing:
                   add_column(FISH_BIOMASS_KG_HA_EXCAVATOR_L95 = 0,
                              FISH_BIOMASS_KG_HA_EXCAVATOR_MEAN = 0,
                              FISH_BIOMASS_KG_HA_EXCAVATOR_SD = 0,
                              FISH_BIOMASS_KG_HA_EXCAVATOR_SE = 0,
                              FISH_BIOMASS_KG_HA_EXCAVATOR_U95 = 0,
                              FISH_DENSITY_ABUNDANCE_HA_EXCAVATOR_L95 = 0,
                              FISH_DENSITY_ABUNDANCE_HA_EXCAVATOR_MEAN = 0,
                              FISH_DENSITY_ABUNDANCE_HA_EXCAVATOR_SD = 0,
                              FISH_DENSITY_ABUNDANCE_HA_EXCAVATOR_SE = 0,
                              FISH_DENSITY_ABUNDANCE_HA_EXCAVATOR_U95 = 0,
                              FISH_EROSION_KG_M2_YR_EXCAVATOR_L95 = 0,
                              FISH_EROSION_KG_M2_YR_EXCAVATOR_MEAN = 0,
                              FISH_EROSION_KG_M2_YR_EXCAVATOR_SD = 0,
                              FISH_EROSION_KG_M2_YR_EXCAVATOR_SE = 0,
                              FISH_EROSION_KG_M2_YR_EXCAVATOR_U95 = 0) %>%
                   select(REGION:FISH_BIOMASS_KG_HA_ALL_U95, FISH_BIOMASS_KG_HA_EXCAVATOR_L95:FISH_BIOMASS_KG_HA_EXCAVATOR_U95, FISH_BIOMASS_KG_HA_OTHER_L95:FISH_BIOMASS_KG_HA_SCRAPER_U95, 
                          FISH_DENSITY_ABUNDANCE_HA_ALL_L95:FISH_DENSITY_ABUNDANCE_HA_ALL_U95, FISH_DENSITY_ABUNDANCE_HA_EXCAVATOR_L95:FISH_DENSITY_ABUNDANCE_HA_EXCAVATOR_U95, 
                          FISH_DENSITY_ABUNDANCE_HA_OTHER_L95:FISH_DENSITY_ABUNDANCE_HA_SCRAPER_U95, FISH_EROSION_KG_M2_YR_ALL_L95:FISH_EROSION_KG_M2_YR_ALL_U95, 
                          FISH_EROSION_KG_M2_YR_EXCAVATOR_L95:FISH_EROSION_KG_M2_YR_EXCAVATOR_U95, FISH_EROSION_KG_M2_YR_OTHER_L95:FISH_EROSION_KG_M2_YR_SCRAPER_U95),
                 fish_strs_spc))
  }
  
  if (sites_associated == "MARIAN") {
    return(rbind(fish_belt$fish_erosion_site, fish_fixed_spc, fish_strs_spc)
    )
  }

  
}


