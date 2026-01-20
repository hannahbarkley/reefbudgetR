#' Process parrotfish erosion rates from belt and spc data
#'
#'@author Rebecca Weible
#'
#'@param fish_data Fish spc survey data.
#'@param method_type Type of SPC data. "Fixed" (fixed only), "StRS" (stratified random), "Mean StRS" (average of nearby stratified random), or "Belt" (belt data only)
#'@param fixed_metadata Dataframe with fixed site info
#'@param dbase_type Erosion rates database to use. Choose either Indo-Pacific
#'ReefBudget ("dbase_type = "IPRB") or U.S. Pacific Islands rates developed
#'by Tye Kindinger, NOAA PIFSC ("dbase_type = "Kindinger").
#'@param subset_distance_m Assigned associated site distances, in meters, from fixed SPC/OCC site to all other fish SPC sites, based on parrotfish foraging boundaries.
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
#'fish_belt <- process_fish(spc_data = fish_data_spc,
#'dbase_type = "Kindinger", belt_data = fish_data_belt, subset_distance_m = 2000)

process_fish <- function(data,
                         method_type = "Fixed",                         
                         fixed_metadata = NULL,
                         dbase_type = "Kindinger",
                         subset_distance_m = 6000) {
  
  if (dbase_type == "Kindinger") {
    rates_dbase <- fish_erosion_dbase_kindinger
  } else{
    rates_dbase <- fish_erosion_dbase_iprb
  }
  
  
  # if you have SPC data then....
  if (method_type == "Fixed") {
      fish_fixed_spc <- calc_fish_fixed_spc(data = data, rates_dbase = rates_dbase)
      
      return(fish_fixed_spc)
  }
    
  
  if (method_type == "StRS") {
      fish_strs_spc_ <- calc_fish_strs_spc(data = data,
                                           rates_dbase = rates_dbase,
                                           method_type = "StRS")
      
      return(fish_strs_spc_)
  }
  
  if (method_type == "Mean StRS") {
    fish_strs_spc_ <- calc_fish_strs_spc(data = data,
                                         fixed_metadata = fixed_metadata,
                                         rates_dbase = rates_dbase,
                                         subset_distance_m,
                                         method_type = "Mean StRS")
    
    return(fish_strs_spc_)
  }
    
    
  if (method_type == "Belt") {
      fish_belt_ <- calc_fish_belt(data = data, 
                                    rates_dbase = rates_dbase, 
                                    full_summary = TRUE)
      return(fish_belt_)
    }
    
}
