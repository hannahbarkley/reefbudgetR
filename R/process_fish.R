#' Process parrotfish erosion rates from belt and spc data
#'
#'@author Rebecca Weible
#'
#'#'@param data_spc fish spc survey data.
#'@param dbase_type Erosion rates database to use. Choose either Indo-Pacific
#'ReefBudget ("dbase_type = "IPRB") or U.S. Pacific Islands rates developed
#'by Tye Kindinger, NOAA PIFSC ("dbase_type = "Kindinger").
#'@param data_belt fish belt survey data.
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

process_fish <- function(spc_data= data_spc,
                         dbase_type = c("Kindinger", "IPRB"),
                         belt_data = data_belt,
                         subset_distance_m) {
  
  if(dbase_type == "Kindinger") {
    rates_dbase <- fish_erosion_dbase_kindinger
  } else{
    rates_dbase <- fish_erosion_dbase_iprb
  }
  
  # if you have SPC data then....
  if(!missing(spc_data)){
  
    fish_fixed_spc <- calc_fish_fixed_spc(
                        data = data_spc, 
                        rates_dbase = rates_dbase)
    
    
    fish_strs_spc_ <- calc_fish_strs_spc(
                      data = data_spc, 
                      rates_dbase = rates_dbase,
                      subset_distance_m)
    fish_strs_spc <- fish_strs_spc_$calc_strs_ero
    
    
    if (!missing(belt_data)) {
      
      fish_belt_ <- calc_fish_belt(
                    data = data_belt, 
                    rates_dbase = rates_dbase, 
                    full_summary = TRUE)
      
        return(list(dat = rbind(fish_belt_$fish_erosion_site, 
                                  fish_fixed_spc,
                                  fish_strs_spc),
                    assoc_site_count = fish_strs_spc_$assoc_survey_count))
    }
    
    
    
    if (missing(belt_data)) {
        return(list(dat = rbind(fish_fixed_spc,
                                  fish_strs_spc),
                    assoc_site_count = fish_strs_spc_$assoc_survey_count))
    }
  }
  
  # If missing SPC data then...
  else {
    
    fish_belt_ <- calc_fish_belt(
      data = data_belt, 
      rates_dbase = rates_dbase, 
      full_summary = TRUE)
    
    return(fish_belt_)
    
  }

  
}


