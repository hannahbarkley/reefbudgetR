#' Process parrotfish erosion rates from belt and spc data
#'
#'@author Rebecca Weible
#'
#'#'@param data_spc fish spc survey data.
#'@param dbase_type Erosion rates database to use. Choose either Indo-Pacific
#'ReefBudget ("dbase_type = "IPRB") or U.S. Pacific Islands rates developed
#'by Tye Kindinger, NOAA PIFSC ("dbase_type = "Kindinger").
#'@param data_belt fish belt survey data.
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
#'fish_belt <- process_fish(data_spc = fish_data_spc, 
#'dbase_type = "Kindinger", data_belt = fish_data_belt)

process_fish <- function(data_spc,
                         dbase_type = c("IPRB", "Kindinger"),
                         data_belt = TRUE) {
  
  ifelse(dbase_type == "IPRB", rates_dbase <- fish_erosion_dbase_iprb, rates_dbase <- fish_erosion_dbase_kindinger)
  

  fish_belt <- calc_fish_belt(
                data = data_belt, 
                dbase_type = dbase_type, 
                full_summary = TRUE)
  
  
  fish_fixed_spc <- calc_fish_fixed_spc(
                      data = data_spc, 
                      dbase_type = dbase_type,
                      shape_file = fish_pacific_islands_shapefile)
  
  
  fish_strs_spc <- calc_fish_strs_spc(
                    data = data_spc, 
                    dbase_type = dbase_type,
                    shape_file = fish_pacific_islands_shapefile)
  
  
  if (data_belt == TRUE) {
      return(rbind(fish_belt$fish_erosion_site, 
                   fish_fixed_spc,
                   fish_strs_spc))
  }
  
  
  
  if (data_belt == FALSE) {
      return(rbind(fish_fixed_spc,
                   fish_strs_spc))
  }

  
}


