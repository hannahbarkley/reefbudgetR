#' Process parrotfish erosion rates from belt and spc data
#'
#'@author Rebecca Weible
#'
#'@param spc_data Fish SPC survey data (default is NULL).
#'@param belt_data Fish Belt survey data (default is NULL).
#'@param method_type Type of survey data to process. "Fixed" (fixed SPC only), "StRS" (stratified random SPC), "Mean StRS" (average of nearby stratified random SPC), or "Belt" (belt data only)
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
#'fish_data_belt <- read.csv("CB_FishBelt_alldata.csv", na = "", check.names = FALSE)
#'
#'fish_belt <- process_fish(spc_data = fish_data_spc,
#'dbase_type = "Kindinger", belt_data = fish_data_belt, subset_distance_m = 2000)


process_fish <- function(spc_data = NULL,
                         belt_data = NULL,
                         method_type = c("Fixed", "StRS", "Mean StRS", "Belt"),
                         fixed_metadata = NULL,
                         dbase_type = c("Kindinger", "IPRB"),
                         subset_distance_m = 6000) {
  
  method_type <- match.arg(method_type)
  dbase_type  <- match.arg(dbase_type)
  
  rates_dbase <- if (dbase_type == "Kindinger") {
    if(exists("fish_erosion_dbase_kindinger")) fish_erosion_dbase_kindinger else stop("Database 'fish_erosion_dbase_kindinger' not found.")
  } else {
    if(exists("fish_erosion_dbase_iprb")) fish_erosion_dbase_iprb else stop("Database 'fish_erosion_dbase_iprb' not found.")
  }
  
  result <- switch(method_type,
                   "Fixed" = {
                     if (is.null(spc_data)) stop("spc_data is missing. Please provide spc_data for the 'Fixed' method.")
                     calc_fish_fixed_spc(
                       data = spc_data, 
                       rates_dbase = rates_dbase
                     )
                   },
                   "StRS" = {
                     if (is.null(spc_data)) stop("spc_data is missing. Please provide spc_data for the 'StRS' method.")
                     calc_fish_strs_spc(
                       data = spc_data,
                       rates_dbase = rates_dbase,
                       method_type = "StRS"
                     )
                   },
                   "Mean StRS" = {
                     if (is.null(spc_data)) stop("spc_data is missing. Please provide spc_data for the 'Mean StRS' method.")
                     calc_fish_strs_spc(
                       data = spc_data,
                       fixed_metadata = fixed_metadata,
                       rates_dbase = rates_dbase,
                       subset_distance_m = subset_distance_m,
                       method_type = "Mean StRS"
                     )
                   },
                   "Belt" = {
                     # Fallback to spc_data if a user accidentally passes belt data into the spc_data argument
                     if (is.null(belt_data)) {
                       if (!is.null(spc_data)) {
                         belt_data <- spc_data
                       } else {
                         stop("belt_data is missing. Please provide belt_data for the 'Belt' method.")
                       }
                     }
                     
                     calc_fish_belt(
                       data = belt_data, 
                       rates_dbase = rates_dbase, 
                       full_summary = TRUE
                     )
                   }
  )
  
  return(result)
}
