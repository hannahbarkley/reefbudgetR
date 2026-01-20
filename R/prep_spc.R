#' Reformat NCRMP stratified random SPC data 
#'
#'@author Hannah Barkley
#'
#'@param data_spc Fish spc survey data.
#'@param fixed_metadata Data frame with fixed site metadata
#'
#'
#'@import dplyr
#'@import tidyr

#'
#'@export prep_spc

prep_spc <- function(data_spc, fixed_metadata = NULL) {
  
  names(data_spc)[names(data_spc) == "MISSIONID"] <- "CRUISE_ID"
  names(data_spc)[names(data_spc) == "REGION_NAME"] <- "REGION"
  names(data_spc)[names(data_spc) == "ISLAND"] <- "LOCATION"
  names(data_spc)[names(data_spc) == "SITE"] <- "REA_SITEID"
  names(data_spc)[names(data_spc) == "DATE_"] <- "LOCALDATE"
  names(data_spc)[names(data_spc) == "REGION_NAME"] <- "REGION"
  
  names(data_spc)[names(data_spc) == "ISLAND"] <- "LOCATION"
  names(data_spc)[names(data_spc) == "OBS_YEAR"] <- "YEAR"
  names(data_spc)[names(data_spc) == "DATE_"] <- "LOCALDATE"
  
  names(data_spc)[names(data_spc) == "CB_METHOD"] <- "METHOD"
  names(data_spc)[names(data_spc) == "TAXON_CODE"] <- "SPECIES"
  names(data_spc)[names(data_spc) == "TAXON_NAME"] <- "TAXONNAME"
  
  
  if (is.na(ymd(data_spc$LOCALDATE)[1]) == TRUE) {
    data_spc$LOCALDATE <- dmy(data_spc$LOCALDATE)
  } else {
    data_spc$LOCALDATE <- ymd(data_spc$LOCALDATE)
  }
  
  
  if (is.null(fixed_metadata) == FALSE) {
    data_spc$REGIONCODE <- fixed_metadata$REGIONCODE[match(data_spc$REGION, fixed_metadata$REGION)]
    data_spc$LOCATIONCODE <- fixed_metadata$LOCATIONCODE[match(data_spc$LOCATION, fixed_metadata$LOCATION)]
  } else{
    data_spc$REGIONCODE <- NA
    data_spc$LOCATIONCODE <- NA
    
  }
  
  data_spc <- data_spc %>% relocate(REGIONCODE, .after = REGION) %>% relocate(YEAR, .after = REGIONCODE)  %>% relocate(LOCATIONCODE, .after = LOCATION)
  
  return(data_spc)
  
  
}
