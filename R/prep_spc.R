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
  
  data_spc <- data_spc %>%
    dplyr::rename(any_of(c(
      "CRUISE_ID"  = "MISSIONID",
      "REGION"     = "REGION_NAME",
      "LOCATION"   = "ISLAND",
      "REA_SITEID" = "SITE",
      "LOCALDATE"  = "DATE_",
      "YEAR"       = "OBS_YEAR",
      "METHOD"     = "CB_METHOD",
      "SPECIES"    = "TAXON_CODE",
      "TAXONNAME"  = "TAXON_NAME"
    )))
  
  data_spc$LOCALDATE <- lubridate::parse_date_time(
    data_spc$LOCALDATE, 
    orders = c("ymd", "dmy", "mdy"), 
    quiet = TRUE
  ) %>% as_date()
  
  if (!is.null(fixed_metadata)) {

    meta_lookup <- fixed_metadata %>%
      dplyr::select(REGION, REGIONCODE, LOCATION, LOCATIONCODE) %>%
      distinct() 
    
    data_spc <- data_spc %>%
      left_join(meta_lookup, by = c("REGION", "LOCATION"))
    
  } else {

    data_spc$REGIONCODE <- NA_character_
    data_spc$LOCATIONCODE <- NA_character_
  }
  
  data_spc <- data_spc %>%
    relocate(REGIONCODE, .after = REGION) %>%
    relocate(YEAR, .after = REGIONCODE) %>%
    relocate(LOCATIONCODE, .after = LOCATION)
  
  return(data_spc)
}
