#' Format SfM production data for processing
#'
#'@author Hannah Barkley
#'
#'@param data SfM data set, output of `combine_sfm`.
#'
#'@export format_sfm

format_sfm <- function(data) {
  
  # Standardize Numeric Columns
  data$LINEAR_METER <- as.numeric(data$Shape_Length)
  data$SUBSTRATE_COVER_CM <- as.numeric(data$SLength) * 100
  data$OBS_ID <- data$OBJECTID
  data$CB_METHOD <- "SfM"
  data$TRANSECT_OBSID <- data$OBJECTID..
  
  # Clean Factors and Substrate Codes 
  data$MORPHOLOGYCODE[data$MORPHOLOGYCODE == "<Null>"] <- NA
  data$MORPHOLOGYCODE[data$MORPHOLOGYCODE == "EN"] <- "EM"

  data$SUBSTRATE_CODE <- dplyr::recode(data$SUBSTRATE_CODE,
                                       "CCAR"  = "CCA",
                                       "CCAH"  = "CCA",
                                       "TURFR" = "TURF",
                                       "TURFH" = "TURF",
                                       "PAD"   = "MA"
  )

  data$OCC_SITEID_TRANSECT <- paste0(data$OCC_SITEID, "-", data$CB_TRANSECTID)
  
  # Calculate Transect Metrics
  data <- data %>%
    group_by(OCC_SITEID_TRANSECT) %>%
    mutate(
      TRANSECT_PLANAR_LENGTH_M = round(sum(LINEAR_METER, na.rm = TRUE), 2)
    ) %>%
    ungroup()
  
  data <- data[c(
    "REGION", "REGIONCODE", "YEAR", "CRUISE_ID", "LOCATION", "LOCATIONCODE",
    "OCC_SITEID", "SITEVISITID", "LATITUDE", "LONGITUDE", "SITE_DEPTH_M",
    "LOCALDATE", "CB_TRANSECTID", "OCC_SITEID_TRANSECT", 
    "TRANSECT_PLANAR_LENGTH_M", "SUBSTRATE_CODE", "MORPHOLOGYCODE", 
    "SUBSTRATE_COVER_CM"
  )]
  
  return(data)
}
