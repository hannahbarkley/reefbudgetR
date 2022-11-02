#' Format SfM production data for processing
#'
#'@author Hannah Barkley
#'
#'@param data SfM data set, output of `combine_sfm`.
#'
#'@export format_sfm

format_sfm <- function(data) {
  data$OCC_SITEID <- data$Site
  data$CB_TRANSECTID <- data$Transect
  data$SUBSTRATE_CODE <- data$Taxon
  data$MORPHOLOGYCODE <- data$Morph
  data$LINEAR_METER <- as.numeric(data$Shape_Length)
  data$SUBSTRATE_COVER_CM <- as.numeric(data$SLength) * 100
  data$OBS_ID <- data$OBJECTID
  data$CB_METHOD <- "SfM"
  data$MORPHOLOGYCODE[data$MORPHOLOGYCODE == "<Null>"] <- NA
  data$TRANSECT_OBSID <- data$OBJECTID..
  data$MORPHOLOGYCODE[data$MORPHOLOGYCODE == "EN"] <- "EM"
  data$SUBSTRATE_CODE[data$SUBSTRATE_CODE %in% c("CCAR", "CCAH")] <-
    "CCA"
  data$SUBSTRATE_CODE[data$SUBSTRATE_CODE %in% c("TURFR", "TURFH")] <-
    "TURF"
  data <- data[c(
    "CRUISEID",
    "LOCALDATE",
    "REGIONCODE",
    "LOCATIONCODE",
    "OCC_SITEID",
    "CB_METHOD",
    "CB_TRANSECTID",
    "LINEAR_METER",
    "SUBSTRATE_CODE",
    "MORPHOLOGYCODE",
    "SUBSTRATE_COVER_CM"
  )]
  return(data)
}
