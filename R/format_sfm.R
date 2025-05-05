#' Format SfM production data for processing
#'
#'@author Hannah Barkley
#'
#'@param data SfM data set, output of `combine_sfm`.
#'
#'@export format_sfm

format_sfm <- function(data) {
  data$LINEAR_METER <- as.numeric(data$Shape_Length)
  data$SUBSTRATE_COVER_CM <- as.numeric(data$SLength) * 100
  data$OBS_ID <- data$OBJECTID
  data$CB_METHOD <- "SfM"
  data$MORPHOLOGYCODE[data$MORPHOLOGYCODE == "<Null>"] <- NA
  data$TRANSECT_OBSID <- data$OBJECTID..
  data$MORPHOLOGYCODE[data$MORPHOLOGYCODE == "EN"] <- "EM"
  data$SUBSTRATE_CODE[data$SUBSTRATE_CODE %in% c("CCAR", "CCAH")] <
    "CCA"
  data$SUBSTRATE_CODE[data$SUBSTRATE_CODE %in% c("TURFR", "TURFH")] <-
    "TURF"
  data$SUBSTRATE_CODE[data$SUBSTRATE_CODE %in% c("PAD")] <-
    "MA"

  data$OCC_SITEID_TRANSECT <-
    paste(data$OCC_SITEID,
          data$CB_TRANSECTID,
          sep = "-")

  transect_summary <-
    data %>%
    dplyr::group_by(.data$OCC_SITEID,
                    .data$CB_TRANSECTID) %>%
    summarize(
      TRANSECT_PLANAR_LENGTH_M = (sum(.data$LINEAR_METER)),
      TRANSECT_TOTAL_SUBSTRATE_COVER_M  =
        sum(.data$SUBSTRATE_COVER_CM / 100, na.rm = TRUE),
      TRANSECT_RUGOSITY =
        .data$TRANSECT_TOTAL_SUBSTRATE_COVER_M /
        .data$TRANSECT_PLANAR_LENGTH_M
    )

  transect_summary$OCC_SITEID_TRANSECT <-
    paste(transect_summary$OCC_SITEID,
          transect_summary$CB_TRANSECTID,
          sep = "-")


  data$TRANSECT_PLANAR_LENGTH_M <-
    transect_summary$TRANSECT_PLANAR_LENGTH_M[match(data$OCC_SITEID_TRANSECT, transect_summary$OCC_SITEID_TRANSECT)]

  data$TRANSECT_PLANAR_LENGTH_M <- round(data$TRANSECT_PLANAR_LENGTH_M, 2)

  data <- data[c(
    "REGION",
    "REGIONCODE",
    "YEAR",
    "CRUISE_ID",
    "LOCATION",
    "LOCATIONCODE",
    "OCC_SITEID",
    "LATITUDE",
    "LONGITUDE",
    "SITE_DEPTH_M",
    "LOCALDATE",
    "CB_TRANSECTID",
    "OCC_SITEID_TRANSECT",
    "TRANSECT_PLANAR_LENGTH_M",
    "SUBSTRATE_CODE",
    "MORPHOLOGYCODE",
    "SUBSTRATE_COVER_CM"
  )]
  return(data)
}
