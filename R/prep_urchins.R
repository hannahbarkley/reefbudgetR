#' Prep urchin census data
#'
#'@author Hannah Barkley
#'
#'@param data Urchin observation data set.
#'@param sites_metadata Data frame containing information about sites
#'
#'
#'@export prep_urchins

prep_urchins <- function(data, sites_metadata) {
  
  if (is.null(data$CBURCHINID) == TRUE) {
    data$CBURCHINID <- NA
  }
  
  data$SITEVISITID[is.null(data$SITEVISITID) == TRUE] <- NA
  data$TAXON_NAME[is.null(data$TAXON_NAME) == TRUE] <- NA
  data$TAXON_CODE[is.null(data$TAXON_CODE) == TRUE] <- NA
  
  data$TEST_SIZE_BIN_161_180_MM[is.null(data$TEST_SIZE_BIN_161_180_MM) == TRUE] <- 0
  data$TEST_SIZE_BIN_181_200_MM[is.null(data$TEST_SIZE_BIN_181_200_MM) == TRUE] <- 0

  data <- data %>% 
    mutate(
      across(TEST_SIZE_BIN_0_20_MM:TEST_SIZE_BIN_181_200_MM, ~replace_na(.x, 0))
    )
  
  names(data)[names(data) == "DEPTH_M"] <- "SITE_DEPTH_M"
  
  data$YEAR <- as.factor(data$YEAR)
  
  data <- data[c(
    "SITEVISITID",
    "CBURCHINID",
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
    "TRANSECT_LENGTH_M",
    "URCH_OBS_TF",
    "TAXON_NAME",
    "TAXON_CODE",
    "TEST_SIZE_BIN_0_20_MM",
    "TEST_SIZE_BIN_21_40_MM",
    "TEST_SIZE_BIN_41_60_MM",
    "TEST_SIZE_BIN_61_80_MM",
    "TEST_SIZE_BIN_81_100_MM",
    "TEST_SIZE_BIN_101_120_MM",
    "TEST_SIZE_BIN_121_140_MM",
    "TEST_SIZE_BIN_141_160_MM",
    "TEST_SIZE_BIN_161_180_MM",
    "TEST_SIZE_BIN_181_200_MM"
  )]
  
  return(data)
}
