#' Prep urchin census data
#'
#'@author Hannah Barkley
#'
#'@param data Urchin observation data set.
#'@param sites_metadata Data frame containing information about sites
#'
#'
#'@export prep_urchins
#'
prep_urchins <- function(data, sites_metadata) {
  
  # Standardize Depth Column Name
  if ("DEPTH_M" %in% names(data)) {
    names(data)[names(data) == "DEPTH_M"] <- "SITE_DEPTH_M"
  }
  
  #  Define Expected Columns
  # Columns to fill with NA if missing
  cols_na_default <- c("SITEVISITID", "CBURCHINID", "TAXON_NAME", "TAXON_CODE", "URCH_OBS_TF")
  
  # Columns to fill with 0 if missing 
  cols_zero_default <- c(
    "TEST_SIZE_BIN_0_20_MM", "TEST_SIZE_BIN_21_40_MM", "TEST_SIZE_BIN_41_60_MM",
    "TEST_SIZE_BIN_61_80_MM", "TEST_SIZE_BIN_81_100_MM", "TEST_SIZE_BIN_101_120_MM",
    "TEST_SIZE_BIN_121_140_MM", "TEST_SIZE_BIN_141_160_MM", "TEST_SIZE_BIN_161_180_MM",
    "TEST_SIZE_BIN_181_200_MM"
  )
  
  # Bulk Create Missing Columns 
  missing_na_cols <- setdiff(cols_na_default, names(data))
  if (length(missing_na_cols) > 0) {
    data[missing_na_cols] <- NA
  }
  
  missing_zero_cols <- setdiff(cols_zero_default, names(data))
  if (length(missing_zero_cols) > 0) {
    data[missing_zero_cols] <- 0
  }
  
  # Clean Data & Reorder
  data <- data %>%
    mutate(
      across(all_of(cols_zero_default), ~ replace_na(as.numeric(.x), 0)),
      YEAR = as.factor(YEAR)
    ) %>%
    select(
      SITEVISITID, CBURCHINID, REGION, REGIONCODE, YEAR, CRUISE_ID,
      LOCATION, LOCATIONCODE, OCC_SITEID, LATITUDE, LONGITUDE,
      SITE_DEPTH_M, LOCALDATE, CB_TRANSECTID, TRANSECT_LENGTH_M,
      URCH_OBS_TF, TAXON_NAME, TAXON_CODE, 
      all_of(cols_zero_default)
    )
  
  return(data)
}
