#' Format net production data for NCEI submission
#'
#'@author Hannah Barkley
#'
#'@param data output of `process_net`.
#'
#'@import dplyr
#'@import tools
#'@import tidyr
#'@import reshape2
#'
#'@export format_4ncei
#'

format_4ncei <- function(data) {
  
  # Define NCEI Column Schema
  cols_ncei <- c(
    # Metadata
    "REGION", "REGIONCODE", "YEAR", "CRUISE_ID", "LOCATION", "LOCATIONCODE",
    "OCC_SITEID", "SITEVISITID", "LATITUDE", "LONGITUDE", "SITE_DEPTH_M", "LOCALDATE",
    
    # Benthic Metrics (Mean + SE)
    "RUGOSITY_MEAN", "RUGOSITY_SE",
    "HARD_CORAL_COVER_PCT_MEAN", "HARD_CORAL_COVER_PCT_SE",
    "CCA_COVER_PCT_MEAN", "CCA_COVER_PCT_SE",
    
    # Production Metrics
    "HARD_CORAL_CARB_PROD_KG_M2_YR_MEAN", "HARD_CORAL_CARB_PROD_KG_M2_YR_SE",
    "CCA_CARB_PROD_KG_M2_YR_MEAN", "CCA_CARB_PROD_KG_M2_YR_SE",
    "GROSS_CARB_PROD_KG_M2_YR_MEAN", "GROSS_CARB_PROD_KG_M2_YR_SE",
    
    # Bioerosion Metrics
    "MACROBIOEROSION_KG_M2_YR_MEAN", "MACROBIOEROSION_KG_M2_YR_SE",
    "MICROBIOEROSION_KG_M2_YR_MEAN", "MICROBIOEROSION_KG_M2_YR_SE",
    "BIOEROSION_KG_M2_YR_MEAN", "BIOEROSION_KG_M2_YR_SE",
    "URCHIN_DENSITY_NO_M2_MEAN", "URCHIN_DENSITY_NO_M2_SE",
    "URCHIN_EROSION_KG_M2_YR_MEAN", "URCHIN_EROSION_KG_M2_YR_SE",
    "FISH_BIOMASS_KG_HA_ALL_MEAN", "FISH_BIOMASS_KG_HA_ALL_SE",
    "FISH_DENSITY_ABUNDANCE_HA_ALL_MEAN", "FISH_DENSITY_ABUNDANCE_HA_ALL_SE",
    "FISH_EROSION_KG_M2_YR_ALL_MEAN", "FISH_EROSION_KG_M2_YR_ALL_SE",
    
    # Net Budget Metrics
    "GROSS_EROSION_KG_M2_YR_MEAN", "GROSS_EROSION_KG_M2_YR_SE",
    "NET_CARB_PROD_KG_M2_YR_MEAN", "NET_CARB_PROD_KG_M2_YR_SE"
  )
  
  # Identify which NCEI columns are missing from the input data
  missing_cols <- setdiff(cols_ncei, names(data))
  
  # Bulk create missing columns as NA
  if (length(missing_cols) > 0) {
    data[missing_cols] <- NA
  }
}
