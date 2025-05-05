#' Combine and format individual SfM benthic data transect files
#'
#'@author Hannah Barkley
#'
#'@param site_file CSV file with site info.
#'@param sfm_folder Directory containing SfM transect files.
#'
#'@import lubridate tidyr
#'@export format_benthic_sfm

format_benthic_sfm <- function(site_file, sfm_folder){
  
  site_info <- read.csv(site_file)
  
  benthic_sfm_data <-
    combine_sfm(sfm_folder) %>%
    format_sfm()
  
  benthic_sfm_data$TRANSECT_PLANAR_LENGTH_M <-
    round(benthic_sfm_data$TRANSECT_PLANAR_LENGTH_M, 2)
  
  benthic_sfm <-
    right_join(site_info, benthic_sfm_data, by = "OCC_SITEID")
  
  benthic_sfm$SUBSTRATE_CLASS <-
    prod_dbase_ncrmp$SUBSTRATE_CLASS[match(benthic_sfm$SUBSTRATE_CODE, prod_dbase_ncrmp$SUBSTRATE_CODE)]
  
  benthic_sfm$SUBSTRATE_NAME <-
    prod_dbase_ncrmp$SUBSTRATE_NAME[match(benthic_sfm$SUBSTRATE_CODE, prod_dbase_ncrmp$SUBSTRATE_CODE)]
  
  benthic_sfm$MORPHOLOGY <-
    prod_dbase_ncrmp$MORPHOLOGY[match(benthic_sfm$MORPHOLOGYCODE, prod_dbase_ncrmp$MORPHOLOGYCODE)]
  
  benthic_sfm_export <- benthic_sfm[c(
    "REGIONCODE",
    "YEAR",
    "LOCATION",
    "OCC_SITEID",
    "LATITUDE",
    "LONGITUDE",
    "SITE_DEPTH_M",
    "CB_METHOD",
    "CB_TRANSECTID",
    "OCC_SITEID_TRANSECT",
    "TRANSECT_PLANAR_LENGTH_M",
    "SUBSTRATE_CLASS",
    "SUBSTRATE_NAME",
    "SUBSTRATE_CODE",
    "MORPHOLOGY",
    "MORPHOLOGYCODE",
    "SUBSTRATE_COVER_CM"
  )]
  
}

