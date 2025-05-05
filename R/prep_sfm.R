#' Prepare SfM production data for processing
#'
#'@author Hannah Barkley
#'
#'@param sfm_folder Directory containing SfM transect files.
#'@param sites_metadata Data frame containing information about sites
#'
#'@export prep_sfm

prep_sfm <- function(sfm_folder, sites_metadata) {
  
  require(abjutils)
  require(reefbudgetR)
  require(tidyverse)
  
  data_comb <-
    combine_sfm(sfm_folder) %>% 
    filter(SUBSTRATE_CODE != "UNKN") 
  
  site_id <- file_sans_ext(sfm_folder) 
  data_comb$REGIONCODE <- sites_metadata$REGIONCODE[match(site_id, sites_metadata$OCC_SITEID)]
  data_comb$REGION <- sites_metadata$REGION[match(site_id, sites_metadata$OCC_SITEID)]
  data_comb$CRUISE_ID <- sites_metadata$CRUISE_ID[match(site_id, sites_metadata$OCC_SITEID)]
  data_comb$YEAR <- sites_metadata$YEAR[match(site_id, sites_metadata$OCC_SITEID)]
  data_comb$LOCATION <- sites_metadata$LOCATION[match(site_id, sites_metadata$OCC_SITEID)]
  data_comb$LOCATIONCODE <- sites_metadata$LOCATIONCODE[match(site_id, sites_metadata$OCC_SITEID)]
  data_comb$LATITUDE <- sites_metadata$LATITUDE[match(site_id, sites_metadata$OCC_SITEID)]
  data_comb$LONGITUDE <- sites_metadata$LONGITUDE[match(site_id, sites_metadata$OCC_SITEID)]
  data_comb$SITE_DEPTH_M <- sites_metadata$SITE_DEPTH_M[match(site_id, sites_metadata$OCC_SITEID)]
  data_comb$LOCALDATE <- sites_metadata$LOCALDATE[match(site_id, sites_metadata$OCC_SITEID)]
  
  data_comb_format <- format_sfm(data_comb)
  
  data_comb_format$SUBSTRATE_CLASS <-
    prod_dbase_ncrmp$SUBSTRATE_CLASS[match(data_comb_format$SUBSTRATE_CODE, prod_dbase_ncrmp$SUBSTRATE_CODE)]
  data_comb_format$SUBSTRATE_NAME <-
    prod_dbase_ncrmp$SUBSTRATE_NAME[match(data_comb_format$SUBSTRATE_CODE, prod_dbase_ncrmp$SUBSTRATE_CODE)]
  data_comb_format$MORPHOLOGY <-
    prod_dbase_ncrmp$MORPHOLOGY[match(data_comb_format$MORPHOLOGYCODE, prod_dbase_ncrmp$MORPHOLOGYCODE)]

  data_comb_format$CB_METHOD <- "SfM"
  
  data_comb_format$MORPHOLOGY[data_comb_format$SUBSTRATE_CLASS != "CORAL"] <- NA
  data_comb_format$MORPHOLOGYCODE[data_comb_format$SUBSTRATE_CLASS != "CORAL"] <- NA
  
  data_export <- data_comb_format[c(
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
  
  
  return(data_export)
  
}
