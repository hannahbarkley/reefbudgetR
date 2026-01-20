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
  
  
  if (is.null(unique(data_comb$OCC_SITEID)) == FALSE) {
    site_i <- unique(data_comb$OCC_SITEID)
  } else {
    site_i <-  file_sans_ext(strsplit(sfm_folder_i, "_")[[1]][2])
  }
  
  if (is.null(data_comb$YEAR) == TRUE) {
    data_comb$YEAR <- unique(as.character(sites_metadata$YEAR[match(site_i, sites_metadata$OCC_SITEID)]))
  }
  
  year_i <- unique(data_comb$YEAR)
  
  site_id <- file_sans_ext(sfm_folder)
  data_comb$REGIONCODE <- sites_metadata$REGIONCODE[match(
    paste0(site_i, "_", year_i),
    paste0(sites_metadata$OCC_SITEID, "_", sites_metadata$YEAR)
  )]
  data_comb$REGION <- sites_metadata$REGION[match(
    paste0(site_i, "_", year_i),
    paste0(sites_metadata$OCC_SITEID, "_", sites_metadata$YEAR)
  )]
  data_comb$CRUISE_ID <- sites_metadata$CRUISE_ID[match(
    paste0(site_i, "_", year_i),
    paste0(sites_metadata$OCC_SITEID, "_", sites_metadata$YEAR)
  )]
  data_comb$LOCATION <- sites_metadata$LOCATION[match(
    paste0(site_i, "_", year_i),
    paste0(sites_metadata$OCC_SITEID, "_", sites_metadata$YEAR)
  )]
  data_comb$LOCATIONCODE <- sites_metadata$LOCATIONCODE[match(
    paste0(site_i, "_", year_i),
    paste0(sites_metadata$OCC_SITEID, "_", sites_metadata$YEAR)
  )]
  data_comb$LATITUDE <- sites_metadata$LATITUDE[match(
    paste0(site_i, "_", year_i),
    paste0(sites_metadata$OCC_SITEID, "_", sites_metadata$YEAR)
  )]
  data_comb$LONGITUDE <- sites_metadata$LONGITUDE[match(
    paste0(site_i, "_", year_i),
    paste0(sites_metadata$OCC_SITEID, "_", sites_metadata$YEAR)
  )]
  data_comb$SITE_DEPTH_M <- sites_metadata$SITE_DEPTH_M[match(
    paste0(site_i, "_", year_i),
    paste0(sites_metadata$OCC_SITEID, "_", sites_metadata$YEAR)
  )]
  data_comb$LOCALDATE <- sites_metadata$LOCALDATE[match(
    paste0(site_i, "_", year_i),
    paste0(sites_metadata$OCC_SITEID, "_", sites_metadata$YEAR)
  )]
  data_comb$SITEVISITID <- sites_metadata$SITEVISITID[match(
    paste0(site_i, "_", year_i),
    paste0(sites_metadata$OCC_SITEID, "_", sites_metadata$YEAR)
  )]
  
  
  
  data_comb_format <- format_sfm(data_comb)
  
  data_comb_format$SUBSTRATE_CLASS <-
    prod_dbase_ncrmp$SUBSTRATE_CLASS[match(data_comb_format$SUBSTRATE_CODE,
                                           prod_dbase_ncrmp$SUBSTRATE_CODE)]
  data_comb_format$SUBSTRATE_NAME <-
    prod_dbase_ncrmp$SUBSTRATE_NAME[match(data_comb_format$SUBSTRATE_CODE,
                                          prod_dbase_ncrmp$SUBSTRATE_CODE)]
  data_comb_format$MORPHOLOGY <-
    prod_dbase_ncrmp$MORPHOLOGY[match(data_comb_format$MORPHOLOGYCODE,
                                      prod_dbase_ncrmp$MORPHOLOGYCODE)]
  
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
    "SITEVISITID",
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
