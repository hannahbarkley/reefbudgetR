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
  
  data_comb <- combine_sfm(sfm_folder) %>%
    filter(SUBSTRATE_CODE != "UNKN")
  
  if (is.null(data_comb$OCC_SITEID)) {
    site_id <- basename(sfm_folder)
    data_comb$OCC_SITEID <- site_id
  }
  
  if (is.null(data_comb$YEAR)) {
    matched_year <- sites_metadata$YEAR[match(data_comb$OCC_SITEID[1], sites_metadata$OCC_SITEID)]
    data_comb$YEAR <- as.character(matched_year)
  }
  
  meta_subset <- sites_metadata %>%
    select(OCC_SITEID, YEAR, REGION, REGIONCODE, CRUISE_ID, 
           LOCATION, LOCATIONCODE, LATITUDE, LONGITUDE, 
           SITE_DEPTH_M, LOCALDATE, SITEVISITID) %>%
    mutate(YEAR = as.character(YEAR))
  
  # Join metadata to data 
  data_comb <- data_comb %>%
    mutate(YEAR = as.character(YEAR)) %>%
    left_join(meta_subset, by = c("OCC_SITEID", "YEAR"))
  
  # Format data
  data_comb_format <- format_sfm(data_comb)
  
  sub_lookup <- prod_dbase_ncrmp %>%
    select(SUBSTRATE_CODE, SUBSTRATE_CLASS, SUBSTRATE_NAME) %>%
    distinct(SUBSTRATE_CODE, .keep_all = TRUE)
  
  morph_lookup <- prod_dbase_ncrmp %>%
    select(MORPHOLOGYCODE, MORPHOLOGY) %>%
    distinct(MORPHOLOGYCODE, .keep_all = TRUE)
  
  data_export <- data_comb_format %>%
    # Join Substrate info
    left_join(sub_lookup, by = "SUBSTRATE_CODE") %>%
    # Join Morphology info
    left_join(morph_lookup, by = "MORPHOLOGYCODE") %>%
    mutate(
      CB_METHOD = "SfM",
      # Apply Coral-Only Logic
      MORPHOLOGY = if_else(SUBSTRATE_CLASS != "CORAL", NA_character_, MORPHOLOGY),
      MORPHOLOGYCODE = if_else(SUBSTRATE_CLASS != "CORAL", NA_character_, MORPHOLOGYCODE)
    ) %>%
    # 6. Final Selection
    select(
      REGION, REGIONCODE, YEAR, CRUISE_ID, LOCATION, LOCATIONCODE,
      OCC_SITEID, SITEVISITID, LATITUDE, LONGITUDE, SITE_DEPTH_M,
      LOCALDATE, CB_METHOD, CB_TRANSECTID, OCC_SITEID_TRANSECT,
      TRANSECT_PLANAR_LENGTH_M, SUBSTRATE_CLASS, SUBSTRATE_NAME,
      SUBSTRATE_CODE, MORPHOLOGY, MORPHOLOGYCODE, SUBSTRATE_COVER_CM
    )
  
  return(data_export)
}
