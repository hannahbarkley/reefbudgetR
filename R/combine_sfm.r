#' Combine individual SfM benthic data transect files
#'
#'@author Hannah Barkley
#'
#'@param sfm_folder Directory containing SfM transect files.
#'@param region_codeSurvey region ("MHI", "MARIAN", ...).
#'@param cruise_id Cruise ID for survey collection ("MP2108", "RA2201", ...).
#'
#'@import lubridate
#'@export combine_sfm
#'
#'@examples
#'combine_sfm(
#'    sfm_folder = "...",
#'    region_code = "MHI",
#'    local_date = "2020-01-01",
#'    cruise_id = "MP2108",
#'    surveyor = "HCB"
#' )

combine_sfm <- function(sfm_folder) {
  
  # Get File List
  sfm_files <- list.files(sfm_folder, full.names = TRUE, recursive = TRUE)
  
  # Define Renaming Map 
  
  process_file <- function(file) {
    df <- read.csv(file, colClasses = "character")
    
    current_names <- names(df)
    
    new_names <- dplyr::case_match(current_names,
                                   c("OBJECTID..", "OID_") ~ "OBJECTID",
                                   "Shape.."               ~ "Shape",
                                   c("SITE", "Site")       ~ "OCC_SITEID",
                                   c("Transect", "CB_TRANSECT") ~ "CB_TRANSECTID",
                                   c("Taxon", "TAXON_ID")  ~ "SUBSTRATE_CODE",
                                   c("MORPH_ID", "Morph", "Moprh") ~ "MORPHOLOGYCODE",
                                   "Year"                  ~ "YEAR",
                                   "DEPTH_M"               ~ "SITE_DEPTH_M",
                                   .default = current_names # Keep original name if no match found
    )
    
    names(df) <- new_names
    
    # Ensure OBJECTID exists 
    if (!"OBJECTID" %in% names(df)) {
      df$OBJECTID <- NA
    }
    
    return(df)
  }
  
  # Read and Combine 
  benthic_sfm <- purrr::map_dfr(sfm_files, process_file)
  
  # Sort Data
  benthic_sfm <- benthic_sfm %>%
    dplyr::arrange(OCC_SITEID, CB_TRANSECTID, OBJECTID)
  
  return(benthic_sfm)
}
