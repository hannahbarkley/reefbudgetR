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

sfm_files <- list.files(sfm_folder, full.names = TRUE, recursive = T)
benthic_sfm <- NULL

for (i in seq_along(sfm_files)) {
  transect_i <- read.csv(sfm_files[i], colClasses = "character")

  names(transect_i)[names(transect_i) == 'OBJECTID..'] <- 'OBJECTID'
  names(transect_i)[names(transect_i) == 'OID_'] <- 'OBJECTID'
  names(transect_i)[names(transect_i) == 'Shape..'] <- 'Shape'
  names(transect_i)[names(transect_i) == 'SITE'] <- 'OCC_SITEID'
  names(transect_i)[names(transect_i) == 'Site'] <- 'OCC_SITEID'    
  names(transect_i)[names(transect_i) == 'Transect'] <- 'CB_TRANSECTID'
  names(transect_i)[names(transect_i) == 'CB_TRANSECT'] <- 'CB_TRANSECTID'
  names(transect_i)[names(transect_i) == 'Taxon'] <- 'SUBSTRATE_CODE'
  names(transect_i)[names(transect_i) == 'TAXON_ID'] <- 'SUBSTRATE_CODE'
  names(transect_i)[names(transect_i) == 'MORPH_ID'] <- 'MORPHOLOGYCODE'
  names(transect_i)[names(transect_i) == 'Morph'] <- 'MORPHOLOGYCODE'
  names(transect_i)[names(transect_i) == 'Moprh'] <- 'MORPHOLOGYCODE'
  names(transect_i)[names(transect_i) == 'Year'] <- 'YEAR'
  names(transect_i)[names(transect_i) == 'DEPTH_M'] <- 'SITE_DEPTH_M'

  benthic_sfm <- rbind(benthic_sfm, transect_i)
}


benthic_sfm <- benthic_sfm[with(
  benthic_sfm,
  order(
    benthic_sfm$OCC_SITEID,
    benthic_sfm$CB_TRANSECTID,
    benthic_sfm$OBJECTID
  )
), ]

return(benthic_sfm)

}
