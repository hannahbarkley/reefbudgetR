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

combine_sfm <- function(sfm_folder,
                        region_code,
                        cruise_id
                       ) {

sfm_files <- list.files(sfm_folder, full.names = TRUE)
benthic_sfm <- NULL

for (i in seq_along(sfm_files)) {
  transect_i <- read.csv(sfm_files[i], colClasses = "character")
  colnames(transect_i) <- c("OBJECTID", "Shape", "OCC_SITEID", "CB_TRANSECTID","SUBSTRATE_CODE", "MORPHOLOGYCODE","Shape_Length", "SLength")
  benthic_sfm <- rbind(benthic_sfm, transect_i)

}

benthic_sfm$REGIONCODE <- region_code
benthic_sfm$LOCATIONCODE <- substring(benthic_sfm$OCC_SITEID,5,7)
benthic_sfm$LOCALDATE <- NA
benthic_sfm$CRUISE_ID <- cruise_id

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
