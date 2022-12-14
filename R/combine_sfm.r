#' Combine individual SfM benthic data transect files
#'
#'@author Hannah Barkley
#'
#'@param sfm_folder Directory containing SfM transect files.
#'@param region_codeSurvey region ("MHI", "MARIAN", ...).
#'@param location_code Location code where surveys were conducted
#'("OAH", "GUA", "SAI", ...).
#'@param local_date Local Date of imagery collection (yyyy-mm-dd).
#'@param cruise_id Cruise ID for survey collection ("MP2108", "RA2201", ...).
#'@param surveyor Initials of data collector.
#'
#'@import lubridate
#'@export combine_sfm
#'
#'@examples
#'combine_sfm(
#'    sfm_folder = "...",
#'    region_code = "MHI",
#'    location_code = "OAH",
#'    local_date = "2020-01-01",
#'    cruise_id = "MP2108",
#'    surveyor = "HCB"
#' )

combine_sfm <- function(sfm_folder,
                        region_code,
                        location_code,
                        local_date,
                        cruise_id,
                        surveyor
                       ) {

sfm_files <- list.files(sfm_folder, full.names = TRUE)
benthic_sfm <- NULL

for (i in seq_along(sfm_files)) {
  transect_i <- read.csv(sfm_files[i], colClasses = "character")
  benthic_sfm <- rbind(benthic_sfm, transect_i)
}

benthic_sfm$REGIONCODE <- region_code
benthic_sfm$LOCATIONCODE <- location_code
benthic_sfm$LOCALDATE <- lubridate::ymd(local_date)
benthic_sfm$CRUISEID <- cruise_id
benthic_sfm$SURVEYOR <- surveyor

benthic_sfm <- benthic_sfm[with(
  benthic_sfm,
  order(
    benthic_sfm$Site,
    benthic_sfm$Transect,
    benthic_sfm$OBJECTID..
  )
), ]

return(benthic_sfm)

}
