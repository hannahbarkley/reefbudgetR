#' Calculate carbonate production rates from benthic census data
#'
#'@author Hannah Barkley
#'
#'@param data Benthic field data set.
#'@param dbase_type Production database to use, either Indo-Pacific ReefBudget ("IPRB")
#'or U.S. Pacific Islands NCRMP-specific database ("NCRMP"). The Indo-Pacific ReefBudget
#'database is derived from "IP Calcification and bioerosion rates database v.1.3",
#'downloaded from https://geography.exeter.ac.uk/reefbudget/indopacific/. Defaults to "NCRMP".
#'@param method_name Transect design by which data were collected ("IPRB", "Chords", or "SfM").
#'@param qc_check Remove transects where rugosity < 1.
#'@param sites_metadata Data frame containing information about sites
#'
#'@import dplyr
#'@import tools
#'@import tidyr
#'@import reshape2
#'
#'@export process_prod
#'
#'@examples
#' benthic_data <- read.csv("ESD_CarbBudget_Benthic_OAHU_2021.csv",
#'    na = "", check.names = FALSE)
#'
#' prod_iprb <- process_prod(
#'    data = benthic_data[benthic_data$CB_METHOD == "IPRB", ],
#'    method_name = "IPRB"
#' )

#' prod_chords <- process_prod(
#'    data = benthic_data[benthic_data$CB_METHOD == "Chords", ],
#'    method_name = "Chords"
#' )

#' prod_sfm <- process_prod(
#'    data = benthic_data[benthic_data$CB_METHOD == "SfM", ],
#'    method_name = "SfM"
#' )

process_prod <- function(data,
                         dbase_type = "NCRMP",
                         method_name = "SfM",
                         full_summary = TRUE,
                         label = NULL,
                         qc_check = FALSE,
                         sites_metadata = NULL,
                         ...) {
  options(dplyr.summarise.inform = FALSE,
          scipen = 999)
  
  test = NULL
  
  test <- tryCatch(
    expr = {
      run_calc_prod(data,
                    dbase_type,
                    method_name)
    },
    error = function(e) {
      print(test)
    }
  )
  
  if(is.null(nrow(test)) == FALSE){
    
    print("CHECK DATA")
    print(test)
    return(errors = test)
    
  }
  
  
  
  calc_prod_output <- run_calc_prod(data,
                                    dbase_type,
                                    method_name)
  data <- calc_prod_output$data
  transect_summary <- calc_prod_output$transect_summary
  
  if(qc_check == TRUE){
  bad_transects <- transect_summary$OCC_SITEID_TRANSECT[transect_summary$TRANSECT_RUGOSITY < 1] 
  
  if (length(bad_transects) > 0) {
    data <- data[!data$OCC_SITEID_TRANSECT %in% bad_transects , ]
    transect_summary <-
      transect_summary[!transect_summary$OCC_SITEID_TRANSECT %in% bad_transects , ]
  }
  }
  
  # Calculate cover, planar production, and carbonate production for each
  # SUBSTRATE_CODE on each TRANSECT
  prod_transect_substratecode <-
    summarize_prod(data,
                   transect_summary,
                   dbase_type,
                   summarize_by = "substrate code",
                   level = "transect")
  
  # Calculate cover, planar production, and carbonate production for each
  # SUBSTRATE_CLASS on each TRANSECT
  prod_transect_substrateclass <-
    summarize_prod(data,
                   transect_summary,
                   dbase_type,
                   summarize_by = "substrate class",
                   level = "transect")
  
  # Calculate cover, planar production, and carbonate production for each
  # CORAL_GROUP on each TRANSECT
  prod_transect_coral <-
    summarize_prod(data,
                   transect_summary,
                   dbase_type,
                   summarize_by = "coral group",
                   level = "transect")
  
  # Calculate cover, planar production, and carbonate production
  # on each TRANSECT
  prod_transect <-
    summarize_prod(data,
                   transect_summary,
                   dbase_type,
                   summarize_by = "overall",
                   level = "transect")
  
  # Calculate cover, planar production, and carbonate production for each
  # SUBSTRATE_CODE on each SITE
  prod_site_substratecode <-
    summarize_prod(data,
                   transect_summary,
                   dbase_type,
                   summarize_by = "substrate code",
                   level = "site")
  
  # Calculate cover, planar production, and carbonate production for each
  # SUBSTRATE_CLASS on each SITE
  
  prod_site_substrateclass <-
    summarize_prod(data,
                   transect_summary,
                   dbase_type,
                   summarize_by = "substrate class",
                   level = "site")
  
  # Calculate cover, planar production, and carbonate production for each
  # CORAL_GROUP on each SITE
  prod_site_coral <-
    summarize_prod(data,
                   transect_summary,
                   dbase_type,
                   summarize_by = "coral group",
                   level = "site")
  
  # Calculate cover, planar production, and carbonate production on each SITE
  prod_site <-
    summarize_prod(data,
                   transect_summary,
                   dbase_type,
                   summarize_by = "overall",
                   level = "site")
  
  if (is.null(label) == FALSE) {
    prod_site$LABEL <- label
    prod_site_substrateclass$LABEL <- label
    prod_site_coral$LABEL <- label
    prod_site_substratecode$LABEL <- label
    prod_transect$LABEL <- label
    prod_transect_substrateclass$LABEL <- label
    prod_transect_coral$LABEL <- label
    prod_transect_substratecode$LABEL <- label
    data$LABEL <- label
  }
  
  if (full_summary == TRUE) {
    data <- data[c(
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
      "SUBSTRATE_COVER_CM",
      "COLONY_PROD_G_YR",
      "COLONY_PROD_G_YR_L95",
      "COLONY_PROD_G_YR_U95"
    )]
    
    
    return(
      list(
        summary_site = prod_site,
        summary_site_substrateclass = prod_site_substrateclass,
        summary_site_coral = prod_site_coral,
        summary_site_substratecode = prod_site_substratecode,
        summary_transect = prod_transect,
        summary_transect_substrateclass =
          prod_transect_substrateclass,
        summary_transect_coral = prod_transect_coral,
        summary_transect_substratecode =
          prod_transect_substratecode,
        data = data
      )
    )
  }
  if (full_summary == FALSE) {
    return(summary_site = prod_site)
  }
}
