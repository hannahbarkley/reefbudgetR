#' Calculate carbonate production rates from benthic census data
#'
#'@author Hannah Barkley
#'
#'@param data Benthic field data set.
#'@param dbase_type Production database to use, either Indo-Pacific ReefBudget ("IPRB"), 
#'the U.S. Pacific Islands NCRMP-specific database ("NCRMP"), or a customized location-specific database ("Custom"). The Indo-Pacific ReefBudget
#'database is derived from "IP Calcification and bioerosion rates database v.1.3",
#'downloaded from https://geography.exeter.ac.uk/reefbudget/indopacific/. Defaults to "NCRMP".
#'@param method_name Transect design by which data were collected ("IPRB", "Chords", or "SfM").
#'@param qc_check Remove transects where rugosity < 1.
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
#'    data = benthic_data,
#' )


process_prod <- function(data,
                         dbase_type = "NCRMP",
                         macro_rates = "IPRB",
                         micro_rates = "IPRB",
                         full_summary = TRUE,
                         prod_dbase_custom = NULL,
                         qc_check = FALSE,
                         ...) {
  
  options(dplyr.summarise.inform = FALSE, scipen = 999)
  
  # Calculate production rates for each data row
  calc_prod_output <- tryCatch({
    run_calc_prod(data, dbase_type, prod_dbase_custom)
  }, error = function(e) {
    message("Critical failure in run_calc_prod: ", e$message)
    return(NULL)
  })
  
  # Check for Failures
  if (is.null(calc_prod_output)) {
    stop("Process halted: run_calc_prod could not execute.")
  }
  
  # If row-level errors occurred in the database lookup, report them
  if (nrow(calc_prod_output$error_log) > 0) {
    warning(paste(nrow(calc_prod_output$error_log), "rows failed. See error_log in output."))
  }
  
  data <- calc_prod_output$data
  transect_summary <- calc_prod_output$transect_summary
  error_log <- calc_prod_output$error_log
  
  # Quality Control (Rugosity Check)
  if (qc_check == TRUE) {
    # Flag transects where rugosity is physically impossible (< 1)
    bad_transects <- transect_summary$OCC_SITEID_TRANSECT[transect_summary$TRANSECT_RUGOSITY < 1] 
    
    if (length(bad_transects) > 0) {
      data <- data[!data$OCC_SITEID_TRANSECT %in% bad_transects, ]
      transect_summary <- transect_summary[!transect_summary$OCC_SITEID_TRANSECT %in% bad_transects, ]
    }
  }
  
  # Summarize at different levels
  
  run_sum <- function(s_by, lvl) {
    summarize_prod(
      data, 
      transect_summary, 
      dbase_type,
      prod_dbase_custom = prod_dbase_custom, # Pass the object through
      summarize_by = s_by, 
      level = lvl, 
      macro_rates = macro_rates, 
      micro_rates = micro_rates
    )
  }
  
  # TRANSECT LEVEL
  prod_transect_substratecode  <- run_sum("substrate code", "transect")
  prod_transect_substrateclass <- run_sum("substrate class", "transect")
  prod_transect_coral          <- run_sum("coral group", "transect")
  prod_transect                <- run_sum("overall", "transect")
  
  # SITE LEVEL
  prod_site_substratecode  <- run_sum("substrate code", "site")
  prod_site_substrateclass <- run_sum("substrate class", "site")
  prod_site_coral          <- run_sum("coral group", "site")
  prod_site                <- run_sum("overall", "site")
  
  # 5. Output Packaging
  if (full_summary == TRUE) {
    # Keep only essential columns for the final observation table
    obs_cols <- c("REGION", "REGIONCODE", "YEAR", "CRUISE_ID", "LOCATION", "LOCATIONCODE", 
                  "OCC_SITEID", "SITEVISITID", "LATITUDE", "LONGITUDE", "SITE_DEPTH_M", 
                  "LOCALDATE", "CB_METHOD", "CB_TRANSECTID", "OCC_SITEID_TRANSECT", 
                  "TRANSECT_PLANAR_LENGTH_M", "SUBSTRATE_CLASS", "SUBSTRATE_NAME", 
                  "SUBSTRATE_CODE", "MORPHOLOGY", "MORPHOLOGYCODE", "SUBSTRATE_COVER_CM", 
                  "COLONY_PROD_G_YR", "COLONY_PROD_G_YR_L95", "COLONY_PROD_G_YR_U95")
    
    data <- data[, intersect(names(data), obs_cols)]
    
    # Generate Site-Level Metadata
    sites_metadata <- transect_summary %>%
      group_by(OCC_SITEID) %>%
      summarize(across(.cols = any_of(c("REGION", "REGIONCODE", "YEAR", "CRUISE_ID", 
                                        "LOCALDATE", "LOCATION", "LOCATIONCODE", 
                                        "LATITUDE", "LONGITUDE", "SITE_DEPTH_M", "SITEVISITID")), 
                       .fns = ~first(.x)))
    
    return(list(
      summary_site = prod_site,
      summary_site_substrateclass = prod_site_substrateclass,
      summary_site_coral = prod_site_coral,
      summary_site_substratecode = prod_site_substratecode,
      summary_transect = prod_transect,
      summary_transect_substrateclass = prod_transect_substrateclass,
      summary_transect_coral = prod_transect_coral,
      summary_transect_substratecode = prod_transect_substratecode,
      data = data,
      transect_summary = transect_summary,
      sites_metadata = sites_metadata,
      error_log = error_log # Crucial for auditing
    ))
  } else {
    return(prod_site)
  }
}
