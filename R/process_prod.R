#' Calculate carbonate production rates from benthic census data
#'
#'@author Hannah Barkley
#'
#'@param data Benthic field data set.
#'@param transect_id String of transect names (e.g., ("A1", "A2", "A3", "B1", "B2", "B3")).
#'@param transect_length String of transect lengths in meters (e.g., c(10, 10, 10, 10, 10, 10)).
#'@param dbase_type Production database to use, either Indo-Pacific ReefBudget ("IPRB")
#'or U.S. Pacific Islands NCRMP-specific database ("NCRMP"). The Indo-Pacific ReefBudget
#'database is derived from "IP Calcification and bioerosion rates database v.1.3",
#'downloaded from https://geography.exeter.ac.uk/reefbudget/indopacific/.
#'@param method_name Transect design by which data were collected ("IPRB" or "Chords").
#'@param data_type Type of data collection ("In water" or "SfM").
#'
#'@import dplyr
#'@import tools
#'@import tidyr
#'@import reshape2
#'
#'@export process_prod


process_prod <- function(data,
                         transect_id = c("A1", "A2", "A3", "B1", "B2", "B3"),
                         transect_length = c(10, 10, 10, 10, 10, 10),
                         dbase_type = c("IPRB", "NCRMP"),
                         method_name = c("IPRB", "Chords"),
                         data_type = c("In water", "SfM"),
                         full_summary = FALSE,
                         label = NULL,
                         ...) {

  options(dplyr.summarise.inform = FALSE,
          scipen = 999)

  calc_prod_output <- run_calc_prod(data,
                                    transect_id,
                                    transect_length,
                                    dbase_type,
                                    data_type)
  data <- calc_prod_output$data
  transect_summary <- calc_prod_output$transect_summary
  prod_dbase <- calc_prod_output$prod_dbase

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
      "OCC_SITENAME",
      "LATITUDE",
      "LONGITUDE",
      "DEPTH_M",
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
