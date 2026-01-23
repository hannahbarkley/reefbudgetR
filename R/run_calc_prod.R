#' Run calc_prod function
#'
#'@author Hannah Barkley
#'
#'@param data Benthic data set to process
#'@param dbase_type Production database to use, either Indo-Pacific ReefBudget ("IPRB"), the
#'U.S. Pacific Islands NCRMP-specific database ("NCRMP"), or upload a customized database with location-specific rates ("Custom").
#'The Indo-Pacific ReefBudget database is derived from "IP Calcification and bioerosion rates database v.1.3",
#'downloaded from https://geography.exeter.ac.uk/reefbudget/indopacific/.
#'
#'@import dplyr
#'@importFrom sjmisc seq_row
#'
#'@export run_calc_prod
#'
#' @examples
#' benthic_data <- read.csv("ESD_CarbBudget_Benthic_OAHU_2021.csv",
#'     na = "", check.names = FALSE)
#'
#' calc_prod_output <- run_calc_prod(
#'     data = benthic_data,
#'     dbase_type = "NCRMP"
#'     )


run_calc_prod <- function(data, dbase_type, prod_dbase_custom = NULL, ...) {
  if (dbase_type == "IPRB") {
    data$SUBSTRATE_CODE <- data$SUBSTRATE_CODE_IPRB
    prod_dbase <- prod_dbase_iprb
  }
  
  if (dbase_type == "NCRMP") {
    prod_dbase <- prod_dbase_ncrmp
  }
  
  if (dbase_type == "Custom") {
    prod_dbase <- prod_dbase_custom
  }
  
  # Create a transect summary
  transect_summary <- data %>%
    group_by(
      CRUISE_ID,
      REGION,
      REGIONCODE,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      LOCALDATE,
      YEAR,
      SITEVISITID,
      LATITUDE,
      LONGITUDE,
      SITE_DEPTH_M,
      CB_METHOD,
      CB_TRANSECTID
    ) %>%
    summarize(
      TRANSECT_PLANAR_LENGTH_M = mean(as.numeric(TRANSECT_PLANAR_LENGTH_M), na.rm = TRUE),
      TRANSECT_TOTAL_SUBSTRATE_COVER_M = sum(SUBSTRATE_COVER_CM / 100, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      TRANSECT_RUGOSITY = TRANSECT_TOTAL_SUBSTRATE_COVER_M / TRANSECT_PLANAR_LENGTH_M,
      OCC_SITEID_TRANSECT = str_c(OCC_SITEID, CB_TRANSECTID, sep = "-")
    )
  
  # Add any missing metadata to benthic obs using prod_dbase
  data <- data %>%
    left_join(
      prod_dbase %>%
        select(
          SUBSTRATE_CODE,
          SUBSTRATE_NAME,
          SUBSTRATE_CLASS,
          TAXA_LEVEL,
          CORAL_GROUP,
          CORAL_GROUP_NAME
        ) %>%
        distinct(SUBSTRATE_CODE, .keep_all = TRUE),
      by = "SUBSTRATE_CODE",
      suffix = c("", ".db")
    ) %>%
    mutate(
      SUBSTRATE_NAME   = if ("SUBSTRATE_NAME.db" %in% names(.))
        coalesce(SUBSTRATE_NAME, SUBSTRATE_NAME.db)
      else
        SUBSTRATE_NAME,
      SUBSTRATE_CLASS  = if ("SUBSTRATE_CLASS.db" %in% names(.))
        coalesce(SUBSTRATE_CLASS, SUBSTRATE_CLASS.db)
      else
        SUBSTRATE_CLASS,
      TAXA_LEVEL       = if ("TAXA_LEVEL.db" %in% names(.))
        coalesce(TAXA_LEVEL, TAXA_LEVEL.db)
      else
        TAXA_LEVEL,
      CORAL_GROUP      = if ("CORAL_GROUP.db" %in% names(.))
        coalesce(CORAL_GROUP, CORAL_GROUP.db)
      else
        CORAL_GROUP,
      CORAL_GROUP_NAME = if ("CORAL_GROUP_NAME.db" %in% names(.))
        coalesce(CORAL_GROUP_NAME, CORAL_GROUP_NAME.db)
      else
        CORAL_GROUP_NAME
    ) %>%
    select(-ends_with(".db")) %>%
    mutate(
      MORPHOLOGYCODE = case_when(
        SUBSTRATE_CODE %in% c("PRUS", "PMRC") ~ "LC",!is.na(MORPHOLOGYCODE) ~ MORPHOLOGYCODE,
        TAXA_LEVEL == "SPECIES" ~ prod_dbase$MORPHOLOGYCODE[match(SUBSTRATE_CODE, prod_dbase$SUBSTRATE_CODE)],
        TRUE ~ NA_character_
      )
    )
  
  # Throw an error if GENUS is missing Morphology. This stops the script and names the offending substrate codes.
  {
    genus_errors <- data %>%
      filter(is.na(MORPHOLOGYCODE),
             TAXA_LEVEL == "GENUS",
             SUBSTRATE_CLASS == "CORAL")
    
    if (nrow(genus_errors) > 0) {
      stop(paste(
        "Error: GENUS entries missing MORPHOLOGYCODE for codes:",
        paste(unique(genus_errors$SUBSTRATE_CODE), collapse = ", ")
      ))
    }
    }
  
  # Use the resolved MORPHOLOGYCODE to get CORAL_GROUP and CORAL_GROUP_NAME
  data <- data %>%
    left_join(
      prod_dbase %>%
        select(
          SUBSTRATE_CODE,
          MORPHOLOGYCODE,
          CORAL_GROUP,
          CORAL_GROUP_NAME,
          MORPHOLOGY
        ) %>%
        distinct(),
      by = c("SUBSTRATE_CODE", "MORPHOLOGYCODE"),
      suffix = c("", ".db")
    ) %>%
    mutate(
      CORAL_GROUP = if_else(
        SUBSTRATE_CODE %in% c("PRUS", "PMRC"),
        "POLC",
        if ("CORAL_GROUP.db" %in% names(.))
          coalesce(CORAL_GROUP, CORAL_GROUP.db)
        else
          CORAL_GROUP
      ),
      CORAL_GROUP_NAME = if_else(
        SUBSTRATE_CODE %in% c("PRUS", "PMRC"),
        "Laminar columnar Porites",
        if ("CORAL_GROUP_NAME.db" %in% names(.))
          coalesce(CORAL_GROUP_NAME, CORAL_GROUP_NAME.db)
        else
          CORAL_GROUP_NAME
      ),
      MORPHOLOGY = if_else(
        SUBSTRATE_CODE %in% c("PRUS", "PMRC"),
        "Laminar columnar",
        if ("MORPHOLOGY.db" %in% names(.))
          coalesce(MORPHOLOGY, MORPHOLOGY.db)
        else
          MORPHOLOGY
      )
    ) %>%
    select(-ends_with(".db"))
  
  data$COLONY_PROD_G_YR     <- NA_real_
  data$COLONY_PROD_G_YR_L95 <- NA_real_
  data$COLONY_PROD_G_YR_U95 <- NA_real_
  
  # Create an empty list to store details about rows that failed
  error_log <- data.frame(
    row_index = integer(),
    substrate_code = character(),
    error_message = character()
  )
  
  for (i in seq_len(nrow(data))) {
    
    loc_code <- if (dbase_type == "Custom") data$LOCATIONCODE[i] else NULL
    
    calc_i <- tryCatch({
      calc_prod(
        substrate_class    = data$SUBSTRATE_CLASS[i],
        substrate_code     = data$SUBSTRATE_CODE[i],
        morphology_code    = data$MORPHOLOGYCODE[i],
        substrate_cover_cm = data$SUBSTRATE_COVER_CM[i],
        location_code      = loc_code,
        dbase_type         = dbase_type,
        prod_dbase         = prod_dbase
      )
    }, error = function(e) {
      # If an error occurs, log it and return NULL
      error_log <<- rbind(
        error_log,
        data.frame(
          row_index = i,
          substrate_code = data$SUBSTRATE_CODE[i],
          error_message = as.character(e$message)
        )
      )
      return(NULL)
    })
    
    # Only assign results if calc_prod succeeded (didn't return NULL)
    if (!is.null(calc_i)) {
      data$COLONY_PROD_G_YR[i]     <- calc_i$cp_i
      data$COLONY_PROD_G_YR_L95[i] <- calc_i$cp_l95_i
      data$COLONY_PROD_G_YR_U95[i] <- calc_i$cp_u95_i
    }
  }
  
  data <- data %>%
    mutate(OCC_SITEID_TRANSECT = paste(OCC_SITEID, CB_TRANSECTID, sep = "-")) %>%
    left_join(
      transect_summary %>%
        select(OCC_SITEID_TRANSECT, TRANSECT_TOTAL_SUBSTRATE_COVER_M),
      by = "OCC_SITEID_TRANSECT"
    )
  
  # Return the results and the error log
  return(list(
    data = data,
    transect_summary = transect_summary,
    error_log = error_log
  ))
}
