#' Calculate parrotfish erosion rates from fixed spc data
#'
#'@author Rebecca Weible
#'
#'@param data fish survey data.
#'@param rates_dbase Erosion rates database to use.
#'
#'@import Rmisc
#'@import sf
#'@import tidyverse
#'@import dplyr
#'
#'@export calc_fish_fixed_spc
#'
#'@examples
#'fish_data <- read.csv("CB_FishBelt_alldata.csv", na = "", check.names = FALSE)
#'
#'fish_fixed_spc <- calc_fish_fixed_spc(data = fish_data, rates_dbase = rates_dbase)

calc_fish_fixed_spc <- function(data,
                                rates_dbase = rates_dbase) {

  
    # FOR FIXED SPC DATA ----------------------------------------------------------------
    
    format_spc_output <- format_fish_spc(data,                                                 
                                         method = "fSPC",
                                         rates_dbase = rates_dbase)
    
    summary_spc_erosion  <- format_spc_output %>%
      pivot_longer(
        cols = c(
          SUM_BIOMASS_PER_FISH_KG_HECTARE,
          SUM_DENSITY_PER_FISH_HECTARE,
          SUM_EROSION_PER_FISH_KG_M2_YR
        ),
        names_to = "METRIC",
        values_to = "VALUE"
      ) %>%
      group_by(SITEVISITID, REA_SITEID, REP, COMMONFAMILYALL, FXN_GRP, METRIC) %>%
      summarize(
        MEAN = mean(VALUE, na.rm = TRUE),
        SD   = coalesce(sd(VALUE, na.rm = TRUE), 0),
        .groups = "drop"
      ) %>%
      mutate(
        SE  = SD / sqrt(2),
        L95 = pmax(0, MEAN - (SE * 1.97)),
        U95 = MEAN + (SE * 1.97)
      )

    format_fixed_spc_erosion <- summary_spc_erosion %>%
      filter(
        COMMONFAMILYALL != "NOTPARROTFISH",
        !is.na(REA_SITEID),
        !is.na(FXN_GRP)
      ) %>%
      pivot_longer(
        cols = -c(SITEVISITID, REA_SITEID, REP, COMMONFAMILYALL, FXN_GRP, METRIC),
        names_to = "calc",
        values_to = "value"
      ) %>%
      mutate(
        METRIC = case_match(METRIC,
                            "SUM_BIOMASS_PER_FISH_KG_HECTARE" ~ "FISH_BIOMASS_KG_HA",
                            "SUM_DENSITY_PER_FISH_HECTARE"    ~ "FISH_DENSITY_ABUNDANCE_HA",
                            "SUM_EROSION_PER_FISH_KG_M2_YR"   ~ "FISH_EROSION_KG_M2_YR",
                            .default = METRIC
        )
      ) %>%
      select(-SITEVISITID, -COMMONFAMILYALL, -REP) 
    
    options_fxn_grp <- c("ALL", "OTHER", "SCRAPER", "EXCAVATOR")
    
    format_fixed_spc_erosion2 <- format_fixed_spc_erosion %>%
      complete(
        nesting(REA_SITEID, METRIC, calc), # Keep the valid Site/Metric/Calc combinations observed in data
        FXN_GRP = options_fxn_grp,         # Expand this column to your specific list
        fill = list(value = 0)             # Auto-fill new rows with 0
      )
    

    # Create Master Metadata Table
    site_metadata <- data %>%
      filter(TRAINING_YN != "-1") %>%
      select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, 
             REA_SITEID, OCC_SITEID, LATITUDE, LONGITUDE) %>%
      distinct() %>%
      mutate(across(c(LATITUDE, LONGITUDE), ~round(., 5)))
    
    # 2. Reshape, Join, and Fill
    summary_fixed_spc_erosion <- format_fixed_spc_erosion2 %>%
      unite("metric_name", c(METRIC, FXN_GRP, calc), sep = "_") %>%
      # Pivot to Wide Format
      pivot_wider(
        id_cols = REA_SITEID, 
        names_from = metric_name, 
        values_from = value
      ) %>%
      right_join(site_metadata, by = "REA_SITEID") %>%
      mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
      filter(!is.na(OCC_SITEID)) %>%
      mutate(METHOD = as.factor("Fixed SPC")) %>%
      mutate(across(c(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE), as.factor)) %>%
      select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, 
             OCC_SITEID, LATITUDE, LONGITUDE, METHOD, everything(), -REA_SITEID)
  
    
    return(summary_fixed_spc_erosion)
  
}
