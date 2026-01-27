#' Calculate parrotfish erosion rates from associated spc data
#'
#'@author Rebecca Weible
#'
#'@param data fish survey data.
#'@param rates_dbase Erosion rates database to use.
#'@param fixed_metadata Dataframe with fixed site info
#'@param subset_distance_m Assigned associated site distances, in meters, from fixed SPC/OCC site to all other fish SPC sites, based on parrotfish foraging boundaries.
#'@param method_type Calculate erosion rates for each StRS survey ("StRS") or average all surveys within a specified distance ("Mean StRS")
#'
#'@import Rmisc
#'@import sf
#'@import tidyverse
#'@import dplyr
#'
#'@export calc_fish_strs_spc
#'
#'@examples
#'fish_data <- read.csv("CB_FishBelt_alldata.csv", na = "", check.names = FALSE)
#'
#'fish_strs_spc <- calc_fish_strs_spc(data = fish_data, rates_dbase = rates_dbase, subset_distance_m = 2000)

calc_fish_strs_spc <- function(data,
                               rates_dbase,
                               fixed_metadata = NULL,
                               subset_distance_m = NULL,
                               method_type) {
  
  # ============================================================================
  # StRS SPC DATA (Stratified Random Sampling)
  # ============================================================================
  
  if (method_type == "StRS") {
    
    # Format and Summary Stats
    summary_strsspc <- format_fish_spc(data, method = "nSPC", rates_dbase = rates_dbase) %>%
      tidyr::pivot_longer(
        cols = c(SUM_BIOMASS_PER_FISH_KG_HECTARE, SUM_DENSITY_PER_FISH_HECTARE, SUM_EROSION_PER_FISH_KG_M2_YR),
        names_to = "METRIC",
        values_to = "VALUE"
      ) %>%
      dplyr::group_by(SITEVISITID, REA_SITEID, REP, COMMONFAMILYALL, FXN_GRP, METRIC) %>%
      dplyr::summarise(
        MEAN = mean(VALUE, na.rm = TRUE),
        SD   = sd(VALUE, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        SD  = dplyr::coalesce(SD, 0),
        SE  = SD / sqrt(2),
        L95 = pmax(0, MEAN - (SE * 1.97)),
        U95 = MEAN + (SE * 1.97)
      )
    
    # Filter, Rename, Pivot
    calc_strs_spc <- summary_strsspc %>%
      dplyr::filter(COMMONFAMILYALL != "NOTPARROTFISH", !is.na(FXN_GRP)) %>%
      tidyr::pivot_longer(
        cols = c(MEAN, SD, SE, L95, U95),
        names_to = "calc",
        values_to = "value"
      ) %>%
      dplyr::select(-COMMONFAMILYALL) %>%
      dplyr::mutate(
        METRIC = dplyr::case_match(
          METRIC,
          "SUM_BIOMASS_PER_FISH_KG_HECTARE" ~ "FISH_BIOMASS_KG_HA",
          "SUM_DENSITY_PER_FISH_HECTARE"    ~ "FISH_DENSITY_ABUNDANCE_HA",
          "SUM_EROSION_PER_FISH_KG_M2_YR"   ~ "FISH_EROSION_KG_M2_YR",
          .default = METRIC
        )
      )
    
    # Complete Missing Data
    target_sites <- unique(data$REA_SITEID) 
    options_fxn_grp <- c("ALL", "OTHER", "SCRAPER", "EXCAVATOR")
    
    calc_strs_spc_filled <- calc_strs_spc %>%
      dplyr::select(-any_of(c("SITEVISITID", "REP"))) %>%
      tidyr::complete(
        REA_SITEID = target_sites,
        FXN_GRP = options_fxn_grp,
        METRIC,
        calc,
        fill = list(value = 0)
      )
    
    # Join Metadata and Format Final
    meta_cols <- data %>%
      dplyr::select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, REA_SITEID, LATITUDE, LONGITUDE) %>%
      dplyr::distinct()
    
    calc_final <- calc_strs_spc_filled %>%
      dplyr::left_join(meta_cols, by = "REA_SITEID") %>%
      dplyr::mutate(
        LATITUDE = round(LATITUDE, 5),
        LONGITUDE = round(LONGITUDE, 5),
        METHOD = "StRS SPC"
      ) %>%
      tidyr::pivot_wider(
        names_from = c(METRIC, FXN_GRP, calc),
        names_sep = "_",
        values_from = value,
        values_fill = 0
      ) %>%
      dplyr::mutate(
        dplyr::across(c(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, METHOD), as.factor),
        dplyr::across(where(is.numeric), as.numeric)
      ) %>%
      dplyr::select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, REA_SITEID, LATITUDE, LONGITUDE, METHOD, everything())
    
    return(list(calc_strs_ero = calc_final))
  }
  
  
  # ============================================================================
  # MEAN StRS SPC DATA (Associated Sites)
  # ============================================================================
  
  if (method_type == "Mean StRS") {
    
    # Base Calculations
    format_strsspc_output <- format_fish_spc(data, method = "nSPC", rates_dbase = rates_dbase)
    
    # Create Associations
    sites_assoc <- suppressWarnings(create_fish_assoc_sites(data, fixed_metadata, subset_distance_m))
    sites_assoc_db <- sites_assoc$output
    survey_sample_size <- sites_assoc$surveysamplesize
    
    # Initial Summary (Site Level)
    summary_site_level <- format_strsspc_output %>%
      tidyr::pivot_longer(
        cols = -c(SITEVISITID:FXN_GRP, REPLICATEID),
        names_to = "METRIC",
        values_to = "VALUE"
      ) %>%
      dplyr::mutate(VALUE = as.numeric(VALUE)) %>%
      dplyr::group_by(SITEVISITID, REA_SITEID, REP, COMMONFAMILYALL, FXN_GRP, METRIC) %>%
      dplyr::summarise(
        MEAN = mean(VALUE, na.rm = TRUE),
        SD   = sd(VALUE, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        SD  = dplyr::coalesce(SD, 0),
        SE  = SD / sqrt(2),
        L95 = pmax(0, MEAN - (SE * 1.97)),
        U95 = MEAN + (SE * 1.97)
      )
    
    # Associate and Filter
    format_assoc <- summary_site_level %>%
      dplyr::select(-SD, -SE, -L95, -U95) %>%
      dplyr::inner_join(
        sites_assoc_db %>%
          dplyr::filter(METHOD == "nSPC") %>%
          dplyr::filter(!value %in% c("0", "1", 0, 1) & !is.na(value)) %>%
          dplyr::select(REA_SITEID, ASSOC_OCCSITEID),
        by = "REA_SITEID",
        relationship = "many-to-many"
      ) %>%
      dplyr::select(-SITEVISITID, -REA_SITEID, -REP)
    
    # Secondary Aggregation (Associated Site Level)
    calc_assoc_level <- format_assoc %>%
      dplyr::group_by(ASSOC_OCCSITEID, COMMONFAMILYALL, FXN_GRP, METRIC) %>%
      dplyr::summarise(
        n = n(),
        MEAN_val = mean(MEAN, na.rm = TRUE),
        SD = sd(MEAN, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        SD  = dplyr::coalesce(SD, 0),
        SE  = SD / sqrt(n),
        L95 = pmax(0, MEAN_val - (SE * 1.97)),
        U95 = MEAN_val + (SE * 1.97)
      ) %>%
      dplyr::rename(MEAN = MEAN_val) %>%
      dplyr::filter(COMMONFAMILYALL != "NOTPARROTFISH") %>%
      dplyr::filter(!is.na(ASSOC_OCCSITEID) & !is.na(FXN_GRP)) %>%
      tidyr::pivot_longer(
        cols = c(MEAN, SD, SE, L95, U95),
        names_to = "calc",
        values_to = "value"
      ) %>%
      dplyr::select(-COMMONFAMILYALL, -n) %>%
      dplyr::mutate(
        METRIC = dplyr::case_when(
          METRIC == "SUM_BIOMASS_PER_FISH_KG_HECTARE" ~ "FISH_BIOMASS_KG_HA",
          METRIC == "SUM_DENSITY_PER_FISH_HECTARE"    ~ "FISH_DENSITY_ABUNDANCE_HA",
          METRIC == "SUM_EROSION_PER_FISH_KG_M2_YR"     ~ "FISH_EROSION_KG_M2_YR",
          TRUE ~ METRIC
        )
      )
    
    # Complete Missing Data
    all_fixed_sites <- fixed_metadata %>%
      dplyr::filter(LOCATIONCODE %in% unique(data$LOCATIONCODE)) %>%
      dplyr::pull(OCC_SITEID) %>%
      unique()
    
    options_fxn_grp <- c("ALL", "OTHER", "SCRAPER", "EXCAVATOR")
    
    calc_filled <- calc_assoc_level %>%
      dplyr::ungroup() %>%
      tidyr::complete(
        ASSOC_OCCSITEID = all_fixed_sites, # Adds missing sites
        FXN_GRP = options_fxn_grp,         # Adds missing groups
        METRIC,
        calc,
        fill = list(value = 0)
      ) %>%
      dplyr::filter(!is.na(METRIC)) %>%
      dplyr::distinct()
    
    # Final Formatting and Metadata Join
    calc_final_mean <- calc_filled %>%
      dplyr::left_join(
        fixed_metadata %>%
          dplyr::select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID, LATITUDE, LONGITUDE) %>%
          dplyr::distinct(),
        by = c("ASSOC_OCCSITEID" = "OCC_SITEID")
      ) %>%
      dplyr::mutate(
        LATITUDE = round(LATITUDE, 5),
        LONGITUDE = round(LONGITUDE, 5),
        METHOD = "StRS SPC",
        OCC_SITEID = ASSOC_OCCSITEID
      ) %>%
      tidyr::pivot_wider(
        names_from = c(METRIC, FXN_GRP, calc),
        names_sep = "_",
        values_from = value,
        values_fill = 0
      ) %>%
      dplyr::mutate(
        dplyr::across(c(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, METHOD), as.factor),
        dplyr::across(where(is.numeric), as.numeric)
      ) %>%
      dplyr::select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID, LATITUDE, LONGITUDE, METHOD, everything())
    
    return(list(calc_strs_ero = calc_final_mean, assoc_survey_count = survey_sample_size))
  }
}