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
                               rates_dbase = rates_dbase,
                               fixed_metadata = NULL,
                               subset_distance_m = NULL,
                               method_type = method_type) {
  # FOR StRS SPC DATA ----------------------------------------------------------------
  
  if (method_type == "StRS") {
    strs_metadata <- data %>%
      dplyr::select(
        CRUISE_ID,
        REGION,
        REGIONCODE,
        LOCATION,
        LOCATIONCODE,
        LATITUDE,
        LONGITUDE,
        YEAR,
        METHOD,
        REA_SITEID
      ) %>%
      distinct()
    
    format_strsspc_output <- format_fish_spc(data, method = "nSPC", rates_dbase = rates_dbase)
    
    summary_strsspc_erosion <- format_strsspc_output %>%
      tidyr::pivot_longer(
        cols = c(
          SUM_BIOMASS_PER_FISH_KG_HECTARE,
          SUM_DENSITY_PER_FISH_HECTARE,
          SUM_EROSION_PER_FISH_KG_M2_YR
        ),
        names_to = "METRIC",
        values_to = "VALUE"
      ) %>%
      dplyr::group_by(SITEVISITID,
                      REA_SITEID,
                      REP,
                      COMMONFAMILYALL,
                      FXN_GRP,
                      METRIC) %>%
      dplyr::summarize(
        MEAN = mean(VALUE, na.rm = TRUE),
        SD   = dplyr::coalesce(sd(VALUE, na.rm = TRUE), 0),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        SE  = SD / sqrt(2),
        L95 = pmax(0, MEAN - (SE * 1.97)),
        U95 = MEAN + (SE * 1.97)
      ) %>%
      dplyr::select(-SITEVISITID, -REP)
    
    calc_strs_spc_erosion <- summary_strsspc_erosion %>%
      dplyr::mutate(L95 = pmax(0, L95)) %>%
      dplyr::filter(COMMONFAMILYALL != "NOTPARROTFISH", !is.na(FXN_GRP)) %>%
      tidyr::pivot_longer(
        cols = -c(REA_SITEID:METRIC),
        # Keeps your negative selection logic
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
    
    options_fxn_grp <- c("ALL", "OTHER", "SCRAPER", "EXCAVATOR")
    
    target_sites <- strs_metadata %>%
      filter(LOCATIONCODE %in% data$LOCATIONCODE) %>%
      pull(REA_SITEID) %>%
      unique()
    
    calc_strs_spc_erosion3 <- calc_strs_spc_erosion %>%
      dplyr::select(-any_of("SITEVISITID")) %>%
      tidyr::complete(
        REA_SITEID = target_sites,
        FXN_GRP = options_fxn_grp,
        METRIC,
        calc,
        fill = list(value = 0)
      ) %>%
      dplyr::filter(!is.na(METRIC), !is.na(calc)) %>%
      dplyr::distinct(.)
    
    calc_strs_spc_erosion4 <- calc_strs_spc_erosion3 %>%
      dplyr::left_join(
        strs_metadata %>%
          dplyr::select(
            REGION,
            REGIONCODE,
            CRUISE_ID,
            LOCATION,
            LOCATIONCODE,
            REA_SITEID,
            LATITUDE,
            LONGITUDE
          ),
        by = "REA_SITEID",
        relationship = "many-to-one"
      ) %>%
      tidyr::pivot_wider(
        id_cols = c(
          REGION,
          REGIONCODE,
          CRUISE_ID,
          LOCATION,
          LOCATIONCODE,
          REA_SITEID,
          LATITUDE,
          LONGITUDE
        ),
        names_from = c(METRIC, FXN_GRP, calc),
        values_from = value,
        names_sep = "_",
        values_fill = 0
      ) %>%
      dplyr::mutate(
        LATITUDE = round(LATITUDE, 5),
        LONGITUDE = round(LONGITUDE, 5),
        METHOD = as.factor("StRS SPC")
      ) %>%
      dplyr::mutate(across(
        c(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE),
        as.factor
      )) %>%
      
      dplyr::select(
        REGION,
        REGIONCODE,
        CRUISE_ID,
        LOCATION,
        LOCATIONCODE,
        REA_SITEID,
        LATITUDE,
        LONGITUDE,
        METHOD,
        everything()
      )
    
    return(list(calc_strs_ero = calc_strs_spc_erosion4))
    
  }
  
  # FOR MEAN StRS SPC DATA ----------------------------------------------------------------
  
  if (method_type == "Mean StRS") {
    format_strsspc_output <- format_fish_spc(data, method = "nSPC", rates_dbase = rates_dbase)
    
    # created associated SPC sites to each OCC fixed site SPC
    sites_associated_dbase_ <- suppressWarnings(create_fish_assoc_sites(data, fixed_metadata, subset_distance_m))
    sites_associated_dbase <- sites_associated_dbase_$output
    survey_sample_size <- sites_associated_dbase_$surveysamplesize
    
    summary_strsspc_erosion <- format_strsspc_output %>%
      tidyr::pivot_longer(
        cols = -c(SITEVISITID:FXN_GRP, REPLICATEID),
        names_to = "METRIC",
        values_to = "VALUE"
      ) %>%
      dplyr::mutate(VALUE = as.numeric(VALUE)) %>%
      dplyr::group_by(SITEVISITID,
                      REA_SITEID,
                      REP,
                      COMMONFAMILYALL,
                      FXN_GRP,
                      METRIC) %>%
      dplyr::summarise(MEAN = mean(VALUE, na.rm = TRUE),
                       SD = sd(VALUE, na.rm = TRUE),) %>%
      dplyr::mutate(
        SD = dplyr::coalesce(SD, 0),
        SE = SD / sqrt(2),
        L95 = pmax(0, MEAN - (SE * 1.97)),
        U95 = MEAN + (SE * 1.97)
      )
    
    format_strs_spc_erosion <- summary_strsspc_erosion %>%
      dplyr::select(-c(SD:U95)) %>%
      dplyr::inner_join(
        sites_associated_dbase %>%
          dplyr::filter(METHOD == "nSPC") %>%
          dplyr::filter(!value %in% c("0", "1", 0, 1) &
                          !is.na(value)) %>%
          dplyr::select(REA_SITEID, ASSOC_OCCSITEID),
        by = "REA_SITEID",
        relationship = "many-to-many"
      ) %>%
      dplyr::group_by(ASSOC_OCCSITEID) %>%
      dplyr::mutate(TRANSECT = dplyr::dense_rank(REA_SITEID)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-SITEVISITID, -REA_SITEID, -REP)
    
    calc_strs_spc_erosion <- format_strs_spc_erosion %>%
      dplyr::group_by(ASSOC_OCCSITEID, COMMONFAMILYALL, FXN_GRP, METRIC) %>%
      dplyr::summarise(
        n = n(),
        MEAN_val = mean(MEAN, na.rm = TRUE),
        SD = sd(MEAN, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        SD = dplyr::coalesce(SD, 0),
        # Replace NA SD with 0
        SE = SD / sqrt(n),
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
          METRIC == "SUM_DENSITY_PER_FISH_HECTARE" ~ "FISH_DENSITY_ABUNDANCE_HA",
          METRIC == "SUM_EROSION_PER_FISH_KG_M2_YR" ~ "FISH_EROSION_KG_M2_YR",
          TRUE ~ METRIC
        )
      )
    
    calc_strs_spc_erosion2 <- calc_strs_spc_erosion %>%
      dplyr::ungroup() %>%
      dplyr::mutate(FXN_GRP = factor(FXN_GRP, levels = c(
        "ALL", "OTHER", "SCRAPER", "EXCAVATOR"
      ))) %>%
      tidyr::complete(ASSOC_OCCSITEID, 
                      METRIC,
                      calc, 
                      FXN_GRP, 
                      fill = list(value = 0)) %>%
      dplyr::mutate(FXN_GRP = as.character(FXN_GRP)) %>%
      dplyr::filter(!is.na(ASSOC_OCCSITEID)) %>%
      dplyr::distinct()
    
    
    all_fixed_sites <- fixed_metadata %>%
      dplyr::filter(LOCATIONCODE %in% unique(data$LOCATIONCODE)) %>%
      dplyr::pull(OCC_SITEID) %>%
      unique()
    
    calc_strs_spc_erosion3 <- calc_strs_spc_erosion2 %>%
      tidyr::complete(ASSOC_OCCSITEID = all_fixed_sites,
                      FXN_GRP,
                      METRIC,
                      calc,
                      fill = list(value = 0)) %>%
      dplyr::filter(!is.na(FXN_GRP) & !is.na(METRIC)) %>%
      dplyr::distinct()
    
    
    # Continue formatting dataframe to match BELT FINAL BIOEROSION OUTPUT
    calc_strs_spc_erosion4 <- calc_strs_spc_erosion3 %>%
      dplyr::left_join(
        fixed_metadata %>%
          dplyr::select(
            REGION,
            REGIONCODE,
            CRUISE_ID,
            LOCATION,
            LOCATIONCODE,
            OCC_SITEID,
            LATITUDE,
            LONGITUDE
          ) %>%
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
        values_fill = 0 #
      ) %>%
      dplyr::select(
        REGION,
        REGIONCODE,
        CRUISE_ID,
        LOCATION,
        LOCATIONCODE,
        OCC_SITEID,
        LATITUDE,
        LONGITUDE,
        METHOD,
        dplyr::everything()
      ) %>%
      dplyr::mutate(dplyr::across(
        c(
          REGION,
          REGIONCODE,
          CRUISE_ID,
          LOCATION,
          LOCATIONCODE,
          METHOD
        ),
        as.factor
      ),
      dplyr::across(where(is.numeric), as.numeric))
    
    return(
      list(calc_strs_ero = calc_strs_spc_erosion4, assoc_survey_count = survey_sample_size)
    )
    
    
  }
}
