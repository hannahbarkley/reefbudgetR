#' Calculate parrotfish erosion rates from associated spc data
#'
#'@author Rebecca Weible
#'
#'@param data fish survey data.
#'@param rates_dbase Erosion rates database to use.
#'@param subset_distance_m Assigned associated site distances, in meters, from fixed SPC/OCC site to all other fish SPC sites, based on parrotfish foraging boundaries. 
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
                               subset_distance_m) {
  
  
    # FOR StRS SPC DATA ----------------------------------------------------------------
    
    format_strsspc_output <- format_fish_spc(data,                                                 
                                             method = "nSPC",
                                             rates_dbase = rates_dbase)
  
    # created associated SPC sites to each OCC fixed site SPC
    sites_associated_dbase_ <- create_fish_assoc_sites(data, subset_distance_m)
    sites_associated_dbase <- sites_associated_dbase_$output
    survey_sample_size <- sites_associated_dbase_$surveysamplesize # n of surveys that were assigned as associated or not (subset_distance_m from fixed site)
    
    summary_strsspc_erosion <- format_strsspc_output %>%
      # convert REPLICATEID values to Transect '1' and '2'
      dplyr::group_by(SITEVISITID) %>% 
      dplyr::mutate(TRANSECT = match(REPLICATEID, unique(REPLICATEID))) %>%
      tidyr::gather("METRIC", "VALUE", -c(SITEVISITID:FXN_GRP, TRANSECT)) %>% #set up data frame by spreading by transects (or the former REPLICATEID)
      dplyr::select(-REPLICATEID) %>%
      tidyr::spread(TRANSECT, VALUE) %>%
      dplyr::mutate_at(dplyr::vars(`1`, `2`), as.numeric) %>%
      # Average replicates (n=2)
      dplyr::group_by(SITEVISITID, REA_SITEID, REP, COMMONFAMILYALL, FXN_GRP, METRIC) %>% #calculate mean, sd, se, CI95 of data by averaging 2 divers (REPLICATEID)
      dplyr::rowwise() %>%
      dplyr::mutate(MEAN = mean(dplyr::c_across(`1`:`2`), na.rm = TRUE)) %>%
      dplyr::mutate(SD = sd(dplyr::c_across(`1`:`2`), na.rm = TRUE)) %>%
      dplyr::mutate(SD = dplyr::coalesce(SD, 0)) %>% # reaplce NA SD values with 0 because SD of one transect value is 0
      dplyr::mutate(SE = SD / sqrt(2)) %>% #2 is the number of total Transects
      dplyr::mutate(L95 = MEAN - (SE * 1.97)) %>%
      dplyr::mutate(U95 = MEAN + (SE * 1.97)) %>%
      dplyr::select(-c(`1`:`2`)) %>% #remove unnecessary columns
      dplyr::ungroup(.) %>%
      dplyr::mutate(L95 = dplyr::case_when(L95 < 0 ~ 0,
                                    TRUE ~ as.numeric(L95)))
    
    format_strs_spc_erosion <- summary_strsspc_erosion %>%
      dplyr::select(-c(SD:U95)) %>%
      #add ASSOC_OCCSITE before averaging by replicate
      merge(., sites_associated_dbase %>% 
                         dplyr::filter(CB_METHOD %in% "nSPC") %>%
                         dplyr::select(REA_SITEID, ASSOC_OCCSITEID, value), by = "REA_SITEID") %>% # attach associated occ site ID 
      dplyr::mutate(ASSOC_OCCSITEID = dplyr::case_when(value == 0 ~ "",
                                       TRUE ~ ASSOC_OCCSITEID)) %>% # create assoc_occsite ID names...NOTE: the NA in assoc_occsite represents actual fixed site
      dplyr::filter(!value %in% c("0", "1")) %>%
      dplyr::filter(!is.na(value)) %>%
      dplyr::select(-value) %>%
      # convert SITE values to Transect '1' through `n`
      dplyr::group_by(ASSOC_OCCSITEID) %>% 
      dplyr::mutate(TRANSECT = match(REA_SITEID, unique(REA_SITEID))) %>%
      dplyr::select(-SITEVISITID, -REA_SITEID, -REP)
    
    calc_strs_spc_erosion <- format_strs_spc_erosion %>%
      # create column for total transects per associated site ID because they vary
      dplyr::left_join(., format_strs_spc_erosion %>% 
                  dplyr::group_by(ASSOC_OCCSITEID) %>%
                  dplyr::summarise(n = dplyr::n_distinct(TRANSECT)), 
                by = "ASSOC_OCCSITEID") %>%      
      tidyr::spread(TRANSECT, MEAN) %>%
      dplyr::mutate_at(dplyr::vars(n:length(.)), as.numeric) %>%
      dplyr::ungroup(.) %>%
      # Average transects (n = variable) per associated site ID
      dplyr::group_by(ASSOC_OCCSITEID, COMMONFAMILYALL, FXN_GRP, METRIC) %>% #calculate mean, sd, se, CI95 of data by averaging 2 divers (REPLICATEID)
      dplyr::rowwise() %>%
      dplyr::mutate(MEAN = mean(dplyr::c_across(`1`:(length(.)-5)), na.rm = TRUE)) %>%
      dplyr::mutate(SD = sd(dplyr::c_across(`1`:(length(.)-6)), na.rm = TRUE)) %>%
      dplyr::mutate(SD = dplyr::coalesce(SD, 0)) %>% # reaplce NA SD values with 0 because SD of one transect value is 0
      dplyr::mutate(SE = SD / sqrt(n)) %>% #2 is the number of total Transects
      dplyr::mutate(L95 = MEAN - (SE * 1.97)) %>%
      dplyr::mutate(U95 = MEAN + (SE * 1.97)) %>%
      dplyr::ungroup(.) %>%
      dplyr::select(c(COMMONFAMILYALL:ASSOC_OCCSITEID, MEAN:U95)) %>% #remove unnecessary columns
      dplyr::mutate(L95 = dplyr::case_when(L95 < 0 ~ 0,
                                    TRUE ~ as.numeric(L95))) %>%
      dplyr::filter(!COMMONFAMILYALL %in% "NOTPARROTFISH") %>% # removed non-parrotfish species group
      dplyr::filter(!ASSOC_OCCSITEID %in% NA) %>% # remove NA sites...these correspond to the fixed site data as well
      dplyr::filter(!FXN_GRP %in% NA) %>% # remove NA Graz_types (these are also all non-parrotfish species group)
      # Format dataframe to match BELT FINAL BIOEROSION OUTPUT
      tidyr::gather(., "calc", "value", -c(COMMONFAMILYALL:ASSOC_OCCSITEID)) %>%
      dplyr::select(-COMMONFAMILYALL) %>% #remove n and REP
      dplyr::mutate(METRIC = dplyr::case_when(METRIC == "SUM_BIOMASS_PER_FISH_KG_HECTARE" ~ "FISH_BIOMASS_KG_HA",
                                              METRIC == "SUM_DENSITY_PER_FISH_HECTARE" ~ "FISH_DENSITY_ABUNDANCE_HA",
                                              METRIC == "SUM_EROSION_PER_FISH_KG_M2_YR" ~ "FISH_EROSION_KG_M2_YR")) 
    
    # pause here to complete FXN_GRP column if values are missing
    options_fxn_grp <- c("ALL", "OTHER", "SCRAPER", "EXCAVATOR") # all four functional groups should be present
    missing_fxn_grp <- options_fxn_grp[!options_fxn_grp %in% unique(calc_strs_spc_erosion$FXN_GRP)] #isolate missing functional groups
    
    # Add back in rows of missing functional groups and assign value to 0 with completing the dataframe
    if(length(missing_fxn_grp) >= 1) {
      for(i in 1:length(missing_fxn_grp)) {
        
        calc_strs_spc_erosion <- add_row(calc_strs_spc_erosion, FXN_GRP = missing_fxn_grp[i])
        
      }
      
      calc_strs_spc_erosion2 <- calc_strs_spc_erosion %>%
        select(-SITEVISITID) %>%
        complete(REA_SITEID, METRIC, calc, FXN_GRP,
                 fill = list(value = 0)
        ) %>%
        filter_all(all_vars(!is.na(.))) %>%
        distinct(.)
      
    } else {calc_strs_spc_erosion2 <- calc_strs_spc_erosion}
    
    
    # Add back in occ fixed sites that are missing and assign value to 0 with completing the dataframe
    iprb_data <- data %>% dplyr::filter(CB_METHOD %in% "IPRB")
    missing_sites <- unique(iprb_data$OCC_SITEID)[!unique(iprb_data$OCC_SITEID) %in% unique(calc_strs_spc_erosion2$ASSOC_OCCSITEID)]
    
    if(length(missing_sites) >= 1){
      for(i in 1:length(missing_sites)){
        
        calc_strs_spc_erosion2 <- dplyr::add_row(calc_strs_spc_erosion2, ASSOC_OCCSITEID = missing_sites[i])
        
      }
      
      calc_strs_spc_erosion3 <- calc_strs_spc_erosion2 %>%
        tidyr::complete(FXN_GRP, METRIC, ASSOC_OCCSITEID, calc,
                        fill = list(value = 0)) %>%
        dplyr::filter_all(dplyr::all_vars(!is.na(.))) %>%
        dplyr::distinct(.)
      
      
    } else {calc_strs_spc_erosion3 <- calc_strs_spc_erosion2}
    
    
    
    # Continue formatting dataframe to match BELT FINAL BIOEROSION OUTPUT
    calc_strs_spc_erosion4 <-   
      calc_strs_spc_erosion3 %>%
      tidyr::unite("METRIC", METRIC, FXN_GRP) %>% 
      tidyr::unite("METRIC", METRIC, calc) %>%
      # fill in missing metadata info
      dplyr::left_join(., data %>% dplyr::select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID, LATITUDE, LONGITUDE, CB_METHOD), 
                       by = c("ASSOC_OCCSITEID" = "OCC_SITEID")) %>% #join important metadata that was lost during averaging replicates
      dplyr::distinct(.) %>% #left_join creates duplicates, so remove duplicates
      dplyr::mutate(LATITUDE = round(LATITUDE, 5)) %>%
      dplyr::mutate(LONGITUDE = round(LONGITUDE, 5)) %>%
      # rename(REGION = REGION_NAME,
      #        CRUISE_ID = MISSIONID,
      #        LOCATION = ISLAND,
      #        CB_METHOD = METHOD) %>%
      # mutate(REGIONCODE = loc,
      #        LOCATIONCODE = str_sub(ASSOC_OCCSITE, 5,7)) %>%
      dplyr::rename(OCC_SITEID = ASSOC_OCCSITEID) %>%
      dplyr::left_join(., data %>% dplyr::select(OCC_SITEID, OCC_SITENAME), by = "OCC_SITEID") %>%
      dplyr::distinct(.) %>%
      dplyr::mutate(CB_METHOD = "StRS SPC") %>%
      tidyr::spread(., METRIC, value, fill = "0") %>%
      dplyr::select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID, OCC_SITENAME, LATITUDE, LONGITUDE, CB_METHOD, everything(.)) %>%
      dplyr::mutate_at(dplyr::vars(REGION:CB_METHOD), as.factor) %>%
      dplyr::mutate_at(dplyr::vars(FISH_BIOMASS_KG_HA_ALL_L95:FISH_EROSION_KG_M2_YR_SCRAPER_U95), as.numeric)
    
    return(list(calc_strs_ero <- calc_strs_spc_erosion4, assoc_survey_count <- survey_sample_size))
    
  
}


