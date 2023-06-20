#' Calculate parrotfish erosion rates from associated spc data
#'
#'@author Rebecca Weible
#'
#'@param data fish survey data.
#'@param dbase_type Erosion rates database to use. Choose either Indo-Pacific
#'ReefBudget ("dbase_type = "IPRB") or U.S. Pacific Islands rates developed
#'by Tye Kindinger, NOAA PIFSC ("dbase_type = "Kindinger").
#'#'@param sites_associated Location data was collected. Choose either Oahu ("sites_associated = "OAH"),
#'or Mariana Islands ("sites_associated = "MARIAN").
#'
#'@import Rmisc
#'@import tidyverse
#'@import dplyr
#'
#'@export process_fish
#'
#'@examples
#'fish_data <- read.csv("CB_FishBelt_alldata.csv", na = "", check.names = FALSE)
#'
#'fish_belt <- process_fish(data = fish_data, dbase_type = "Kindinger",
#'sites_associated = "OAH")

calc_fish_strs_spc <- function(data,
                               dbase_type = c("IPRB", "Kindinger"),
                               sites_associated = c("OAH", "MARIAN")) {
  
  
  ifelse(dbase_type == "IPRB", rates_dbase <- fish_erosion_dbase_iprb, rates_dbase <- fish_erosion_dbase_kindinger)
  ifelse(sites_associated == "OAH", sites_associated_dbase <- fish_assoc_sites_oahu, sites_associated_dbase <- fish_assoc_sites_marian)
  ifelse(sites_associated == "OAH", loc <- "MHI", loc <- "MARIAN")
  
    
    # FOR StRS SPC DATA ----------------------------------------------------------------
    
    format_strsspc_output <- format_fish_spc(data,                                                 
                                             method = "nSPC",
                                             rates_dbase = rates_dbase)
    
    summary_strsspc_erosion <- format_strsspc_output %>%
      # convert REPLICATEID values to Transect '1' and '2'
      group_by(SITEVISITID) %>% 
      mutate(TRANSECT = match(REPLICATEID, unique(REPLICATEID))) %>%
      gather("METRIC", "VALUE", -c(SITEVISITID:FXN_GRP, TRANSECT)) %>% #set up data frame by spreading by transects (or the former REPLICATEID)
      select(-REPLICATEID) %>%
      spread(TRANSECT, VALUE) %>%
      mutate_at(vars(`1`, `2`), as.numeric) %>%
      # Average replicates (n=2)
      dplyr::group_by(SITEVISITID, SITE, REP, COMMONFAMILYALL, FXN_GRP, METRIC) %>% #calculate mean, sd, se, CI95 of data by averaging 2 divers (REPLICATEID)
      rowwise() %>%
      dplyr::mutate(MEAN = mean(c_across(`1`:`2`), na.rm = TRUE)) %>%
      dplyr::mutate(SD = sd(c_across(`1`:`2`), na.rm = TRUE)) %>%
      dplyr::mutate(SD = coalesce(SD, 0)) %>% # reaplce NA SD values with 0 because SD of one transect value is 0
      dplyr::mutate(SE = SD / sqrt(2)) %>% #2 is the number of total Transects
      dplyr::mutate(L95 = MEAN - (SE * 1.97)) %>%
      dplyr::mutate(U95 = MEAN + (SE * 1.97)) %>%
      select(-c(`1`:`2`)) %>% #remove unnecessary columns
      ungroup(.) %>%
      dplyr::mutate(L95 = case_when(L95 < 0 ~ 0,
                                    TRUE ~ as.numeric(L95)))
    
    format_strs_spc_erosion <- summary_strsspc_erosion %>%
      select(-c(SD:U95)) %>%
      #add ASSOC_OCCSITE before averaging by replicate
      left_join(., sites_associated_dbase %>% 
                  select(-LOCATION), by = "SITE") %>% # attach associated occ site ID 
      mutate(ASSOC_OCCSITE = case_when(value == 0 ~ "",
                                       TRUE ~ ASSOC_OCCSITE)) %>% # create assoc_occsite ID names...NOTE: the NA in assoc_occsite represents actual fixed site
      filter(!value %in% c("0", "1")) %>%
      filter(!is.na(value)) %>%
      select(-value) %>%
      # convert SITE values to Transect '1' through `n`
      group_by(ASSOC_OCCSITE) %>% 
      mutate(TRANSECT = match(SITE, unique(SITE))) %>%
      select(-SITEVISITID, -SITE, -REP)
    
    calc_strs_spc_erosion <- format_strs_spc_erosion %>%
      # create column for total transects per associated site ID because they vary
      left_join(., format_strs_spc_erosion %>% 
                  group_by(ASSOC_OCCSITE) %>%
                  summarise(n = n_distinct(TRANSECT)), 
                by = "ASSOC_OCCSITE") %>%      
      spread(TRANSECT, MEAN) %>%
      mutate_at(vars(n:length(.)), as.numeric) %>%
      ungroup(.) %>%
      # Average transects (n = variable) per associated site ID
      dplyr::group_by(ASSOC_OCCSITE, COMMONFAMILYALL, FXN_GRP, METRIC) %>% #calculate mean, sd, se, CI95 of data by averaging 2 divers (REPLICATEID)
      rowwise() %>%
      dplyr::mutate(MEAN = mean(c_across(`1`:(length(.)-5)), na.rm = TRUE)) %>%
      dplyr::mutate(SD = sd(c_across(`1`:(length(.)-6)), na.rm = TRUE)) %>%
      dplyr::mutate(SD = coalesce(SD, 0)) %>% # reaplce NA SD values with 0 because SD of one transect value is 0
      dplyr::mutate(SE = SD / sqrt(n)) %>% #2 is the number of total Transects
      dplyr::mutate(L95 = MEAN - (SE * 1.97)) %>%
      dplyr::mutate(U95 = MEAN + (SE * 1.97)) %>%
      ungroup(.) %>%
      select(c(COMMONFAMILYALL:ASSOC_OCCSITE, MEAN:U95)) %>% #remove unnecessary columns
      dplyr::mutate(L95 = case_when(L95 < 0 ~ 0,
                                    TRUE ~ as.numeric(L95))) %>%
      filter(!COMMONFAMILYALL %in% "NOTPARROTFISH") %>% # removed non-parrotfish species group
      filter(!ASSOC_OCCSITE %in% NA) %>% # remove NA sites...these correspond to the fixed site data as well
      filter(!FXN_GRP %in% NA) %>% # remove NA Graz_types (these are also all non-parrotfish species group)
      # Format dataframe to match BELT FINAL BIOEROSION OUTPUT
      gather(., "calc", "value", -c(COMMONFAMILYALL:ASSOC_OCCSITE)) %>%
      select(-COMMONFAMILYALL) %>% #remove n and REP
      mutate(METRIC = case_when(METRIC == "SUM_BIOMASS_PER_FISH_KG_HECTARE" ~ "FISH_BIOMASS_KG_HA",
                                METRIC == "SUM_DENSITY_PER_FISH_HECTARE" ~ "FISH_DENSITY_ABUNDANCE_HA",
                                METRIC == "SUM_EROSION_PER_FISH_KG_M2_YR" ~ "FISH_EROSION_KG_M2_YR")) %>%
      unite("METRIC", METRIC, FXN_GRP) %>% 
      unite("METRIC", METRIC, calc) %>%
      # fill in missing metadata info
      left_join(., data %>% select(REGION_NAME, MISSIONID, ISLAND, OCC_SITEID, LATITUDE, LONGITUDE, METHOD), by = c("ASSOC_OCCSITE" = "OCC_SITEID")) %>% #join important metadata that was lost during averaging replicates
      distinct(.) %>% #left_join creates duplicates, so remove duplicates
      mutate(LATITUDE = round(LATITUDE, 5)) %>%
      mutate(LONGITUDE = round(LONGITUDE, 5)) %>%
      rename(REGION = REGION_NAME,
             CRUISE_ID = MISSIONID,
             LOCATION = ISLAND,
             CB_METHOD = METHOD) %>%
      mutate(REGIONCODE = loc,
             LOCATIONCODE = str_sub(ASSOC_OCCSITE, 5,7)) %>%
      rename(OCC_SITEID = ASSOC_OCCSITE) %>%
      left_join(., data %>% select(OCC_SITEID, OCC_SITENAME), by = "OCC_SITEID") %>%
      distinct(.) %>%
      mutate(CB_METHOD = "StRS SPC") %>%
      spread(., METRIC, value, fill = "0") %>%
      select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID, OCC_SITENAME, LATITUDE, LONGITUDE, CB_METHOD, everything(.)) %>%
      filter(!is.na(.$REGION)) %>%
      mutate_at(vars(REGION:CB_METHOD), as.factor) %>%
      mutate_at(vars(FISH_BIOMASS_KG_HA_ALL_L95:FISH_EROSION_KG_M2_YR_SCRAPER_U95), as.numeric)
    
    return(calc_strs_spc_erosion)
    
  
}

