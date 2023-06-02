#' Calculate parrotfish erosion rates from belt data
#'
#'@author Rebecca Weible
#'
#'@param data fish survey data.
#'@param method survey method. Choose either Fixed Belt ("method = "IPRB"), 
#'Fixed Stationary Point Count ("method = "Fixed SPC"), or Associated Stationary 
#'Point Count ("method = "StRS SPC").
#'@param rates_dbase Erosion rates database to use. Choose either Indo-Pacific
#'ReefBudget ("rates_dbase = "IPRB") or U.S. Pacific Islands rates developed
#'by Tye Kindinger, NOAA PIFSC ("rates_dbase = "Kindinger").
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
#'fish_belt <- process_fish(data = fish_data, method = "Fixed SPC", rates_dbase = "Kindinger",
#'full_summary = TRUE)

process_fish <- function(data,
                         method = "IPRB",
                         rates_dbase = c("IPRB", "Kindinger"),
                         sites_associated = "OAH",
                         full_summary = TRUE) {

  ifelse(sites_associated == "OAH", sites_associated_dbase <- fish_assoc_sites_oahu, sites_associated_dbase <- fish_assoc_sites_marian)
  ifelse(rates_dbase == "IPRB", rates_dbase <- fish_erosion_dbase_iprb, rates_dbase <- fish_erosion_dbase_kindinger)
  
  if (method == "IPRB") {  
        
        # Calculate erosion rates per fish -------------------------------------------
      
        calc_eros_fish_output <- calc_eros_fish(data,
                                                rates_dbase = rates_dbase)
      
        # Calculate bioerosion metrics per grazing type per site ---------------------
      
        density_average <- suppressWarnings(
          summarize_fish_metrics(
            data = calc_eros_fish_output,
            metric = "density",
            level = "transect",
            summarize_by = "species"
          )
        ) %>%
          add_column(METRIC = "FISH_DENSITY_ABUNDANCE_HA")
      
        biomass_average <- suppressWarnings(
          summarize_fish_metrics(
            data = calc_eros_fish_output,
            metric = "biomass",
            level = "transect",
            summarize_by = "species"
          )
        ) %>%
          add_column(METRIC = "FISH_BIOMASS_KG_HA")
      
        bioerosion_average <- suppressWarnings(
          summarize_fish_metrics(
            data = calc_eros_fish_output,
            metric = "bioerosion",
            level = "transect",
            summarize_by = "species"
          )
        ) %>%
          add_column(METRIC = "FISH_EROSION_KG_M2_YR")
      
      
        species_table <-
          rbind(density_average, biomass_average, bioerosion_average)
      
      
        summary_belt_erosion <- summarize_fish_erosion(species_table, full_summary)
      
        return(summary_belt_erosion)
      
        if (full_summary == TRUE) {
          return(list(
            fish_erosion_transect = summary_belt_erosion$fish_erosion_transect,
            fish_erosion_site = summary_belt_erosion$fish_erosion_site)
          )
        }
      
        if (full_summary == FALSE) {
          return(summary_belt_erosion$fish_erosion_site)
        }
  }
  
  
  if (method == "Fixed SPC") {  
    
    format_spc_output <- format_fish_spc(data, rates_dbase = "Kindinger")
    
    summary_spc_erosion <- format_spc_output %>%
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
    
    format_fixed_spc_erosion <- summary_spc_erosion %>%
      filter(!COMMONFAMILYALL %in% "NOTPARROTFISH") %>% # removed non-parrotfish species group
      filter(!SITE %in% NA) %>% # remove NA sites
      filter(!FXN_GRP %in% NA) %>% # remove NA Functional Groups (these are also all non-parrotfish species group)
      
      # Format dataframe to match BELT FINAL BIOEROSION OUTPUT
      gather(., "calc", "value", -c(SITEVISITID:METRIC)) %>%
      select(-COMMONFAMILYALL, -REP) %>% #remove n and REP
      mutate(METRIC = case_when(METRIC == "SUM_BIOMASS_PER_FISH_KG_HECTARE" ~ "FISH_BIOMASS_KG_HA",
                                METRIC == "SUM_DENSITY_PER_FISH_HECTARE" ~ "FISH_DENSITY_ABUNDANCE_HA",
                                METRIC == "SUM_EROSION_PER_FISH_KG_M2_YR" ~ "FISH_EROSION_KG_M2_YR")) %>%
      unite("METRIC", METRIC:FXN_GRP) %>% 
      unite("METRIC", METRIC:calc) %>%
      # fill in missing metadata info
      left_join(., data %>% select(REGION_NAME, MISSIONID, ISLAND, SITEVISITID, SITE, LATITUDE, LONGITUDE, METHOD), by = c("SITEVISITID", "SITE")) %>% #join important metadata that was lost during averaging replicates
      distinct(.) %>% #left_join creates duplicates, so remove duplicates
      mutate(LATITUDE = round(LATITUDE, 5)) %>%
      mutate(LONGITUDE = round(LONGITUDE, 5)) %>%
      rename(REGION = REGION_NAME,
             CRUISE_ID = MISSIONID,
             LOCATION = ISLAND,
             CB_METHOD = METHOD) %>%
      mutate(REGIONCODE = sites_associated,
             LOCATIONCODE = str_extract(SITE, "(\\w+)")) %>%
      #add OCC_SITEID column
      left_join(., sites_associated_dbase %>% 
                  filter(value == "1") %>% # select only fixed OCC sites (do not include associated)
                  select(ASSOC_OCCSITE, SITE, LOCATION) %>%
                  rename(OCC_SITENAME = LOCATION), 
                by = "SITE") %>%
      rename(OCC_SITEID = ASSOC_OCCSITE) %>%
      spread(., METRIC, value, fill = 0) %>%
      select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID, OCC_SITENAME, LATITUDE, LONGITUDE, CB_METHOD, everything(.), -SITEVISITID, -SITE)
    
    
    
    summary_fixed_spc_erosion <- format_fixed_spc_erosion %>%
                                  filter(!is.na(OCC_SITEID)) %>% # remove all non fixed SPC site data
                                  mutate(CB_METHOD = "Fixed SPC")
    
    return(summary_fixed_spc_erosion)
    }
    
    
    
  if (method == "StRS SPC") {  
    
    format_spc_output <- format_fish_spc(data, rates_dbase = "Kindinger")
    
    summary_spc_erosion <- format_spc_output %>%
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
    
    format_strs_spc_erosion <- summary_spc_erosion %>%
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
      left_join(., summary_fixed_spc_erosion %>% select(c(REGION:LONGITUDE)), 
                by = c("ASSOC_OCCSITE" = "OCC_SITEID")) %>% #join important metadata that was lost during averaging replicates
      distinct(.) %>% #left_join creates duplicates, so remove duplicates
      rename(OCC_SITEID = ASSOC_OCCSITE) %>%
      mutate(CB_METHOD = "StRS SPC") %>%
      spread(., METRIC, value, fill = "0") %>%
      select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID, OCC_SITENAME, LATITUDE, LONGITUDE, CB_METHOD, everything(.))
      
    return(calc_strs_spc_erosion)
      
    }
}


