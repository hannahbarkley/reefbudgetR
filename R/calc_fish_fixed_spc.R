#' Calculate parrotfish erosion rates from fixed spc data
#'
#'@author Rebecca Weible
#'
#'@param data fish survey data.
#'@param dbase_type Erosion rates database to use. Choose either Indo-Pacific
#'ReefBudget ("dbase_type = "IPRB") or U.S. Pacific Islands rates developed
#'by Tye Kindinger, NOAA PIFSC ("dbase_type = "Kindinger").
#'
#'@import Rmisc
#'@import sf
#'@import rgeos
#'@import tidyverse
#'@import dplyr
#'
#'@export calc_fish_fixed_spc
#'
#'@examples
#'fish_data <- read.csv("CB_FishBelt_alldata.csv", na = "", check.names = FALSE)
#'
#'fish_fixed_spc <- calc_fish_fixed_spc(data = fish_data, dbase_type = "Kindinger")

calc_fish_fixed_spc <- function(data,
                                dbase_type = c("IPRB", "Kindinger")) {
  
  
  ifelse(dbase_type == "IPRB", rates_dbase <- fish_erosion_dbase_iprb, rates_dbase <- fish_erosion_dbase_kindinger)
  
  
    # FOR FIXED SPC DATA ----------------------------------------------------------------
    
    format_spc_output <- format_fish_spc(data,                                                 
                                         method = "IPRB",
                                         rates_dbase = rates_dbase)
  
    # created associated SPC sites to each OCC fixed site SPC
    sites_associated_dbase <- create_fish_assoc_sites(data)
    
    summary_spc_erosion <- format_spc_output %>%
      # convert REPLICATEID values to Transect '1' and '2'
      group_by(SITEVISITID) %>% 
      mutate(TRANSECT = match(REPLICATEID, unique(REPLICATEID))) %>%
      gather("METRIC", "VALUE", -c(SITEVISITID:FXN_GRP, TRANSECT)) %>% #set up data frame by spreading by transects (or the former REPLICATEID)
      select(-REPLICATEID) %>%
      spread(TRANSECT, VALUE) %>%
      mutate_at(vars(`1`, `2`), as.numeric) %>%
      # Average replicates (n=2)
      dplyr::group_by(SITEVISITID, REA_SITEID, REP, COMMONFAMILYALL, FXN_GRP, METRIC) %>% #calculate mean, sd, se, CI95 of data by averaging 2 divers (REPLICATEID)
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
      filter(!REA_SITEID %in% NA) %>% # remove NA sites
      filter(!FXN_GRP %in% NA) %>% # remove NA Functional Groups (these are also all non-parrotfish species group)
      
      # Format dataframe to match BELT FINAL BIOEROSION OUTPUT
      gather(., "calc", "value", -c(SITEVISITID:METRIC)) %>%
      select(-COMMONFAMILYALL, -REP) %>% #remove n and REP
      mutate(METRIC = case_when(METRIC == "SUM_BIOMASS_PER_FISH_KG_HECTARE" ~ "FISH_BIOMASS_KG_HA",
                                METRIC == "SUM_DENSITY_PER_FISH_HECTARE" ~ "FISH_DENSITY_ABUNDANCE_HA",
                                METRIC == "SUM_EROSION_PER_FISH_KG_M2_YR" ~ "FISH_EROSION_KG_M2_YR")) 
    
    # pause here to complete FXN_GRP column if values are missing
    options_fxn_grp <- c("ALL", "OTHER", "SCRAPER", "EXCAVATOR") # all four functional groups should be present
    missing_fxn_grp <- options_fxn_grp[!options_fxn_grp %in% unique(format_fixed_spc_erosion$FXN_GRP)] #isolate missing functional groups
      
    # Add back in rows of missing functional groups and assign value to 0 with completing the dataframe
    if(length(missing_fxn_grp) >= 1) {
      for(i in 1:length(missing_fxn_grp)) {
        
        format_fixed_spc_erosion <- add_row(format_fixed_spc_erosion, FXN_GRP = missing_fxn_grp[i])
        
      }
      
      format_fixed_spc_erosion2 <- format_fixed_spc_erosion %>%
                                      select(-SITEVISITID) %>%
                                      complete(REA_SITEID, METRIC, calc, FXN_GRP,
                                        fill = list(value = 0)
                                      ) %>%
                                      filter(!METRIC %in% NA) %>%
                                      filter(!calc %in% NA) %>%
                                      filter_all(all_vars(!is.na(.))) %>%
                                      distinct(.)
      
      } else {format_fixed_spc_erosion2 <- format_fixed_spc_erosion}
      
    # Continue formatting dataframe to match BELT FINAL BIOEROSION OUTPUT
    format_fixed_spc_erosion3 <-   
      format_fixed_spc_erosion2 %>%
      unite("METRIC", METRIC,FXN_GRP) %>% 
      unite("METRIC", METRIC:calc) %>%
      # fill in missing metadata info
      left_join(., data %>% select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, REA_SITEID, LATITUDE, LONGITUDE, CB_METHOD), by = "REA_SITEID") %>% #join important metadata that was lost during averaging replicates
      distinct(.) %>% #left_join creates duplicates, so remove duplicates
      mutate(LATITUDE = round(LATITUDE, 5)) %>%
      mutate(LONGITUDE = round(LONGITUDE, 5)) %>%
      #mutate(REGIONCODE = sites_associated_dbase$REGIONCODE[1],
      #       LOCATIONCODE = str_extract(REA_SITEID, "(\\w+)")) %>%
      #add OCC_SITEID column
      left_join(., data %>% select(REA_SITEID, OCC_SITEID, OCC_SITENAME), by = "REA_SITEID") %>%
      distinct(.) %>%
      spread(., METRIC, value, fill = 0) %>%
      select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID, OCC_SITENAME, LATITUDE, LONGITUDE, CB_METHOD, everything(.), -REA_SITEID)
    
    
    summary_fixed_spc_erosion <- format_fixed_spc_erosion3 %>%
      filter(!is.na(OCC_SITEID)) %>% # remove all non fixed SPC site data
      mutate(CB_METHOD = "Fixed SPC") %>%
      mutate_at(vars(REGION:CB_METHOD), as.factor) %>%
      mutate_at(vars(FISH_BIOMASS_KG_HA_ALL_L95:FISH_EROSION_KG_M2_YR_SCRAPER_U95), as.numeric)
    
    return(summary_fixed_spc_erosion)
  
  
  
 }