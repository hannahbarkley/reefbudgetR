#' Calculate parrotfish erosion rates from belt data
#'
#'@author Rebecca Weible
#'
#'@param data fish survey data.
#'@param method survey method. Choose either Fixed Belt ("method = "Fixed Belt"), 
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
                         method = c("Fixed Belt", "Fixed SPC", "StRS SPC"),
                         rates_dbase = c("IPRB", "Kindinger"),
                         sites_associated = c("OAH", "MARIAN"),
                         full_summary = TRUE) {

  ifelse(sites_associated == "OAH", sites_associated <- fish_assoc_occ_oahu, sites_associated <- fish_assoc_occ_marian)
  
  if (method == "Fixed Belt") {  
        
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
      group_by(SITEVISITID) %>% # convert REPLICATEID values to Transect '1' and '2'
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
                                    TRUE ~ as.numeric(L95))) %>%
      filter(!COMMONFAMILYALL %in% "NOTPARROTFISH") %>% # removed non-parrotfish species group
      filter(!SITE %in% NA) %>% # remove NA sites
      filter(!FXN_GRP %in% NA) %>% # remove NA Functional Groups (these are also all non-parrotfish species group)
      
      # Format dataframe to match BELT FINAL BIOEROSION OUTPUT
      gather(., "calc", "value", -c(SITEVISITID:METRIC)) %>%
      select(-COMMONFAMILYALL, -REP) %>% #remove n and REP
      mutate(METRIC = case_when(METRIC == "SUM_BIOMASS_PER_FISH_KG_HECTARE" ~ "FISH_BIOMASS_KG_HA",
                                METRIC == "SUM_DESNITY_PER_FISH_HECTARE" ~ "FISH_DENSITY_ABUNDANCE_HA",
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
      mutate(REGIONCODE = "MARIAN",
             LOCATIONCODE = str_extract(SITE, "(\\w+)")) %>%
      #add OCC_SITEID column
      left_join(., sites_associated %>% 
                  filter(value == "1") %>% # select only fixed OCC sites (do not include associated)
                  select(ASSOC_OCCSITE, SITE, LOCATION) %>%
                  rename(OCC_SITENAME = LOCATION), 
                by = "SITE") %>%
      rename(OCC_SITEID = ASSOC_OCCSITE) %>%
      spread(., METRIC, value, fill = 0) %>%
      select(REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID, OCC_SITENAME, LATITUDE, LONGITUDE, CB_METHOD, everything(.), -SITEVISITID, -SITE)
    
    
    
    summary_fixed_spc_erosion <- summary_spc_erosion %>%
                                  filter(!is.na(OCC_SITEID)) %>% # remove all non fixed SPC site data
                                  mutate(CB_METHOD = "Fixed SPC")
    
    return(summary_fixed_spc_erosion)
  }


}

