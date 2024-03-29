#' Format Stationary Point Count data
#'
#'@author Rebecca Weible
#'
#'@param data all stationary point count data.
#'@param method type of SPC survey conducted. Choose either Fixed site SPC 
#'("method = "CbB") or stratified random sampling SPC ("method = "nSPC").
#'@param rates_dbase Erosion rates database to use.
#'
#'@import tidyverse
#'@import dplyr
#'
#'@export format_fish_spc
#'


format_fish_spc <- function(data, 
                            method = c("IPRB", "nSPC"),
                            rates_dbase = rates_dbase) {
  
  if(method == "IPRB") {method_type = "IPRB"} else {method_type = "nSPC"}
  
  prepdat <- data %>% 
    dplyr::filter(., CB_METHOD %in% method_type) %>%
    dplyr::filter(., !(TRAINING_YN %in% "-1")) %>%
    dplyr::filter(., !(REA_SITEID == "GUA-2587" & CB_METHOD == "nSPC")) %>% # must remove special case "GUA-2587" SPC because it's a fixed site and data was collected twice.
    dplyr::select(SITEVISITID, YEAR, REA_SITEID, LOCATION, REP, REPLICATEID, TAXON_CODE, COUNT, SIZE_, TAXON_NAME, COMMONFAMILYALL, LW_A:LENGTH_CONVERSION_FACTOR) %>% #subset only columns that matter for density, biomass, bioerosion calculation
    mutate(AREA_M2 = pi*(7.5^2)) %>% #calculate Area per m^2 of survey cylinder
    #calculate density below
    mutate(DENSITY_PER_FISH_HECTARE = COUNT/(AREA_M2/10000)) %>% # calculate density by dividing count by area and converting to per hectare with the 10000
    #calculate biomass below
    mutate(BIOMASS_PER_FISH_G = LW_A * ((SIZE_ * LENGTH_CONVERSION_FACTOR) ^ LW_B)) %>% #convert size to wet weight in grams with conversion factor and constants
    mutate(BIOMASS_PER_FISH_KG = (COUNT * BIOMASS_PER_FISH_G) / 1000) %>% #calculate biomass per row with number of fish seen and converted grams to kilograms with /1000
    mutate(BIOMASS_PER_FISH_KG_HECTARE = BIOMASS_PER_FISH_KG/(AREA_M2 / 10000)) %>% #calculate final biomass value of each row by dividing by area and converting m^2 to hectares
    #calculate bioerosion below
    mutate(SIZE_CLASS = case_when(SIZE_ >= 0 & SIZE_ <= 10 ~ "0-10cm", # need to create column with size bin range for joining Tye's bioerosion metrics
                                  SIZE_ >= 11 & SIZE_ <= 20 ~ "11-20cm",
                                  SIZE_ >= 21 & SIZE_ <= 30 ~ "21-30cm",
                                  SIZE_ >= 31 & SIZE_ <= 40 ~ "31-40cm",
                                  SIZE_ >= 41 & SIZE_ <= 50 ~ "41-50cm",
                                  SIZE_ >= 51 & SIZE_ <= 60 ~ "51-60cm")) %>% 
    left_join(., rates_dbase %>% dplyr::select(-PHASE), by = c("TAXON_NAME", "SIZE_CLASS"), relationship = "many-to-many") %>% # join Tye's bioerosion metrics or bioerosion calculation to follow
    distinct(.) %>%
    mutate(EROSION_PER_FISH_KG_M2_YR = (COUNT * EROSION_RATE)/AREA_M2) %>% # calculate bioerosion in kg/m^2/yr
    mutate_at(.vars = "EROSION_PER_FISH_KG_M2_YR", funs(ifelse(EROSION_PER_FISH_KG_M2_YR <= 0, 0, .)))  %>% # change all negative bioerosion values to zero...can use this to change multiple columns to zero based on a single column
    
    # Label COMMONFAMILYALL column as "Parrotfish" and "NotParrotfish"        
    dplyr::select(SITEVISITID:COMMONFAMILYALL, SIZE_CLASS, FXN_GRP, DENSITY_PER_FISH_HECTARE, BIOMASS_PER_FISH_KG_HECTARE, EROSION_PER_FISH_KG_M2_YR) %>% # select only relevant columns
    mutate(COMMONFAMILYALL = case_when(COMMONFAMILYALL != "Parrotfish" ~ "NOTPARROTFISH", # assign non parrotfish species NOTPARROTFISH
                                       TRUE ~ "Parrotfish")) %>%
    
    # Sum FXN_GRP to the REPLICATEID level
    mutate(FXN_GRP = replace(as.character(FXN_GRP), FXN_GRP == "Browser", "Other")) %>%
    dplyr::group_by(SITEVISITID, REA_SITEID, REP, REPLICATEID, COMMONFAMILYALL, FXN_GRP) %>% #now we want to sum each surveyors estimates by Grazing Type (note, REP is not divided by surveyor, but REPLICATEID is)
    dplyr::summarize("SUM_BIOMASS_PER_FISH_KG_HECTARE" = sum(BIOMASS_PER_FISH_KG_HECTARE),
                     "SUM_DENSITY_PER_FISH_HECTARE" = sum(DENSITY_PER_FISH_HECTARE),
                     "SUM_EROSION_PER_FISH_KG_M2_YR" = sum(EROSION_PER_FISH_KG_M2_YR)) %>%
    ungroup(.) %>% # need to do this to make next group_by work properly
    bind_rows(dplyr::filter(.) %>% # create new rows where Excavator, Scraper, and Other are totaled in sum for each site and replicate ID
                group_by(SITEVISITID, REA_SITEID, REP, REPLICATEID, COMMONFAMILYALL) %>%
                summarise(across(SUM_BIOMASS_PER_FISH_KG_HECTARE:SUM_EROSION_PER_FISH_KG_M2_YR, ~(sum(.x, na.rm=T)))) %>% 
                mutate(FXN_GRP = "ALL")) %>%
    #format by adding back in where replicateID and Graz_Type were zero before averaging
    complete(., REPLICATEID, COMMONFAMILYALL, FXN_GRP, fill = list(SUM_BIOMASS_PER_FISH_KG_HECTARE = 0, 
                                                                   SUM_DENSITY_PER_FISH_HECTARE = 0, 
                                                                   SUM_EROSION_PER_FISH_KG_M2_YR = 0)) %>% #...I can complete (fill in missing) SPECIES codes for each REPLICATEID. So we're adding back zeros and this is important becuase of the mean calculations we're about to do.
    dplyr::select(-SITEVISITID, -REA_SITEID, -REP) %>% # remove columns with NA and will join them back in next step
    mutate_at(vars(REPLICATEID), as.integer) %>% # make structure of REPLICATEID the same again for the sake of the join
    left_join(., data %>% dplyr::select(SITEVISITID, REA_SITEID, REP, REPLICATEID), by = "REPLICATEID", relationship = "many-to-many") %>% # join the three columns that produced NAs with function complete
    distinct(.) %>% # remove duplicate rows
    dplyr::select(SITEVISITID, REA_SITEID, REP, everything(.)) %>% # re-order columns for visual effects
    mutate(FXN_GRP = case_when(FXN_GRP == "Other" ~ "OTHER",
                               FXN_GRP == "Scraper" ~ "SCRAPER",
                               FXN_GRP == "Excavator" ~ "EXCAVATOR",
                               TRUE ~ FXN_GRP))
  
  
  return(prepdat)
  
}
