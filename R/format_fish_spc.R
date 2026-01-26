#' Format Stationary Point Count data
#'
#'@author Rebecca Weible
#'
#'@param data all stationary point count data.
#'@param method type of SPC survey conducted. Choose either Fixed site SPC 
#'("method = "fSPC") or stratified random sampling SPC ("method = "nSPC").
#'@param rates_dbase Erosion rates database to use.
#'
#'@import tidyverse
#'@import dplyr
#'
#'@export format_fish_spc


format_fish_spc <- function(data, 
                            method = c("fSPC", "nSPC"),
                            rates_dbase = rates_dbase) {
  
  method_arg <- match.arg(method) 
  method_type <- if(method_arg == "fSPC" || method_arg == "CbB") "fSPC" else "nSPC"
  
  cyl_area_m2 <- pi * (7.5^2)
  cyl_area_ha <- cyl_area_m2 / 10000
  
  formatdat <- data %>% 
    
    dplyr::filter(METHOD %in% method_type, TRAINING_YN != "-1") %>%
    
    dplyr::select(SITEVISITID, YEAR, REA_SITEID, LOCATION, REP, REPLICATEID, 
                  SPECIES, COUNT, SIZE_, TAXONNAME, COMMONFAMILYALL, 
                  LW_A:LENGTH_CONVERSION_FACTOR) %>%
    
    mutate(
      AREA_M2 = cyl_area_m2,
      
      # Density
      DENSITY_PER_FISH_HECTARE = COUNT / cyl_area_ha,
      
      # Calculate grams, then kg, then kg/ha 
      BIOMASS_PER_FISH_G = LW_A * ((SIZE_ * LENGTH_CONVERSION_FACTOR) ^ LW_B),
      BIOMASS_PER_FISH_KG = (COUNT * BIOMASS_PER_FISH_G) / 1000,
      BIOMASS_PER_FISH_KG_HECTARE = BIOMASS_PER_FISH_KG / cyl_area_ha,
      
      SIZE_CLASS = cut(SIZE_, 
                       breaks = c(0, 10, 20, 30, 40, 50, Inf), 
                       labels = c("0-10cm", "11-20cm", "21-30cm", "31-40cm", "41-50cm", "51-60cm"),
                       include.lowest = TRUE, right = TRUE)
    ) %>%
    
    # Join Rates
    left_join(rates_dbase, by = c("TAXONNAME", "SIZE_CLASS"), relationship = "many-to-many") %>%
    distinct() %>%
    
    # Final Bioerosion Math
    mutate(
      EROSION_RATE = coalesce(EROSION_RATE, 0), 
      EROSION_PER_FISH_KG_M2_YR = pmax(0, (COUNT * EROSION_RATE) / cyl_area_m2)
    )
  
  if ("PHASE" %in% names(formatdat)) {
    prep_step <- formatdat %>% filter(PHASE %in% c("I", "T"))
  } else {
    prep_step <- formatdat
  }
  
  prepdat <- prep_step %>%

    mutate(
      COMMONFAMILYALL = if_else(COMMONFAMILYALL == "Parrotfish", "Parrotfish", "NOTPARROTFISH"),
      FXN_GRP = case_when(
        FXN_GRP == "Browser"   ~ "OTHER",
        FXN_GRP == "Scraper"   ~ "SCRAPER",
        FXN_GRP == "Excavator" ~ "EXCAVATOR",
        TRUE ~ "OTHER"
      )
    ) %>%
    

    group_by(SITEVISITID, REA_SITEID, REP, REPLICATEID, COMMONFAMILYALL, FXN_GRP) %>%
    summarize(
      SUM_BIOMASS_PER_FISH_KG_HECTARE = sum(BIOMASS_PER_FISH_KG_HECTARE, na.rm = TRUE),
      SUM_DENSITY_PER_FISH_HECTARE    = sum(DENSITY_PER_FISH_HECTARE, na.rm = TRUE),
      SUM_EROSION_PER_FISH_KG_M2_YR   = sum(EROSION_PER_FISH_KG_M2_YR, na.rm = TRUE),
      .groups = "drop"
    )
  

  sum_all <- prepdat %>%
    group_by(SITEVISITID, REA_SITEID, REP, REPLICATEID, COMMONFAMILYALL) %>%
    summarize(
      SUM_BIOMASS_PER_FISH_KG_HECTARE = sum(SUM_BIOMASS_PER_FISH_KG_HECTARE, na.rm = TRUE),
      SUM_DENSITY_PER_FISH_HECTARE    = sum(SUM_DENSITY_PER_FISH_HECTARE, na.rm = TRUE),
      SUM_EROSION_PER_FISH_KG_M2_YR   = sum(SUM_EROSION_PER_FISH_KG_M2_YR, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(FXN_GRP = "ALL")
  
  prepdat <- bind_rows(prepdat, sum_all) %>%
    complete(
      nesting(SITEVISITID, REA_SITEID, REP, REPLICATEID), 
      COMMONFAMILYALL, 
      FXN_GRP, 
      fill = list(
        SUM_BIOMASS_PER_FISH_KG_HECTARE = 0, 
        SUM_DENSITY_PER_FISH_HECTARE = 0, 
        SUM_EROSION_PER_FISH_KG_M2_YR = 0
      )
    ) %>%
    dplyr::select(SITEVISITID, REA_SITEID, REP, REPLICATEID, COMMONFAMILYALL, FXN_GRP,
                  SUM_BIOMASS_PER_FISH_KG_HECTARE, 
                  SUM_DENSITY_PER_FISH_HECTARE, 
                  SUM_EROSION_PER_FISH_KG_M2_YR)
  
  
  return(prepdat)
  
}


