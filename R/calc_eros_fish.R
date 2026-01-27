#' Calculate parrotfish biomass, density, and bioerosion rates

#'@author Rebecca Weible
#'
#'@param data Parrotfish belt data, including number of fish observed of each
#'species, size class, and phase.
#'@param rates_dbase Erosion rates database to use.
#'
#'@import dplyr
#'@importFrom rlang .data
#'
#'@export calc_eros_fish
#'
#'@examples fish_data <- read.csv("CB_FishBelt_alldata.csv", na = "", check.names = FALSE)
#'@examples calc_eros_fish_output <- calc_eros_fish(fish_data, rates_dbase = rates_dbase)


calc_eros_fish <- function(data, rates_dbase) {
  
  # Setup Lookup Tables -----------------------------
  size_to_length <- c(
    "SIZE_BIN_0_10_CM"  = 7,
    "SIZE_BIN_11_20_CM" = 15,
    "SIZE_BIN_21_30_CM" = 25,
    "SIZE_BIN_31_40_CM" = 35,
    "SIZE_BIN_41_50_CM" = 45,
    "SIZE_BIN_51_60_CM" = 55
  )
  
  rates_size_map <- c(
    "0-10cm"  = "SIZE_BIN_0_10_CM",
    "11-20cm" = "SIZE_BIN_11_20_CM",
    "21-30cm" = "SIZE_BIN_21_30_CM",
    "31-40cm" = "SIZE_BIN_31_40_CM",
    "41-50cm" = "SIZE_BIN_41_50_CM",
    "51-60cm" = "SIZE_BIN_51_60_CM"
  )
  
  # Create Reference Databases ----------------------------------------------------
  
  rates_db <- rates_dbase %>%
    dplyr::rename(any_of(c(TAXONNAME = "TAXON_NAME", SPECIES = "TAXON_CODE"))) %>%
    dplyr::mutate(SIZE_CLASS = rates_size_map[SIZE_CLASS]) %>%
    dplyr::filter(!is.na(SIZE_CLASS)) %>%
    dplyr::select(TAXONNAME, SIZE_CLASS, PHASE, EROSION_RATE, FXN_GRP_DB = FXN_GRP)
  
  species_db <- fish_species_dbase %>%
    dplyr::select(SPECIES, LW_A, LW_B, LENGTH_CONVERSION_FACTOR) %>%
    dplyr::distinct()
  
  # Calculation Pipeline -----------------------------------------------
  
  fish_all <- data %>%
    dplyr::rename(any_of(c(METHOD = "CB_METHOD", TAXONNAME = "TAXON_NAME",SPECIES = "TAXON_CODE"))) %>%
    dplyr::mutate(
      dplyr::across(c(REGION:LONGITUDE, METHOD:VISIBILITY_M, TRANSECT, 
                      HABITAT_CODE:TAXONNAME, PHASE), as.factor),
      dplyr::across(c(TRANSECT_LENGTH_M:AREA_M2, starts_with("SIZE_BIN_")), as.numeric)
    ) %>%
    dplyr::select(
      REGION:TRANSECT, 
      SPECIES:PHASE, 
      starts_with("SIZE_BIN_")
    ) %>%
    tidyr::pivot_longer(
      cols = starts_with("SIZE_BIN_"),
      names_to = "SIZE_CLASS",
      values_to = "COUNT",
      values_drop_na = TRUE
    ) %>%
    dplyr::filter(
      !(PHASE == "J" & SIZE_CLASS != "SIZE_BIN_0_10_CM"),
      !(PHASE %in% c("I", "T") & SIZE_CLASS == "SIZE_BIN_0_10_CM"),
      !(PHASE == "I" & SIZE_CLASS == "SIZE_BIN_51_60_CM")
    ) %>%
    dplyr::mutate(length = size_to_length[SIZE_CLASS]) %>%
    dplyr::left_join(species_db, by = "SPECIES") %>%
    dplyr::left_join(rates_db, by = c("TAXONNAME", "SIZE_CLASS", "PHASE")) %>%
    dplyr::mutate(
      # --- Biomass ---
      length = dplyr::coalesce(length, 0),
      BIOMASS_PER_FISH_G = LW_A * ((length * LENGTH_CONVERSION_FACTOR) ^ LW_B),
      BIOMASS_PER_FISH_KG = (COUNT * BIOMASS_PER_FISH_G) / 1000,
      BIOMASS_PER_FISH_KG = dplyr::coalesce(BIOMASS_PER_FISH_KG, 0),
      # --- Bioerosion ---
      FXN_GRP = FXN_GRP_DB,
      FXN_GRP = dplyr::case_when(
        SPECIES == "PARR" & SIZE_CLASS != "SIZE_BIN_0_10_CM" ~ "Scraper",
        SIZE_CLASS == "SIZE_BIN_0_10_CM" ~ "Browser",
        is.na(FXN_GRP) ~ "Scraper",
        TRUE ~ FXN_GRP
      ),
      FXN_GRP = dplyr::if_else(!FXN_GRP %in% c("Excavator", "Scraper"), "Other", FXN_GRP),
      
      EROSION_RATE = dplyr::coalesce(as.numeric(EROSION_RATE), 0),
      FISH_EROSION_KG_M2_YR = COUNT * EROSION_RATE,
      FISH_EROSION_KG_M2_YR = pmax(0, FISH_EROSION_KG_M2_YR)
    ) %>%
    dplyr::select(-any_of("FXN_GRP_DB"))
  
  return(fish_all)
}
