#' Calculate parrotfish biomass, density, and bioerosion rates

#'@author Rebecca Weible
#'
#'@param data Parrotfish belt data, including number of fish observed of each
#'species, size class, and phase.
#'@param dbase_types Erosion rates database to use. Choose either Indo-Pacific
#'ReefBudget ("dbase_type = "IPRB") or U.S. Pacific Islands rates developed
#'by Tye Kindinger, NOAA PIFSC ("dbase_type = "Kindinger").
#'
#'@import dplyr
#'@importFrom rlang .data
#'
#'@export calc_eros_fish
#'
#'@examples fish_data <- read.csv("CB_FishBelt_alldata.csv", na = "", check.names = FALSE)
#'@examples calc_eros_fish_output <- calc_eros_fish(fish_data, dbase_type = "Kindinger")


calc_eros_fish <- function(data,
                           dbase_types = c("IPRB", "Kindinger")) {
  
  ifelse(dbase_types == "Kindinger", rates_dbase <- fish_erosion_dbase_kindinger, rates_dbase <- fish_erosion_dbase_iprb)
  
  # Format dataframe for biomass and bioerosion calculations below
  data_formatted <- data %>%
    #change characters to factors
    dplyr::mutate_at(vars(
      c(REGION:LONGITUDE,
        CB_METHOD:VISIBILITY_M,
        TRANSECT,
        HABITAT_CODE:TAXON_NAME,
        PHASE
      )
    ), as.factor) %>%
    # convert integers to numbers
    mutate_at(vars(c(
      TRANSECT_LENGTH_M:AREA_M2,
      "SIZE_BIN_0_10_CM":"SIZE_BIN_51_60_CM"
    )), as.numeric) %>%
    # select only the important columns
    select(REGION:TRANSECT,
           TAXON_CODE:PHASE) %>%
    #column names to row values
    gather(.,
           "SIZE_CLASS",
           "COUNT", -c(REGION:TAXON_NAME),
           -PHASE, na.rm = TRUE) %>%
    # make size bins into factors
    mutate_at(vars(SIZE_CLASS), as.factor) %>%
    #remove all size classes above 10cm for J classification
    filter(!(
      SIZE_CLASS %in% c("SIZE_BIN_11_20_CM",
                        "SIZE_BIN_21_30_CM",
                        "SIZE_BIN_31_40_CM",
                        "SIZE_BIN_41_50_CM",
                        "SIZE_BIN_51_60_CM") &
        PHASE == "J"
    )) %>%
    # remove all I and T for 0-10cm classification
    filter(!(SIZE_CLASS %in% "SIZE_BIN_0_10_CM" &
               PHASE %in% c("I", "T"))) %>%
    # remove all sizes above 51cm for initial phases
    filter(!(SIZE_CLASS %in% "SIZE_BIN_51_60_CM" &
               PHASE == "I"))
  
  # rates_dbase have different size class name, so need to match them with our input data
  rates_dbase_formatted <- rates_dbase %>%
    mutate(
      SIZE_CLASS = case_when(
        SIZE_CLASS == "0-10cm" ~ "SIZE_BIN_0_10_CM",
        SIZE_CLASS == "11-20cm" ~ "SIZE_BIN_11_20_CM",
        SIZE_CLASS == "21-30cm" ~ "SIZE_BIN_21_30_CM",
        SIZE_CLASS == "31-40cm" ~ "SIZE_BIN_31_40_CM",
        SIZE_CLASS == "41-50cm" ~ "SIZE_BIN_41_50_CM",
        SIZE_CLASS == "51-60cm" ~ "SIZE_BIN_51_60_CM"
      ))
  
  
  # calculate Biomass
  fish_biomass <- data_formatted %>%
    dplyr::rename(SPECIES = TAXON_CODE) %>%
    #join L-W relationship constants with data to do Biomass Calculations
    left_join(., fish_species_dbase, by = "SPECIES") %>%
    # add column with average size of bin to be length of SIZE_CLASS
    mutate(
      length = case_when(
        SIZE_CLASS == "SIZE_BIN_0_10_CM" ~ "7",
        SIZE_CLASS == "SIZE_BIN_11_20_CM" ~ "15",
        SIZE_CLASS == "SIZE_BIN_21_30_CM" ~ "25",
        SIZE_CLASS == "SIZE_BIN_31_40_CM" ~ "35",
        SIZE_CLASS == "SIZE_BIN_41_50_CM" ~ "45",
        SIZE_CLASS == "SIZE_BIN_51_60_CM" ~ "55",
        TRUE ~ "Other"
      )
    ) %>%
    mutate_at(vars(length), as.numeric) %>%
    # calculate Biomass per fish
    mutate(BIOMASS_PER_FISH_G = LW_A * ((length * LENGTH_CONVERSION_FACTOR) ^
                                          LW_B)) %>%
    # calculate biomass for all fish per row and converted g to kg by /1000
    mutate(BIOMASS_PER_FISH_KG = (COUNT * BIOMASS_PER_FISH_G) / 1000) %>%
    select(REGION:SIZE_CLASS, BIOMASS_PER_FISH_KG)
  
  
  
  fish_bioerosion <- data_formatted %>%
    #combine COUNT data and equations to determine bioerosion rates
    left_join(., rates_dbase_formatted, by = c("TAXON_NAME", "SIZE_CLASS", "PHASE")) %>%
    # clean up
    mutate(FXN_GRP = case_when((TAXON_CODE == "PARR" & SIZE_CLASS != "SIZE_BIN_0_10_CM") ~ "Scraper",
                               SIZE_CLASS == "SIZE_BIN_0_10_CM" ~ "Browser",
                               is.na(FXN_GRP) ~ "Scraper",
                               TRUE ~ FXN_GRP)) %>% # Label PARR as "Scraper"
    mutate(FXN_GRP = case_when((FXN_GRP != "Excavator" & FXN_GRP != "Scraper") ~ "Other",
                               TRUE ~ FXN_GRP)) %>%
    # erosion rates do not include 0-10cm, so need to replace NA with 0
    replace(is.na(.), 0) %>%
    mutate_at(vars(EROSION_RATE), as.numeric) %>%
    # calculate bioerosion value by multiplying COUNT value with bioerosion value
    # happens in summarize_fish_metrics.R script
    mutate(FISH_EROSION_KG_M2_YR = COUNT * EROSION_RATE) %>% 
    select(REGION:TAXON_CODE, FXN_GRP, PHASE, SIZE_CLASS, FISH_EROSION_KG_M2_YR) %>%
    #change all negative bioerosion values to zero...can use this to change
    # multiple columns to zero based on a single column
    mutate_at(.vars = "FISH_EROSION_KG_M2_YR",
              list(~ifelse(FISH_EROSION_KG_M2_YR <= 0, 0, .)))
  
  fish_all <-
    left_join(
      data_formatted,
      fish_biomass,
      by = colnames(data_formatted)[colnames(data_formatted) %in% colnames(fish_biomass)]
    ) %>%
    left_join(
      .,
      fish_bioerosion,
      by = colnames(data_formatted)[colnames(data_formatted) %in% colnames(fish_bioerosion)]
    ) %>%
    mutate(BIOMASS_PER_FISH_KG = replace_na(BIOMASS_PER_FISH_KG, 0)) %>%
    mutate(FISH_EROSION_KG_M2_YR = replace_na(FISH_EROSION_KG_M2_YR, 0))
  
  return(fish_all)
  
  
}

