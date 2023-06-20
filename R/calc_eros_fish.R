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
#'@examples calc_eros_fish_output <- calc_eros_fish(fish_data, rates_dbase)


calc_eros_fish <- function(data,
                           rates_dbase) {
  
  
  # Format dataframe for biomass and bioerosion calculations below
  data_formatted <- data %>%
    select(CRUISE_ID:PHASE) %>%
    #change characters to factors
    dplyr::mutate_at(vars(
      c(
        REGION:VISIBILITY_M,
        TRANSECT,
        HABITAT_TYPE:SPECIES,
        PHASE
      )
    ), as.factor) %>%
    # convert integers to numbers
    mutate_at(vars(c(
      TRANSECT_LENGTH_M:AREA_M2,
      "0-10cm":"51-60cm"
    )), as.numeric) %>%
    # select only the important columns
    select(CRUISE_ID:TRANSECT,
           AREA_M2,
           SPECIES:PHASE) %>%
    #column names to row values
    gather(.,
           "SIZE_CLASS",
           "COUNT", -c(CRUISE_ID:SPECIES),
           -PHASE, na.rm = TRUE) %>%
    # make size bins into factors
    mutate_at(vars(SIZE_CLASS), as.factor) %>%
    #remove all size classes above 10cm for J classification
    filter(!(
      SIZE_CLASS %in% c("11-20cm",
                        "21-30cm",
                        "31-40cm",
                        "41-50cm",
                        "51-60cm") &
        PHASE == "J"
    )) %>%
    # remove all I and T for 0-10cm classification
    filter(!(SIZE_CLASS %in% "0-10cm" &
               PHASE %in% c("I", "T"))) %>%
    # remove all sizes above 51cm for initial phases
    filter(!(SIZE_CLASS %in% "51-60cm" &
               PHASE == "I"))
  
  # calculate Biomass
  fish_biomass <- data_formatted %>%
    #join L-W relationship constants with data to do Biomass Calculations
    left_join(., fish_species_dbase, by = "SPECIES") %>%
    # add column with average size of bin to be length of SIZE_CLASS
    mutate(
      length = case_when(
        SIZE_CLASS == "0-10cm" ~ "7",
        SIZE_CLASS == "11-20cm" ~ "15",
        SIZE_CLASS == "21-30cm" ~ "25",
        SIZE_CLASS == "31-40cm" ~ "35",
        SIZE_CLASS == "41-50cm" ~ "45",
        SIZE_CLASS == "51-60cm" ~ "55",
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
    left_join(., fish_functional_groups %>%
                #join spdata b/c need full scientific name
                select(SPECIES, TAXONNAME), by = "SPECIES") %>% # remove unwanted columns from the join
    dplyr::rename(TAXON_NAME = TAXONNAME) %>%
    #combine COUNT data and equations to determine bioerosion rates
    left_join(., rates_dbase, by = c("TAXON_NAME", "SIZE_CLASS", "PHASE")) %>%
    # clean up
    mutate(FXN_GRP = case_when((SPECIES == "PARR" & SIZE_CLASS != "0-10cm") ~ "Scraper",
                               SIZE_CLASS == "0-10cm" ~ "Browser",
                               is.na(FXN_GRP) ~ "Scraper",
                               TRUE ~ FXN_GRP)) %>% # Label PARR as "Scraper"
    mutate(FXN_GRP = case_when((FXN_GRP != "Excavator" & FXN_GRP != "Scraper") ~ "Other",
                               TRUE ~ FXN_GRP)) %>%
    mutate(TAXON_NAME = case_when(SPECIES == "PARR" ~ "Parrotfish",
                                  TRUE ~ TAXON_NAME)) %>% # Label Taxonname of PARR with Parrotfish
    # erosion rates do not include 0-10cm, so need to replace NA with 0
    replace(is.na(.), 0) %>%
    mutate_at(vars(EROSION_RATE), as.numeric) %>%
    # calculate bioerosion value by multiplying COUNT value with bioerosion value
    # happens in summarize_fish_metrics.R script
    mutate(FISH_EROSION_KG_M2_YR = COUNT * EROSION_RATE) %>% 
    select(REGION:SPECIES, FXN_GRP, PHASE, SIZE_CLASS, FISH_EROSION_KG_M2_YR) %>%
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

