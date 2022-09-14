#' Calculate parrotfish biomass, density, and bioerosion rates

#'@author Rebecca Weible
#'
#'@param data Parrotfish belt data, including number of fish observed of each
#'species, size class, and phase.
#'@param rates_dbase Erosion rates database to use. Choose either Indo-Pacific
#'ReefBudget ("rates_dbase = "IPRB") or U.S. Pacific Islands rates developed
#'by Tye Kindinger, NOAA PIFSC ("rates_dbase = "Kindinger").
#'
#'@import dplyr
#'@importFrom rlang .data
#'
#'@export calc_eros_fish

calc_eros_fish <- function(data,
                             rates_dbase = c("IPRB", "Kindinger")) {
  if (rates_dbase == "IPRB") {
    rates_dbase <- fish_erosion_dbase_iprb
  }

  if (rates_dbase == "Kindinger") {
    rates_dbase <- fish_erosion_dbase_kindinger
  }

  # Format dataframe for biomass and bioerosion calculations below
  data_formatted <- data %>%
    select(.data$REGION:.data$PHASE) %>%
    #change characters to factors
    dplyr::mutate_at(vars(
      c(
        .data$REGION:.data$VISIBILITY_M,
        .data$CB_TRANSECTID,
        .data$HABITAT_TYPE:.data$SPECIES,
        .data$PHASE
      )
    ), as.factor) %>%
    # convert integers to numbers
    mutate_at(vars(c(
      .data$TRANSECT_LENGTH_M:.data$AREA_M2,
      "SIZE_CLASS_0_10_CM":"SIZE_CLASS_51_60_CM"
    )), as.numeric) %>%
    # select only the important columns
    select(.data$REGION:.data$CB_TRANSECTID,
           .data$AREA_M2,
           .data$SPECIES:PHASE) %>%
    #column names to row values
    gather(.,
           "SIZE_CLASS",
           "COUNT", -c(.data$REGION:.data$SPECIES),
           -.data$PHASE, na.rm = TRUE) %>%
    # make size bins into factors
    mutate_at(vars(SIZE_CLASS), as.factor) %>%
    #remove all size classes above 10cm for J classification
    filter(!(
      SIZE_CLASS %in% c("SIZE_CLASS_11_20_CM",
                        "SIZE_CLASS_21_30_CM",
                        "SIZE_CLASS_31_40_CM",
                        "SIZE_CLASS_41_50_CM",
                        "SIZE_CLASS_51_60_CM") &
        PHASE == "J"
    )) %>%
    # remove all I and T for 0-10cm classification
    filter(!(SIZE_CLASS %in% "SIZE_CLASS_0_10_CM" &
               PHASE %in% c("I", "T"))) %>%
    # remove all sizes above 51cm for initial phases
    filter(!(SIZE_CLASS %in% "SIZE_CLASS_51_60_CM" &
               PHASE == "I"))

  data_formatted$SIZE_CLASS <- recode_factor(
    data_formatted$SIZE_CLASS,
    "SIZE_CLASS_0_10_CM" = "0-10cm",
    "SIZE_CLASS_11_20_CM" = "11-20cm",
    "SIZE_CLASS_21_30_CM" = "21-30cm",
    "SIZE_CLASS_31_40_CM" = "31-40cm",
    "SIZE_CLASS_41_50_CM" = "41-50cm",
    "SIZE_CLASS_51_60_CM" = "51-60cm"
  )

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
    mutate(BIOMASS_KG_HA = (COUNT * BIOMASS_PER_FISH_G) / 1000) %>%
    select(REGION:SIZE_CLASS, BIOMASS_KG_HA)


  fish_bioerosion <- data_formatted %>%
    left_join(., fish_grazing_types %>%
                #join spdata b/c need full scientific name
                select(SPECIES, FISH_sciname), by = "SPECIES") %>%
    # remove unwanted columns from the join
    dplyr::rename(TAXON_NAME = FISH_sciname) %>%
    #combine COUNT data and equations to determine bioerosion rates
    left_join(., rates_dbase, by = c("TAXON_NAME", "SIZE_CLASS", "PHASE")) %>%
    # erosion rates do not include 0-10cm, so need to replace NA with 0
    replace(is.na(.), 0) %>%
    mutate_at(vars(EROSION_RATE), as.numeric) %>%
    # calculate bioerosion value by multiplying COUNT value with
    # bioerosion rate value, dividing by 10000 to cancel out the /10000 that
    # happens in the next function
    mutate(FISH_EROSION_KG_M2_YR = (COUNT * EROSION_RATE) / 10000) %>%
    select(REGION:SPECIES, PHASE, SIZE_CLASS, FISH_EROSION_KG_M2_YR) %>%
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
    mutate(BIOMASS_KG_HA = replace_na(BIOMASS_KG_HA, 0)) %>%
    mutate(FISH_EROSION_KG_M2_YR = replace_na(FISH_EROSION_KG_M2_YR, 0))

  return(fish_all)

}
