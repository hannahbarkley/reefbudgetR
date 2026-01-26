#' Calculate urchin erosion rates from urchin census data
#'
#'@author Hannah Barkley
#'
#'@param data Urchin observation data set.
#'
#'@import tidyr
#'@import dplyr
#'@importFrom sjmisc seq_row
#'
#'@export process_urchins
#'
#'@examples
#' urchin_data <- read.csv("ESD_CarbBudget_Urchins_OAHU_2021.csv",
#'   na = "", check.names = FALSE)
#'
#' urch<- process_urchins(
#'   data = urchin_data
#' )

process_urchins <- function(data,
                            transect_length = NULL,
                            full_summary = TRUE) {
  
  options(dplyr.summarise.inform = FALSE)
  
  # Create transect identifier
  data$OCC_SITEID_TRANSECT <- paste0(data$OCC_SITEID, "-", data$CB_TRANSECTID)
  
  # Extract site metadata
  site_meta <- data %>%
    distinct(OCC_SITEID, .keep_all = TRUE) %>%
    select(OCC_SITEID, REGION, REGIONCODE, YEAR, CRUISE_ID, LOCATION, 
           LOCATIONCODE, LATITUDE, LONGITUDE, SITE_DEPTH_M, LOCALDATE)
  
  # Create long data frame of urchin counts
  data_counts <- data %>%
    select(OCC_SITEID, CB_TRANSECTID, TAXON_CODE, TAXON_NAME, starts_with("TEST_SIZE_BIN")) %>%
    filter(!is.na(TAXON_CODE)) %>%
    pivot_longer(
      cols = starts_with("TEST_SIZE_BIN"),
      names_to = "TEST_SIZE_BIN_MM",
      values_to = "COUNT"
    )
  
  # Set up full reference grid with all urchins and size classes
  all_urchs <- c("ECST", "ECTH", "DISP", "ECSP", "ECMA", "ECOB", 
                 "PAGR", "EUME", "TRGR", "HEMA", "MEGL")
  
  all_bins <- c("TEST_SIZE_BIN_0_20_MM", "TEST_SIZE_BIN_21_40_MM", "TEST_SIZE_BIN_41_60_MM",
                "TEST_SIZE_BIN_61_80_MM", "TEST_SIZE_BIN_81_100_MM", "TEST_SIZE_BIN_101_120_MM",
                "TEST_SIZE_BIN_121_140_MM", "TEST_SIZE_BIN_141_160_MM", "TEST_SIZE_BIN_161_180_MM",
                "TEST_SIZE_BIN_181_200_MM")
  
  full_grid <- expand.grid(
    OCC_SITEID = unique(data$OCC_SITEID),
    CB_TRANSECTID = unique(data$CB_TRANSECTID),
    TAXON_CODE = all_urchs,
    TEST_SIZE_BIN_MM = all_bins,
    stringsAsFactors = FALSE
  )
  
  # Populate the full urchin data with observed counts
  data_full <- full_grid %>%
    # attach metadata
    left_join(site_meta, by = "OCC_SITEID") %>%
    # attach counts
    left_join(data_counts, by = c("OCC_SITEID", "CB_TRANSECTID", "TAXON_CODE", "TEST_SIZE_BIN_MM")) %>%
    mutate(
      COUNT = replace_na(COUNT, 0),
      TAXON_NAME = dplyr::recode_factor(coalesce(TAXON_NAME, TAXON_CODE), # Handle missing names
                                        "ECST" = "Echinostrephus sp.", "ECTH" = "Echinothrix sp.", "DISP" = "Diadema sp.",
                                        "ECSP" = "Echinometra sp.", "ECMA" = "Echinometra mathaei", "ECOB" = "Echinometra oblonga",
                                        "PAGR" = "Parasalenia gratiosa", "EUME" = "Eucidaris metularia", "TRGR" = "Tripneustes gratilla",
                                        "HEMA" = "Heterocentrotus mammilatus", "MEGL" = "Mespilia globulus"
      ),
      OCC_SITEID_TRANSECT = paste0(OCC_SITEID, "-", CB_TRANSECTID)
    ) %>%
    # attach transect length
    left_join(
      distinct(data, OCC_SITEID_TRANSECT, TRANSECT_LENGTH_M), 
      by = "OCC_SITEID_TRANSECT"
    )
  
  # calculate erosion rates
  data_eroders <- data_full %>%
    filter(!TAXON_CODE %in% c("TRGR", "HEMA", "MEGL")) %>%
    mutate(
      TEST_SIZE_BIN_MM = factor(TEST_SIZE_BIN_MM, levels = all_bins),

      TEST_SIZE_MEDIAN_MM = case_when(
        TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_0_20_MM"   ~ 10,
        TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_21_40_MM"  ~ 30,
        TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_41_60_MM"  ~ 50,
        TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_61_80_MM"  ~ 70,
        TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_81_100_MM" ~ 90,
        TRUE ~ 110 
      ),
      
      TRANSECT_URCHIN_DENSITY_NO_M2 = COUNT / TRANSECT_LENGTH_M,
      
      # Single vectorized case_when for erosion math
      TRANSECT_EROSION_G_M2_YR = case_when(
        TAXON_CODE %in% c("ECMA", "ECOB", "ECSP") ~ 
          0.0003 * (TEST_SIZE_MEDIAN_MM ^ 1.9671) * TRANSECT_URCHIN_DENSITY_NO_M2 * 365,
        
        TAXON_CODE %in% c("DISP", "ECTH") ~ 
          0.000003 * (TEST_SIZE_MEDIAN_MM ^ 3.2887) * TRANSECT_URCHIN_DENSITY_NO_M2 * 365,
        
        TAXON_CODE %in% c("ECST", "PAGR", "EUME") ~ 
          0.00004 * (TEST_SIZE_MEDIAN_MM ^ 2.6025) * TRANSECT_URCHIN_DENSITY_NO_M2 * 365,
        
        TRUE ~ 0
      )
    )
  
  # Summarize data
  
  # Transect Density by Taxon
  transect_density_taxon <- data_eroders %>%
    group_by(REGION, REGIONCODE, YEAR, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID,
             LATITUDE, LONGITUDE, SITE_DEPTH_M, LOCALDATE, CB_TRANSECTID,
             TAXON_CODE, TAXON_NAME, TEST_SIZE_BIN_MM, TRANSECT_LENGTH_M, 
             TEST_SIZE_MEDIAN_MM) %>%
    summarize(
      TRANSECT_URCHIN_ABUNDANCE_NO = sum(COUNT),
      TRANSECT_URCHIN_DENSITY_NO_M2 = sum(COUNT) / mean(TRANSECT_LENGTH_M, na.rm = TRUE),
      TRANSECT_EROSION_G_M2_YR = sum(TRANSECT_EROSION_G_M2_YR),
      .groups = "drop"
    )
  
  # Site Density by Taxon
  site_density_taxon <- transect_density_taxon %>%
    group_by(REGION, REGIONCODE, YEAR, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID,
             LATITUDE, LONGITUDE, SITE_DEPTH_M, LOCALDATE, TAXON_CODE, TAXON_NAME,
             TEST_SIZE_BIN_MM, TRANSECT_LENGTH_M) %>%
    summarize(
      URCHIN_ABUNDANCE_NO = sum(TRANSECT_URCHIN_ABUNDANCE_NO),
      URCHIN_DENSITY_NO_M2_MEAN = mean(TRANSECT_URCHIN_DENSITY_NO_M2),
      URCHIN_DENSITY_NO_M2_SD = sd(TRANSECT_URCHIN_DENSITY_NO_M2),
      URCHIN_DENSITY_NO_M2_SE = sd(TRANSECT_URCHIN_DENSITY_NO_M2) / sqrt(length(TRANSECT_URCHIN_DENSITY_NO_M2)),
      URCHIN_DENSITY_NO_M2_N = length(TRANSECT_URCHIN_DENSITY_NO_M2),
      .groups = "drop"
    )
  
  # Transect Erosion Total
  transect_erosion <- transect_density_taxon %>%
    group_by(REGION, REGIONCODE, YEAR, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID,
             LATITUDE, LONGITUDE, SITE_DEPTH_M, LOCALDATE, CB_TRANSECTID, TRANSECT_LENGTH_M) %>%
    summarize(
      URCHIN_ABUNDANCE_NO = sum(TRANSECT_URCHIN_ABUNDANCE_NO),
      URCHIN_EROSION_KG_M2_YR = sum(TRANSECT_EROSION_G_M2_YR) / 1000,
      .groups = "drop"
    ) %>%
    mutate(
      OCC_SITEID_TRANSECT = paste(OCC_SITEID, CB_TRANSECTID, sep = "-"),
      URCHIN_DENSITY_NO_M2 = URCHIN_ABUNDANCE_NO / TRANSECT_LENGTH_M
    ) %>%
    select(REGION, REGIONCODE, YEAR, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID,
           LATITUDE, LONGITUDE, SITE_DEPTH_M, LOCALDATE, CB_TRANSECTID, 
           OCC_SITEID_TRANSECT, TRANSECT_LENGTH_M, URCHIN_ABUNDANCE_NO, 
           URCHIN_DENSITY_NO_M2, URCHIN_EROSION_KG_M2_YR)
  
  # Transect Erosion by Taxon
  transect_erosion_taxon <- transect_density_taxon %>%
    group_by(REGION, REGIONCODE, YEAR, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID,
             LATITUDE, LONGITUDE, SITE_DEPTH_M, LOCALDATE, CB_TRANSECTID, TAXON_CODE,
             TAXON_NAME, TRANSECT_LENGTH_M) %>%
    summarize(
      URCHIN_ABUNDANCE_NO = sum(TRANSECT_URCHIN_ABUNDANCE_NO),
      URCHIN_EROSION_KG_M2_YR = sum(TRANSECT_EROSION_G_M2_YR) / 1000,
      .groups = "drop"
    ) %>%
    mutate(URCHIN_DENSITY_NO_M2 = URCHIN_ABUNDANCE_NO / TRANSECT_LENGTH_M)
  
  # Site Erosion Total
  site_erosion <- transect_erosion %>%
    group_by(REGION, REGIONCODE, YEAR, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID,
             LATITUDE, LONGITUDE, SITE_DEPTH_M, LOCALDATE) %>%
    summarize(
      URCHIN_EROSION_KG_M2_YR_MEAN = mean(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE),
      URCHIN_EROSION_KG_M2_YR_SD = sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE),
      URCHIN_EROSION_KG_M2_YR_N = length(URCHIN_EROSION_KG_M2_YR),
      URCHIN_EROSION_KG_M2_YR_SE = sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE) / length(URCHIN_EROSION_KG_M2_YR),
      URCHIN_EROSION_KG_M2_YR_CI = 1.97 * sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE) / length(URCHIN_EROSION_KG_M2_YR),
      URCHIN_DENSITY_NO_M2_MEAN = mean(URCHIN_DENSITY_NO_M2, na.rm = TRUE),
      URCHIN_DENSITY_NO_M2_SD = sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE),
      URCHIN_DENSITY_NO_M2_N = length(URCHIN_DENSITY_NO_M2),
      URCHIN_DENSITY_NO_M2_SE = sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE) / length(URCHIN_DENSITY_NO_M2),
      URCHIN_DENSITY_NO_M2_CI = 1.97 * sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE) / length(URCHIN_DENSITY_NO_M2),
      URCHIN_ABUNDANCE_NO = sum(URCHIN_ABUNDANCE_NO),
      .groups = "drop"
    )
  
  # Site Erosion by Taxon
  site_erosion_taxon <- transect_erosion_taxon %>%
    group_by(REGION, REGIONCODE, YEAR, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID,
             LATITUDE, LONGITUDE, SITE_DEPTH_M, LOCALDATE, TAXON_CODE, TAXON_NAME) %>%
    summarize(
      URCHIN_EROSION_KG_M2_YR_MEAN = mean(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE),
      URCHIN_EROSION_KG_M2_YR_SD = sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE),
      URCHIN_EROSION_KG_M2_YR_N = length(URCHIN_EROSION_KG_M2_YR),
      URCHIN_EROSION_KG_M2_YR_SE = sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE) / length(URCHIN_EROSION_KG_M2_YR),
      URCHIN_EROSION_KG_M2_YR_CI = 1.97 * sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE) / length(URCHIN_EROSION_KG_M2_YR),
      URCHIN_DENSITY_NO_M2_MEAN = mean(URCHIN_DENSITY_NO_M2, na.rm = TRUE),
      URCHIN_DENSITY_NO_M2_SD = sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE),
      URCHIN_DENSITY_NO_M2_N = length(URCHIN_DENSITY_NO_M2),
      URCHIN_DENSITY_NO_M2_SE = sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE) / length(URCHIN_DENSITY_NO_M2),
      URCHIN_DENSITY_NO_M2_CI = 1.97 * sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE) / length(URCHIN_DENSITY_NO_M2),
      URCHIN_ABUNDANCE_NO = sum(URCHIN_ABUNDANCE_NO),
      .groups = "drop"
    )
  
  # 8. Filter Output (Strictly keeping requested columns)
  data_out <- data[c(
    "REGION", "REGIONCODE", "YEAR", "CRUISE_ID", "LOCATION", "LOCATIONCODE",
    "OCC_SITEID", "LATITUDE", "LONGITUDE", "SITE_DEPTH_M", "LOCALDATE",
    "CB_TRANSECTID", "TRANSECT_LENGTH_M", "URCH_OBS_TF", "TAXON_NAME",
    "TAXON_CODE", "TEST_SIZE_BIN_0_20_MM", "TEST_SIZE_BIN_21_40_MM",
    "TEST_SIZE_BIN_41_60_MM", "TEST_SIZE_BIN_61_80_MM", "TEST_SIZE_BIN_81_100_MM",
    "TEST_SIZE_BIN_101_120_MM", "TEST_SIZE_BIN_121_140_MM", "TEST_SIZE_BIN_141_160_MM",
    "TEST_SIZE_BIN_161_180_MM", "TEST_SIZE_BIN_181_200_MM"
  )]
  
  # 9. Return
  if (full_summary == TRUE) {
    return(list(
      site_erosion = site_erosion,
      transect_erosion = transect_erosion,
      site_taxon = site_erosion_taxon,
      transect_taxon = transect_erosion_taxon,
      data = data_out
    ))
  }
  if (full_summary == FALSE) {
    return(site_erosion = site_erosion)
  }
}
