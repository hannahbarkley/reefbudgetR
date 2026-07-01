#' Calculate urchin erosion rates from urchin census data
#'
#' @author Hannah Barkley 
#'
#' @param data Urchin observation data set.
#' @param full_summary Logical. If TRUE, returns a list of dataframes.
#'
#' @import tidyr
#' @import dplyr
#'
#' @export process_urchins

process_urchins <- function(data,
                            full_summary = TRUE) {
  
  options(dplyr.summarise.inform = FALSE)
  
  # --- 0. STANDARDIZE COLUMN NAMES ---
  
  # Rename DEPTH_M to SITE_DEPTH_M if it exists in the dataset
  if ("DEPTH_M" %in% names(data) && !"SITE_DEPTH_M" %in% names(data)) {
    data <- data %>% rename(SITE_DEPTH_M = DEPTH_M)
  }
  
  # --- 1. SETUP LOOKUP TABLES ---
  
  size_dict <- tibble(
    TEST_SIZE_BIN_MM = c(
      "TEST_SIZE_BIN_0_20_MM", "TEST_SIZE_BIN_21_40_MM", "TEST_SIZE_BIN_41_60_MM",
      "TEST_SIZE_BIN_61_80_MM", "TEST_SIZE_BIN_81_100_MM", "TEST_SIZE_BIN_101_120_MM",
      "TEST_SIZE_BIN_121_140_MM", "TEST_SIZE_BIN_141_160_MM", "TEST_SIZE_BIN_161_180_MM",
      "TEST_SIZE_BIN_181_200_MM"
    ),
    TEST_SIZE_MEDIAN_MM = c(10, 30, 50, 70, 90, 110, 130, 150, 170, 190)
  )
  
  erosion_params <- tibble(
    TAXON_CODE = c("ECMA", "ECOB", "ECSP", "DISP", "ECTH", "ECST", "PAGR", "EUME"),
    PARAM_A = c(0.0003, 0.0003, 0.0003, 0.000003, 0.000003, 0.00004, 0.00004, 0.00004),
    PARAM_B = c(1.9671, 1.9671, 1.9671, 3.2887, 3.2887, 2.6025, 2.6025, 2.6025)
  )
  
  all_urchs <- c("ECST", "ECTH", "DISP", "ECSP", "ECMA", "ECOB", 
                 "PAGR", "EUME", "TRGR", "HEMA", "MEGL")
  
  taxon_names <- c(
    "ECST" = "Echinostrephus sp.", "ECTH" = "Echinothrix sp.", "DISP" = "Diadema sp.",
    "ECSP" = "Echinometra sp.", "ECMA" = "Echinometra mathaei", "ECOB" = "Echinometra oblonga",
    "PAGR" = "Parasalenia gratiosa", "EUME" = "Eucidaris metularia", "TRGR" = "Tripneustes gratilla",
    "HEMA" = "Heterocentrotus mammilatus", "MEGL" = "Mespilia globulus"
  )
  
  # --- 2. PREPARE METADATA & KEYS ---
  
  site_meta <- data %>%
    # explicitly keep transects surveyed for urchins (TRUE or FALSE) but drop NA
    filter(URCH_OBS_TF %in% c(TRUE, FALSE)) %>%
    select(OCC_SITEID, REGION, REGIONCODE, YEAR, CRUISE_ID, LOCATION, 
           LOCATIONCODE, LATITUDE, LONGITUDE, SITE_DEPTH_M, LOCALDATE) %>%
    distinct(OCC_SITEID, .keep_all = TRUE)
  
  transect_meta <- data %>%
    # explicitly keep transects surveyed for urchins (TRUE or FALSE) but drop NA
    filter(URCH_OBS_TF %in% c(TRUE, FALSE)) %>%
    mutate(OCC_SITEID_TRANSECT = paste0(OCC_SITEID, "-", CB_TRANSECTID)) %>%
    distinct(OCC_SITEID, CB_TRANSECTID, OCC_SITEID_TRANSECT, TRANSECT_LENGTH_M)
  
  # --- 3. PROCESS COUNTS & EROSION ---
  
  df_long <- data %>%
    select(OCC_SITEID, CB_TRANSECTID, TAXON_CODE, starts_with("TEST_SIZE_BIN")) %>%
    filter(!is.na(TAXON_CODE)) %>%
    pivot_longer(cols = starts_with("TEST_SIZE_BIN"),
                 names_to = "TEST_SIZE_BIN_MM",
                 values_to = "COUNT") %>%
    left_join(size_dict, by = "TEST_SIZE_BIN_MM") %>%
    left_join(erosion_params, by = "TAXON_CODE") %>%
    mutate(
      COUNT = replace_na(COUNT, 0),
      EROSION_VAL = coalesce(PARAM_A * (TEST_SIZE_MEDIAN_MM ^ PARAM_B), 0)
    )
  
  # --- 4. EXPAND GRID (HANDLE ZEROS SAFELY) ---
  
  # The base grid now explicitly contains those URCH_OBS_TF = FALSE transects
  base_grid <- tidyr::expand_grid(
    transect_meta %>% select(OCC_SITEID, CB_TRANSECTID, TRANSECT_LENGTH_M),
    TAXON_CODE = all_urchs,
    TEST_SIZE_BIN_MM = size_dict$TEST_SIZE_BIN_MM
  )
  
  df_complete <- base_grid %>%
    left_join(df_long, by = c("OCC_SITEID", "CB_TRANSECTID", "TAXON_CODE", "TEST_SIZE_BIN_MM")) %>%
    left_join(size_dict, by = "TEST_SIZE_BIN_MM", suffix = c("", "_dict")) %>% 
    mutate(
      # Replace NA counts/erosion with 0 for empty transects that just got rejoined
      COUNT = replace_na(COUNT, 0),
      EROSION_VAL = coalesce(EROSION_VAL, 0),
      TEST_SIZE_MEDIAN_MM = coalesce(TEST_SIZE_MEDIAN_MM, TEST_SIZE_MEDIAN_MM_dict),
      TAXON_NAME = recode_factor(TAXON_CODE, !!!taxon_names),
      TRANSECT_URCHIN_DENSITY_NO_M2 = COUNT / TRANSECT_LENGTH_M,
      TRANSECT_EROSION_G_M2_YR = EROSION_VAL * TRANSECT_URCHIN_DENSITY_NO_M2 * 365
    )
  
  # --- 5. SUMMARIZE DATA ---
  
  join_meta <- function(df) {
    left_join(df, site_meta, by = "OCC_SITEID") %>%
      select(REGION, REGIONCODE, YEAR, CRUISE_ID, LOCATION, LOCATIONCODE, OCC_SITEID,
             LATITUDE, LONGITUDE, SITE_DEPTH_M, LOCALDATE, everything())
  }
  
  transect_density_taxon <- df_complete %>%
    group_by(OCC_SITEID, CB_TRANSECTID, TAXON_CODE, TAXON_NAME, TEST_SIZE_BIN_MM, TEST_SIZE_MEDIAN_MM, TRANSECT_LENGTH_M) %>%
    summarize(
      TRANSECT_URCHIN_ABUNDANCE_NO = sum(COUNT),
      TRANSECT_URCHIN_DENSITY_NO_M2 = sum(COUNT) / mean(TRANSECT_LENGTH_M),
      TRANSECT_EROSION_G_M2_YR = sum(TRANSECT_EROSION_G_M2_YR),
      .groups = "drop"
    ) 
  
  transect_erosion <- transect_density_taxon %>%
    group_by(OCC_SITEID, CB_TRANSECTID, TRANSECT_LENGTH_M) %>%
    summarize(
      URCHIN_ABUNDANCE_NO = sum(TRANSECT_URCHIN_ABUNDANCE_NO),
      URCHIN_EROSION_KG_M2_YR = sum(TRANSECT_EROSION_G_M2_YR) / 1000,
      .groups = "drop"
    ) %>%
    mutate(
      OCC_SITEID_TRANSECT = paste(OCC_SITEID, CB_TRANSECTID, sep = "-"),
      URCHIN_DENSITY_NO_M2 = URCHIN_ABUNDANCE_NO / TRANSECT_LENGTH_M
    ) %>%
    join_meta()
  
  transect_erosion_taxon <- transect_density_taxon %>%
    group_by(OCC_SITEID, CB_TRANSECTID, TAXON_CODE, TAXON_NAME, TRANSECT_LENGTH_M) %>%
    summarize(
      URCHIN_ABUNDANCE_NO = sum(TRANSECT_URCHIN_ABUNDANCE_NO),
      URCHIN_EROSION_KG_M2_YR = sum(TRANSECT_EROSION_G_M2_YR) / 1000,
      .groups = "drop"
    ) %>%
    mutate(URCHIN_DENSITY_NO_M2 = URCHIN_ABUNDANCE_NO / TRANSECT_LENGTH_M) %>%
    join_meta()
  
  transect_erosion_taxon_size <- transect_density_taxon %>%
    mutate(URCHIN_EROSION_KG_M2_YR = TRANSECT_EROSION_G_M2_YR / 1000) %>%
    join_meta() %>%
    select(-TRANSECT_EROSION_G_M2_YR) 
  
  calc_stats <- function(df, group_vars) {
    df %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(
        TRANSECT_COUNT = n(), 
        URCHIN_EROSION_KG_M2_YR_MEAN = mean(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE),
        URCHIN_EROSION_KG_M2_YR_SD = sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE),
        URCHIN_EROSION_KG_M2_YR_SE = sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE) / sqrt(n()),
        URCHIN_DENSITY_NO_M2_MEAN = mean(URCHIN_DENSITY_NO_M2, na.rm = TRUE),
        URCHIN_DENSITY_NO_M2_SD = sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE),
        URCHIN_DENSITY_NO_M2_SE = sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE) / sqrt(n()),
        URCHIN_ABUNDANCE_NO = sum(URCHIN_ABUNDANCE_NO),
        .groups = "drop"
      ) %>%
      join_meta()
  }
  
  site_erosion <- calc_stats(transect_erosion, "OCC_SITEID")
  site_erosion_taxon <- calc_stats(transect_erosion_taxon, c("OCC_SITEID", "TAXON_CODE", "TAXON_NAME"))
  
  site_erosion_taxon_size <- transect_erosion_taxon_size %>%
    group_by(OCC_SITEID, TAXON_CODE, TAXON_NAME, TEST_SIZE_BIN_MM) %>%
    summarize(
      TRANSECT_COUNT = n(), 
      URCHIN_ABUNDANCE_NO = sum(TRANSECT_URCHIN_ABUNDANCE_NO),
      URCHIN_EROSION_KG_M2_YR_MEAN = mean(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE),
      URCHIN_EROSION_KG_M2_YR_SD = sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE),
      URCHIN_EROSION_KG_M2_YR_SE = sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    join_meta()
  
  # --- 6. OUTPUT ---
  
  if (full_summary) {
    return(list(
      site_erosion = site_erosion,
      transect_erosion = transect_erosion,
      site_taxon = site_erosion_taxon,
      transect_taxon = transect_erosion_taxon,
      site_taxon_size = site_erosion_taxon_size,
      transect_taxon_size = transect_erosion_taxon_size,
      data = df_complete %>% join_meta() 
    ))
  } else {
    return(site_erosion)
  }
}