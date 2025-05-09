#' Calculate urchin erosion rates from urchin census data
#'
#'@author Hannah Barkley
#'
#'@param data Urchin observation data set.
#'@param method_name Transect design by which data were collected ("IPRB", "Chords", or "SfM").
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
#' urch_iprb <- process_urchins(
#'   data = urchin_data[urchin_data$CB_METHOD == "IPRB", ],
#'   method_name = "IPRB"
#' )
#'
#' urch_chords <- process_urchins(
#'   data = urchin_data[urchin_data$CB_METHOD == "Chords", ],
#'   method_name = "Chords"
#' )

process_urchins <- function(data,
                            transect_length = NULL,
                            method_name = NULL,
                            full_summary = TRUE) {
  options(dplyr.summarise.inform = FALSE)

  # Create data frame with transect info
  data$OCC_SITEID_TRANSECT <-
    paste0(data$OCC_SITEID, "-", data$CB_TRANSECTID)

  summary_transect <-
    unique(data[c("OCC_SITEID_TRANSECT", "TRANSECT_LENGTH_M")])

  # Convert to long format
  data_long <- data %>%
    dplyr::select(-c(URCH_OBS_TF)) %>%
    filter(is.na(TAXON_CODE) == FALSE) %>%
    pivot_longer(
      cols = !c(
        "SITEVISITID",
        "TRANSECTID",
        "CBURCHINID",
        "REGION",
        "REGIONCODE",
        "YEAR",
        "CRUISE_ID",
        "LOCATION",
        "LOCATIONCODE",
        "OCC_SITEID",
        "LATITUDE",
        "LONGITUDE",
        "DEPTH_M",
        "LOCALDATE",
        "CB_TRANSECTID",
        "OCC_SITEID_TRANSECT",
        "TRANSECT_LENGTH_M",
        "TAXON_NAME",
        "TAXON_CODE"
      ),
      names_to = "TEST_SIZE_BIN_MM",
      values_to = "COUNT"
    )

  # Define full urchin table
  all_urchs <- c("ECST",
                 "ECTH",
                 "DISP",
                 "ECSP",
                 "ECMA",
                 "ECOB",
                 "PAGR",
                 "EUME",
                 "TRGR",
                 "HEMA",
                 "MEGL")
  all_bins <- c(
    "TEST_SIZE_BIN_0_20_MM",
    "TEST_SIZE_BIN_21_40_MM",
    "TEST_SIZE_BIN_41_60_MM",
    "TEST_SIZE_BIN_61_80_MM",
    "TEST_SIZE_BIN_81_100_MM",
    "TEST_SIZE_BIN_101_120_MM",
    "TEST_SIZE_BIN_121_140_MM",
    "TEST_SIZE_BIN_141_160_MM",
    "TEST_SIZE_BIN_161_180_MM",
    "TEST_SIZE_BIN_181_200_MM"
  )

  full_table <- expand.grid(
    OCC_SITEID = unique(data$OCC_SITEID),
    CB_TRANSECTID = unique(data$CB_TRANSECTID),
    TAXON_CODE = all_urchs,
    TEST_SIZE_BIN_MM = all_bins
  )


  # Populate full table
  data_full <- full_join(
    full_table,
    data_long,
    by = c(
      "OCC_SITEID",
      "CB_TRANSECTID",
      "TAXON_CODE",
      "TEST_SIZE_BIN_MM"
    )
  )

  data_full <-
    data_full %>% group_by(OCC_SITEID) %>%
    fill(
      c(
        "REGION",
        "REGIONCODE",
        "YEAR",
        "CRUISE_ID",
        "LOCATION",
        "LOCATIONCODE",
        "OCC_SITEID",
        "LATITUDE",
        "LONGITUDE",
        "DEPTH_M",
        "TRANSECT_LENGTH_M"
      ),
      .direction = "downup"
    ) %>% group_by(OCC_SITEID) %>%
    fill(c("LOCALDATE"),
         .direction = "downup")

  for (i in 1:nrow(data_full)){
    if (is.na(data_full$LOCATIONCODE[i]) == TRUE){
      data_full$REGION[i] <- data$REGION[match(data_full$OCC_SITEID[i], data$OCC_SITEID)]
      data_full$REGIONCODE[i] <- data$REGIONCODE[match(data_full$OCC_SITEID[i], data$OCC_SITEID)]
      data_full$YEAR[i] <- data$YEAR[match(data_full$OCC_SITEID[i], data$OCC_SITEID)]
      data_full$CRUISE_ID[i] <- data$CRUISE_ID[match(data_full$OCC_SITEID[i], data$OCC_SITEID)]
      data_full$LOCATION[i] <- data$LOCATION[match(data_full$OCC_SITEID[i], data$OCC_SITEID)]
      data_full$LOCATIONCODE[i] <- data$LOCATIONCODE[match(data_full$OCC_SITEID[i], data$OCC_SITEID)]
      data_full$LATITUDE[i] <- data$LATITUDE[match(data_full$OCC_SITEID[i], data$OCC_SITEID)]
      data_full$LONGITUDE[i] <- data$LONGITUDE[match(data_full$OCC_SITEID[i], data$OCC_SITEID)]
      data_full$DEPTH_M[i] <- data$DEPTH_M[match(data_full$OCC_SITEID[i], data$OCC_SITEID)]
      data_full$TRANSECT_LENGTH_M[i] <- data$TRANSECT_LENGTH_M[match(data_full$OCC_SITEID[i], data$OCC_SITEID)]
      data_full$LOCALDATE[i] <- data$LOCALDATE[match(data_full$OCC_SITEID[i], data$OCC_SITEID)]


    }
  }

  data_full$TAXON_NAME[is.na(data_full$TAXON_NAME)] <-
    data_full$TAXON_CODE[is.na(data_full$TAXON_NAME)]

  data_full$TAXON_NAME <- recode_factor(
    data_full$TAXON_NAME,
    "ECST" = "Echinostrephus sp.",
    "ECTH" = "Echinothrix sp.",
    "DISP" = "Diadema sp.",
    "ECSP" = "Echinometra sp.",
    "ECMA" = "Echinometra mathaei",
    "ECOB" = "Echinometra oblonga",
    "PAGR" = "Parasalenia gratiosa",
    "EUME" = "Eucidaris metularia",
    "TRGR" = "Tripneustes gratilla",
    "HEMA" = "Heterocentrotus mammilatus",
    "MEGL" =  "Mespilia globulus"
  )

  # Set all NA values to 0
  data_full$COUNT[is.na(data_full$COUNT)] <- 0

  # Add transect length to data frame

  data_full$OCC_SITEID_TRANSECT <-
    paste0(data_full$OCC_SITEID, "-", data_full$CB_TRANSECTID)

  data_full$TRANSECT_LENGTH_M <-
    summary_transect$TRANSECT_LENGTH_M[match(data_full$OCC_SITEID_TRANSECT , summary_transect$OCC_SITEID_TRANSECT)]

  # Remove non-eroding urchin taxa
  data_eroders <- data_full[!data_full$TAXON_CODE %in% c("TRGR", "HEMA", "MEGL"),]

  # Set factor levels and order
  data_eroders$TEST_SIZE_BIN_MM <- factor(
    data_eroders$TEST_SIZE_BIN_MM,
    levels = c(
      "TEST_SIZE_BIN_0_20_MM",
      "TEST_SIZE_BIN_21_40_MM",
      "TEST_SIZE_BIN_41_60_MM",
      "TEST_SIZE_BIN_61_80_MM",
      "TEST_SIZE_BIN_81_100_MM",
      "TEST_SIZE_BIN_101_120_MM",
      "TEST_SIZE_BIN_121_140_MM",
      "TEST_SIZE_BIN_141_160_MM",
      "TEST_SIZE_BIN_161_180_MM",
      "TEST_SIZE_BIN_181_200_MM"
    )
  )

  # Calculate transect-level density and abundance
  transect_density_taxon <- data_eroders %>%
    group_by(
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      LATITUDE,
      LONGITUDE,
      DEPTH_M,
      LOCALDATE,
      CB_TRANSECTID,
      TAXON_CODE,
      TAXON_NAME,
      TEST_SIZE_BIN_MM,
      TRANSECT_LENGTH_M
    ) %>%
    summarize(
      TRANSECT_URCHIN_ABUNDANCE_NO = sum(COUNT),
      TRANSECT_URCHIN_DENSITY_NO_M2 = sum(COUNT) / mean(TRANSECT_LENGTH_M, na.rm = TRUE)
    )

  # Calculate site level density
  site_density_taxon <- transect_density_taxon %>%
    group_by(
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      LATITUDE,
      LONGITUDE,
      DEPTH_M,
      LOCALDATE,
      TAXON_CODE,
      TAXON_NAME,
      TEST_SIZE_BIN_MM,
      TRANSECT_LENGTH_M
    ) %>%
    summarize(
      URCHIN_ABUNDANCE_NO = sum(TRANSECT_URCHIN_ABUNDANCE_NO),
      URCHIN_DENSITY_NO_M2_MEAN = mean(TRANSECT_URCHIN_DENSITY_NO_M2),
      URCHIN_DENSITY_NO_M2_SD = sd(TRANSECT_URCHIN_DENSITY_NO_M2),
      URCHIN_DENSITY_NO_M2_SE = sd(TRANSECT_URCHIN_DENSITY_NO_M2)/
        sqrt(length(sd(TRANSECT_URCHIN_DENSITY_NO_M2))),
     URCHIN_DENSITY_NO_M2_N = length(TRANSECT_URCHIN_DENSITY_NO_M2)
    )

  # Initiate data columns
  transect_density_taxon$TEST_SIZE_MEDIAN_MM <- NA
  transect_density_taxon$TRANSECT_EROSION_G_M2_YR <- NA

  # Set median test size for each test size bin, any urchin >120mm is set to 110mm following IPRB v 1.3
  transect_density_taxon$TEST_SIZE_MEDIAN_MM[transect_density_taxon$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_0_20_MM"] <-
    10
  transect_density_taxon$TEST_SIZE_MEDIAN_MM[transect_density_taxon$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_21_40_MM"] <-
    30
  transect_density_taxon$TEST_SIZE_MEDIAN_MM[transect_density_taxon$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_41_60_MM"] <-
    50
  transect_density_taxon$TEST_SIZE_MEDIAN_MM[transect_density_taxon$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_61_80_MM"] <-
    70
  transect_density_taxon$TEST_SIZE_MEDIAN_MM[transect_density_taxon$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_81_100_MM"] <-
    90
  transect_density_taxon$TEST_SIZE_MEDIAN_MM[transect_density_taxon$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_101_120_MM"] <-
    110
  transect_density_taxon$TEST_SIZE_MEDIAN_MM[transect_density_taxon$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_121_140_MM"] <-
    110
  transect_density_taxon$TEST_SIZE_MEDIAN_MM[transect_density_taxon$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_141_160_MM"] <-
    110
  transect_density_taxon$TEST_SIZE_MEDIAN_MM[transect_density_taxon$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_161_180_MM"] <-
    110
  transect_density_taxon$TEST_SIZE_MEDIAN_MM[transect_density_taxon$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_181_200_MM"] <-
    110

  # Calculate group-specific erosion rate
  for (i in sjmisc::seq_row(transect_density_taxon)) {
    if (transect_density_taxon$TAXON_CODE[i] %in% c("ECMA", "ECOB", "ECSP") == TRUE) {
      transect_density_taxon$TRANSECT_EROSION_G_M2_YR[i] <- 0.0003 *
        (transect_density_taxon$TEST_SIZE_MEDIAN_MM[i] ^ 1.9671) *
        transect_density_taxon$TRANSECT_URCHIN_DENSITY_NO_M2[i] * 365

    }
    if (transect_density_taxon$TAXON_CODE[i] %in% c("DISP", "ECTH") == TRUE) {
      transect_density_taxon$TRANSECT_EROSION_G_M2_YR[i] <- 0.000003 *
        (transect_density_taxon$TEST_SIZE_MEDIAN_MM[i] ^ 3.2887) *
        transect_density_taxon$TRANSECT_URCHIN_DENSITY_NO_M2[i] * 365

    }
    if (transect_density_taxon$TAXON_CODE[i] %in% c("ECST", "PAGR", "EUME") == TRUE) {
      transect_density_taxon$TRANSECT_EROSION_G_M2_YR[i] <- 0.00004 *
        (transect_density_taxon$TEST_SIZE_MEDIAN_MM[i] ^ 2.6025) *
        transect_density_taxon$TRANSECT_URCHIN_DENSITY_NO_M2[i] * 365

    }
  }

  # Summary transect and site level erosion
  transect_erosion <- transect_density_taxon %>%
    group_by(
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      LATITUDE,
      LONGITUDE,
      DEPTH_M,
      LOCALDATE,
      CB_TRANSECTID,
      TRANSECT_LENGTH_M
    ) %>%
    summarize(URCHIN_ABUNDANCE_NO = sum(TRANSECT_URCHIN_ABUNDANCE_NO),
      URCHIN_EROSION_KG_M2_YR = sum(TRANSECT_EROSION_G_M2_YR) / 1000)

  transect_erosion$OCC_SITEID_TRANSECT <-
    paste(transect_erosion$OCC_SITEID,
          transect_erosion$CB_TRANSECTID,
          sep = "-")

  transect_erosion$URCHIN_DENSITY_NO_M2 <- transect_erosion$URCHIN_ABUNDANCE_NO / transect_erosion$TRANSECT_LENGTH_M

  transect_erosion_taxon <- transect_density_taxon %>%
    group_by(
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      LATITUDE,
      LONGITUDE,
      DEPTH_M,
      LOCALDATE,
      CB_TRANSECTID,
      TAXON_CODE,
      TAXON_NAME,
      TRANSECT_LENGTH_M
    ) %>%
    summarize(URCHIN_ABUNDANCE_NO = sum(TRANSECT_URCHIN_ABUNDANCE_NO),
              URCHIN_EROSION_KG_M2_YR = sum(TRANSECT_EROSION_G_M2_YR) / 1000)

  transect_erosion_taxon$URCHIN_DENSITY_NO_M2 <- transect_erosion_taxon$URCHIN_ABUNDANCE_NO / transect_erosion_taxon$TRANSECT_LENGTH_M

  site_erosion <- transect_erosion %>%
    group_by(
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      LATITUDE,
      LONGITUDE,
      DEPTH_M,
      LOCALDATE
    ) %>%
    summarize(
      URCHIN_EROSION_KG_M2_YR_MEAN = mean(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE),
      URCHIN_EROSION_KG_M2_YR_SD = sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE),
      URCHIN_EROSION_KG_M2_YR_N = length(URCHIN_EROSION_KG_M2_YR),
      URCHIN_EROSION_KG_M2_YR_SE =
        sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE) /
        length(URCHIN_EROSION_KG_M2_YR),
      URCHIN_EROSION_KG_M2_YR_CI = 1.97 *
        sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE) /
        length(URCHIN_EROSION_KG_M2_YR),
      URCHIN_DENSITY_NO_M2_MEAN = mean(URCHIN_DENSITY_NO_M2, na.rm = TRUE),
      URCHIN_DENSITY_NO_M2_SD = sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE),
      URCHIN_DENSITY_NO_M2_N = length(URCHIN_DENSITY_NO_M2),
      URCHIN_DENSITY_NO_M2_SE =
        sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE) /
        length(URCHIN_DENSITY_NO_M2),
      URCHIN_DENSITY_NO_M2_CI = 1.97 *
        sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE) /
        length(URCHIN_DENSITY_NO_M2),
      URCHIN_ABUNDANCE_NO = sum(URCHIN_ABUNDANCE_NO)
    )

  site_erosion_taxon <- transect_erosion_taxon %>%
    group_by(
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      LATITUDE,
      LONGITUDE,
      DEPTH_M,
      LOCALDATE,
      TAXON_CODE,
      TAXON_NAME
    ) %>%
    summarize(
      URCHIN_EROSION_KG_M2_YR_MEAN = mean(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE),
      URCHIN_EROSION_KG_M2_YR_SD = sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE),
      URCHIN_EROSION_KG_M2_YR_N = length(URCHIN_EROSION_KG_M2_YR),
      URCHIN_EROSION_KG_M2_YR_SE =
        sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE) /
        length(URCHIN_EROSION_KG_M2_YR),
      URCHIN_EROSION_KG_M2_YR_CI = 1.97 *
        sd(URCHIN_EROSION_KG_M2_YR, na.rm = TRUE) /
        length(URCHIN_EROSION_KG_M2_YR),
      URCHIN_DENSITY_NO_M2_MEAN = mean(URCHIN_DENSITY_NO_M2, na.rm = TRUE),
      URCHIN_DENSITY_NO_M2_SD = sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE),
      URCHIN_DENSITY_NO_M2_N = length(URCHIN_DENSITY_NO_M2),
      URCHIN_DENSITY_NO_M2_SE =
        sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE) /
        length(URCHIN_DENSITY_NO_M2),
      URCHIN_DENSITY_NO_M2_CI = 1.97 *
        sd(URCHIN_DENSITY_NO_M2, na.rm = TRUE) /
        length(URCHIN_DENSITY_NO_M2),
      URCHIN_ABUNDANCE_NO = sum(URCHIN_ABUNDANCE_NO)
    )
  
  data <- data[c(
    "REGION",
    "REGIONCODE",
    "YEAR",
    "CRUISE_ID",
    "LOCATION",
    "LOCATIONCODE",
    "OCC_SITEID",
    "LATITUDE",
    "LONGITUDE",
    "DEPTH_M",
    "LOCALDATE",
    "CB_TRANSECTID",
    "TRANSECT_LENGTH_M",
    "URCH_OBS_TF",
    "TAXON_NAME",
    "TAXON_CODE",
    "TEST_SIZE_BIN_0_20_MM",
    "TEST_SIZE_BIN_21_40_MM",
    "TEST_SIZE_BIN_41_60_MM",
    "TEST_SIZE_BIN_61_80_MM",
    "TEST_SIZE_BIN_81_100_MM",
    "TEST_SIZE_BIN_101_120_MM",
    "TEST_SIZE_BIN_121_140_MM",
    "TEST_SIZE_BIN_141_160_MM",
    "TEST_SIZE_BIN_161_180_MM",
    "TEST_SIZE_BIN_181_200_MM"
  )]


  if (full_summary == TRUE) {
    return(
      list(
        site_erosion = site_erosion,
        transect_erosion = transect_erosion,
        site_taxon = site_erosion_taxon,
        transect_taxon = transect_erosion_taxon,
        data = data
      )
    )
  }
  if (full_summary == FALSE) {
    return(site_erosion = site_erosion)
  }
}
