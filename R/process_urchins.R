#' Calculate urchin erosion rates from urchin census data
#'
#'@author Hannah Barkley
#'
#'@param data Urchin observation data set.
#'@param transect_id String of transect names (e.g., ("A1", "A2", "A3", "B1", "B2", "B3")). Defaults to NULL.
#'@param transect_length String of transect lengths in meters (e.g., c(10, 10, 10, 10, 10, 10)). Defaults to NULL.
#'@param method_name Transect design by which data were collected ("IPRB" or "Chords").
#'
#'@import tidyr
#'@import dplyr
#'@importFrom sjmisc seq_row
#'
#'@export process_urchins

process_urchins <- function(data,
                            transect_id = NULL,
                            transect_length = NULL,
                            method_name = c("IPRB", "Chords"),
                            full_summary = TRUE) {

  options(dplyr.summarise.inform = FALSE)

  # Create data frame with transect info
  if (is.null(transect_id) == TRUE &
      is.null(transect_length) == TRUE) {
    transect_summary_pairs <-
      unique(data[c("CB_TRANSECTID", "TRANSECT_LENGTH_M")])
    transect_id <- unique(transect_summary_pairs$CB_TRANSECTID)
    transect_length <-
      unique(transect_summary_pairs$TRANSECT_LENGTH_M )
  }
    summary_transect <- data.frame(transect_id, transect_length)


  # Convert to long format
  data_long <- data %>%
    select(-c(URCH_OBS_TF)) %>%
    #filter(is.na(TAXON_CODE) == FALSE) %>%
    pivot_longer(
      cols = !c(
        "REGION",
        "REGIONCODE",
        "YEAR",
        "CRUISE_ID",
        "LOCATION",
        "LOCATIONCODE",
        "OCC_SITEID",
        "OCC_SITENAME",
        "LATITUDE",
        "LONGITUDE",
        "DEPTH_M",
        "LOCALDATE",
        "CB_METHOD",
        "CB_TRANSECTID",
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
  all_bins <- c("TEST_SIZE_BIN_0_20_MM",
                "TEST_SIZE_BIN_0_20_MM",
                "TEST_SIZE_BIN_41_60_MM",
                "TEST_SIZE_BIN_61_80_MM",
                "TEST_SIZE_BIN_81_100_MM",
                "TEST_SIZE_BIN_101_120_MM",
                "TEST_SIZE_BIN_121_140_MM",
                "TEST_SIZE_BIN_141_160_MM")

  full_table <- expand.grid(
    OCC_SITEID = unique(data$OCC_SITEID),
    CB_TRANSECTID = transect_id,
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
    ),
    all.x = TRUE
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
        "OCC_SITENAME",
        "LATITUDE",
        "LONGITUDE",
        "DEPTH_M",
        "CB_METHOD",
        "TRANSECT_LENGTH_M"
      ),
      .direction = "downup"
    ) %>% group_by(OCC_SITEID, CB_METHOD) %>%
    fill(c("LOCALDATE"),
         .direction = "downup")



  data_full$TAXON_NAME[is.na(data_full$TAXON_NAME)] <- data_full$TAXON_CODE[is.na(data_full$TAXON_NAME)]

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

  # Assign urchin taxa to equation groups
  data_full$GROUP <- NA

  data_full$GROUP[data_full$TAXON_CODE %in%
                    c("ECMA", "ECOB", "ECSP")] <- "ECMA"
  data_full$GROUP[data_full$TAXON_CODE %in%
                    c("DISP", "ECTH")] <- "DISP-ECTH"
  data_full$GROUP[data_full$TAXON_CODE %in%
                    c("ECST", "PAGR")] <- "OTHER"
  data_full$GROUP[data_full$TAXON_CODE %in%
                    c("TRGR", "HEMA", "EUME", "MEGL")] <- "NON"

  # Add transect length to data frame
  data_full$TRANSECT_LENGTH_M <-
    summary_transect$transect_length[match(
      data_full$CB_TRANSECTID, summary_transect$transect_id)]

  data$TRANSECT_LENGTH_M <-
    summary_transect$transect_length[match(
      data$CB_TRANSECTID, summary_transect$transect_id)]

  # Remove non-eroding urchin taxa
  data_eroders <- subset(data_full, data_full$GROUP != "NON")

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
      "TEST_SIZE_BIN_141_160_MM"
    )
  )

  data_eroders$GROUP <- factor(data_eroders$GROUP,
                               levels = c("ECMA", "DISP-ECTH", "OTHER"))

  # Calculate transect-level density and abundance
  transect_density_group <- data_eroders %>%
    group_by(
      .data$REGION,
      .data$REGIONCODE,
      .data$YEAR,
      .data$CRUISE_ID,
      .data$LOCATION,
      .data$LOCATIONCODE,
      .data$OCC_SITEID,
      .data$OCC_SITENAME,
      .data$LATITUDE,
      .data$LONGITUDE,
      .data$DEPTH_M,
      .data$LOCALDATE,
      .data$CB_METHOD,
      .data$CB_TRANSECTID,
      .data$GROUP,
      .data$TEST_SIZE_BIN_MM,
    ) %>%
    summarize(
      TRANSECT_LENGTH_M = mean(.data$TRANSECT_LENGTH_M),
      TRANSECT_URCHIN_ABUNDANCE = sum(.data$COUNT),
      TRANSECT_URCHIN_DENSITY = sum(.data$COUNT) / mean(.data$TRANSECT_LENGTH_M)
    )

  # Calculate site level density
  site_density_group <- transect_density_group %>%
    group_by(
      .data$REGION,
      .data$REGIONCODE,
      .data$YEAR,
      .data$CRUISE_ID,
      .data$LOCATION,
      .data$LOCATIONCODE,
      .data$OCC_SITEID,
      .data$OCC_SITENAME,
      .data$LATITUDE,
      .data$LONGITUDE,
      .data$DEPTH_M,
      .data$LOCALDATE,
      .data$CB_METHOD,
      .data$GROUP,
      .data$TEST_SIZE_BIN_MM
    ) %>%
    summarize(
      SITE_URCHIN_DENSITY_MEAN = mean(.data$TRANSECT_URCHIN_DENSITY),
      SITE_URCHIN_DENSITY_SD = sd(.data$TRANSECT_URCHIN_DENSITY)
    )

  # Initiate data columns
  transect_density_group$TEST_SIZE_MEDIAN <- NA
  transect_density_group$TRANSECT_EROSION <- NA

  # Set median test size for each test size bin
  transect_density_group$TEST_SIZE_MEDIAN[
    transect_density_group$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_0_20_MM"] <- 10
  transect_density_group$TEST_SIZE_MEDIAN[
    transect_density_group$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_21_40_MM"] <- 30
  transect_density_group$TEST_SIZE_MEDIAN[
    transect_density_group$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_41_60_MM"] <- 50
  transect_density_group$TEST_SIZE_MEDIAN[
    transect_density_group$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_61_80_MM"] <- 70
  transect_density_group$TEST_SIZE_MEDIAN[
    transect_density_group$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_81_100_MM"] <- 90
  transect_density_group$TEST_SIZE_MEDIAN[
    transect_density_group$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_101_120_MM"] <- 110
  transect_density_group$TEST_SIZE_MEDIAN[
    transect_density_group$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_121_140_MM"] <- 130
  transect_density_group$TEST_SIZE_MEDIAN[
    transect_density_group$TEST_SIZE_BIN_MM == "TEST_SIZE_BIN_141_160_MM"] <- 150


  # Calculate group-specific erosion rate
  for (i in sjmisc::seq_row(transect_density_group)) {
    if (transect_density_group$GROUP[i] == "ECMA") {
      transect_density_group$TRANSECT_EROSION[i] <- 0.0003 *
        (transect_density_group$TEST_SIZE_MEDIAN[i] ^ 1.9671) *
        transect_density_group$TRANSECT_URCHIN_DENSITY[i] * 365

    }
    if (transect_density_group$GROUP[i] == "DISP-ECTH") {
      transect_density_group$TRANSECT_EROSION[i] <- 0.000003 *
        (transect_density_group$TEST_SIZE_MEDIAN[i] ^ 3.2887) *
        transect_density_group$TRANSECT_URCHIN_DENSITY[i] * 365

    }
    if (transect_density_group$GROUP[i] == "OTHER") {
      transect_density_group$TRANSECT_EROSION[i] <- 0.00004 *
        (transect_density_group$TEST_SIZE_MEDIAN[i] ^ 2.6025) *
        transect_density_group$TRANSECT_URCHIN_DENSITY[i] * 365

    }
  }

  # Summary transect and site level erosion
  transect_density_sum <- transect_density_group %>%
    group_by(
      .data$REGION,
      .data$REGIONCODE,
      .data$YEAR,
      .data$CRUISE_ID,
      .data$LOCATION,
      .data$LOCATIONCODE,
      .data$OCC_SITEID,
      .data$OCC_SITENAME,
      .data$LATITUDE,
      .data$LONGITUDE,
      .data$DEPTH_M,
      .data$LOCALDATE,
      .data$CB_METHOD,
      .data$CB_TRANSECTID,
    ) %>%
    summarize(URCHIN_EROSION_KG_M2_YR = sum(.data$TRANSECT_EROSION) / 1000)

  transect_density_sum$OCC_SITEID_TRANSECT <-
    paste(transect_density_sum$OCC_SITEID,
          transect_density_sum$CB_TRANSECTID,
          sep = "-")

  site_erosion <- transect_density_sum %>%
    group_by(  .data$REGION,
               .data$REGIONCODE,
               .data$YEAR,
               .data$CRUISE_ID,
               .data$LOCATION,
               .data$LOCATIONCODE,
               .data$OCC_SITEID,
               .data$OCC_SITENAME,
               .data$LATITUDE,
               .data$LONGITUDE,
               .data$DEPTH_M,
               .data$LOCALDATE,
               .data$CB_METHOD) %>%
    summarize(
      URCHIN_EROSION_KG_M2_YR_MEAN = mean(.data$URCHIN_EROSION_KG_M2_YR),
      URCHIN_EROSION_KG_M2_YR_SD = sd(.data$URCHIN_EROSION_KG_M2_YR),
      URCHIN_EROSION_KG_M2_YR_N = length(.data$URCHIN_EROSION_KG_M2_YR),
      URCHIN_EROSION_KG_M2_YR_SE =
        sd(.data$URCHIN_EROSION_KG_M2_YR) /
        length(.data$URCHIN_EROSION_KG_M2_YR),
      URCHIN_EROSION_KG_M2_YR_CI = 1.97 *
        sd(.data$URCHIN_EROSION_KG_M2_YR) /
        length(.data$URCHIN_EROSION_KG_M2_YR)
    )

  data <- data[c(
  "REGION",
  "REGIONCODE",
  "YEAR",
  "CRUISE_ID",
  "LOCATION",
  "LOCATIONCODE",
  "OCC_SITEID",
  "OCC_SITENAME",
  "LATITUDE",
  "LONGITUDE",
  "DEPTH_M",
  "LOCALDATE",
  "CB_METHOD",
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
  "TEST_SIZE_BIN_141_160_MM")]


  if (full_summary == TRUE) {
    return(
      list(
        site_erosion = site_erosion,
        transect_erosion = transect_density_sum,
        data = data
      )
    )
  }
  if (full_summary == FALSE) {
    return(site_erosion = site_erosion)
  }
}
