#' Summarize carbonate production rates at transect and site level
#'@author Hannah Barkley
#'@param data Data set to summarize; product of run_calc_prod
#'@param transect_summary Transect sumamry; product of run_calc_prod
#'@param dbase_type Production database to use ("IPRB" or "NCRMP")
#'@param summarize_by Grouping factor to summarize by ("substrate code",
#'"substrate class", "coral group",or "overall")
#'@param level Summarize at "transect" or "site" level
#'@param macro_rate Rate of macrobioerosion
#'@param macro_rate_ci Confidence interval for rate of macrobioerosion
#'@param micro_rate Rate of microbioerosion
#'@param micro_rate_ci Confidence interval for rate of microbioerosion
#'@import dplyr
#'@importFrom sjmisc seq_row
#'@importFrom rlang .data
#'@export summarize_prod

summarize_prod <- function(data,
                           transect_summary,
                           dbase_type = c("IPRB", "NCRMP"),
                           summarize_by = c("substrate code",
                                            "substrate class",
                                            "coral group",
                                            "overall"),
                           level = c("transect", "site"),
                           macro_rate = 0.209,
                           macro_rate_ci = 0.129,
                           micro_rate = 0.262,
                           micro_rate_ci = 0.180,
                           ...

) {

  if (dbase_type == "IPRB") {
    data$SUBSTRATE_CODE <- data$SUBSTRATE_CODE_IPRB
    prod_dbase <- prod_dbase_iprb
  }

  if (dbase_type == "NCRMP") {
    prod_dbase <- prod_dbase_ncrmp
  }

  # Summarize production by SUBSTRATE_CODE at TRANSECT level  ----------------

  summary_transect_substratecode <- suppressMessages(
    data %>% dplyr::group_by(
      .data$REGION,
      .data$REGIONCODE,
      .data$YEAR,
      .data$CRUISE_ID,
      .data$LOCATION,
      .data$LOCATIONCODE,
      .data$OCC_SITEID,
      .data$OCC_SITENAME,
      .data$LATITUTE,
      .data$LONGITUDE,
      .data$DEPTH_M,
      .data$LOCALDATE,
      .data$CB_METHOD,
      .data$CB_TRANSECTID,
      .data$OCC_SITEID_TRANSECT,
      .data$SUBSTRATE_CLASS,
      .data$SUBSTRATE_NAME,
      .data$SUBSTRATE_CODE,
      .data$MORPHOLOGY,
      .data$MORPHOLOGYCODE
    ) %>%
      summarize(
        TRANSECT_PLANAR_LENGTH_M =
          first(.data$TRANSECT_PLANAR_LENGTH_M),
        TRANSECT_TOTAL_SUBSTRATE_COVER_M =
          first(.data$TRANSECT_TOTAL_SUBSTRATE_COVER_M),
        SUBSTRATE_COVER_CM =
          sum(.data$SUBSTRATE_COVER_CM, na.rm = TRUE),
        SUBSTRATE_COVER_PCT =
          .data$SUBSTRATE_COVER_CM /
          .data$TRANSECT_TOTAL_SUBSTRATE_COVER_M,
        SUBSTRATE_PLANAR_PROD_KG_CM_YR =
          sum(.data$COLONY_PROD_G_YR, na.rm = TRUE) /
          (.data$TRANSECT_PLANAR_LENGTH_M * 100),
        SUBSTRATE_PLANAR_PROD_KG_CM_YR_L95 =
          sum(.data$COLONY_PROD_G_YR_L95, na.rm = TRUE) /
          (.data$TRANSECT_PLANAR_LENGTH_M * 100),
        SUBSTRATE_PLANAR_PROD_KG_CM_YR_U95 =
          sum(.data$COLONY_PROD_G_YR_U95, na.rm = TRUE) /
          (.data$TRANSECT_PLANAR_LENGTH_M * 100),
        SUBSTRATE_CARB_PROD_KG_M2_YR =
          .data$SUBSTRATE_PLANAR_PROD_KG_CM_YR * 10,
        SUBSTRATE_CARB_PROD_KG_M2_YR_L95 =
          .data$SUBSTRATE_PLANAR_PROD_KG_CM_YR_L95 * 10,
        SUBSTRATE_CARB_PROD_KG_M2_YR_U95 =
          .data$SUBSTRATE_PLANAR_PROD_KG_CM_YR_U95 * 10
      )
  )

  # Summarize production by SUBSTRATE_CLASS at TRANSECT level  ---------------
  transect_substrateclass <- suppressMessages(
    data %>% dplyr::group_by(
      .data$REGION,
      .data$REGIONCODE,
      .data$YEAR,
      .data$CRUISE_ID,
      .data$LOCATION,
      .data$LOCATIONCODE,
      .data$OCC_SITEID,
      .data$OCC_SITENAME,
      .data$LATITUTE,
      .data$LONGITUDE,
      .data$DEPTH_M,
      .data$LOCALDATE,
      .data$CB_METHOD,
      .data$CB_TRANSECTID,
      .data$OCC_SITEID_TRANSECT,
      .data$SUBSTRATE_CLASS
    ) %>%
      summarize(
        TRANSECT_PLANAR_LENGTH_M = first(.data$TRANSECT_PLANAR_LENGTH_M),
        TRANSECT_TOTAL_SUBSTRATE_COVER_M =
          first(.data$TRANSECT_TOTAL_SUBSTRATE_COVER_M),
        SUBSTRATE_COVER_CM =
          sum(.data$SUBSTRATE_COVER_CM, na.rm = TRUE),
        SUBSTRATE_COVER_PCT =
          .data$SUBSTRATE_COVER_CM / .data$TRANSECT_TOTAL_SUBSTRATE_COVER_M,
        SUBSTRATE_PLANAR_PROD_KG_CM_YR =
          sum(.data$COLONY_PROD_G_YR, na.rm = TRUE) /
          (.data$TRANSECT_PLANAR_LENGTH_M * 100),
        SUBSTRATE_PLANAR_PROD_KG_CM_YR_L95 =
          sum(.data$COLONY_PROD_G_YR_L95, na.rm = TRUE) /
          (.data$TRANSECT_PLANAR_LENGTH_M * 100),
        SUBSTRATE_PLANAR_PROD_KG_CM_YR_U95 =
          sum(.data$COLONY_PROD_G_YR_U95, na.rm = TRUE) /
          (.data$TRANSECT_PLANAR_LENGTH_M * 100),
        SUBSTRATE_CARB_PROD_KG_M2_YR =
          .data$SUBSTRATE_PLANAR_PROD_KG_CM_YR * 10,
        SUBSTRATE_CARB_PROD_KG_M2_YR_L95 =
          .data$SUBSTRATE_PLANAR_PROD_KG_CM_YR_L95 * 10,
        SUBSTRATE_CARB_PROD_KG_M2_YR_U95 =
          .data$SUBSTRATE_PLANAR_PROD_KG_CM_YR_U95 * 10
      )
  )

  # Create a table of all substrate classes for each transect
  substrate_class_all <- c("CORAL",
                           "CCA",
                           "CARB",
                           "MA",
                           "RCK",
                           "SAND",
                           "TURF",
                           "OTHER")

  substrate_full_table <- expand.grid(
    OCC_SITEID = unique(transect_summary$OCC_SITEID),
    CB_TRANSECTID = unique(transect_summary$CB_TRANSECTID),
    SUBSTRATE_CLASS = substrate_class_all
  )

  substrate_full_table$OCC_SITEID_TRANSECT <-
    paste(substrate_full_table$OCC_SITEID,
          substrate_full_table$CB_TRANSECTID,
          sep = "-")

  # Populate full table with transect-level data
  summary_transect_substrateclass <-  merge(
    substrate_full_table,
    transect_substrateclass,
    by = c(
      "OCC_SITEID",
      "CB_TRANSECTID",
      "OCC_SITEID_TRANSECT",
      "SUBSTRATE_CLASS"
    ),
    all.x = TRUE
  )

  summary_transect_substrateclass <- summary_transect_substrateclass[c(
    "REGION",
    "REGIONCODE",
    "YEAR",
    "CRUISE_ID",
    "LOCATION",
    "LOCATIONCODE",
    "OCC_SITEID",
    "OCC_SITENAME",
    "LATITUTE",
    "LONGITUDE",
    "DEPTH_M",
    "LOCALDATE",
    "CB_METHOD",
    "CB_TRANSECTID",
    "OCC_SITEID_TRANSECT",
    "SUBSTRATE_CLASS",
    "TRANSECT_PLANAR_LENGTH_M",
    "TRANSECT_TOTAL_SUBSTRATE_COVER_M",
    "SUBSTRATE_COVER_CM",
    "SUBSTRATE_COVER_PCT",
    "SUBSTRATE_PLANAR_PROD_KG_CM_YR",
    "SUBSTRATE_PLANAR_PROD_KG_CM_YR_L95",
    "SUBSTRATE_PLANAR_PROD_KG_CM_YR_U95",
    "SUBSTRATE_CARB_PROD_KG_M2_YR",
    "SUBSTRATE_CARB_PROD_KG_M2_YR_L95",
    "SUBSTRATE_CARB_PROD_KG_M2_YR_U95"
  )]


  summary_transect_substrateclass <-
    summary_transect_substrateclass %>%
    group_by(OCC_SITEID_TRANSECT) %>%
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
        "LATITUTE",
        "LONGITUDE",
        "DEPTH_M",
        "LOCALDATE",
        "CB_METHOD"
      )
    )

  summary_transect_substrateclass[, 17:26][is.na(summary_transect_substrateclass[, 17:26])] <- 0

  # Summarize production by CORAL_GROUP at TRANSECT level  -------------------
  transect_coral_group <- data %>% dplyr::group_by(
    .data$REGION,
    .data$REGIONCODE,
    .data$YEAR,
    .data$CRUISE_ID,
    .data$LOCATION,
    .data$LOCATIONCODE,
    .data$OCC_SITEID,
    .data$OCC_SITENAME,
    .data$LATITUTE,
    .data$LONGITUDE,
    .data$DEPTH_M,
    .data$LOCALDATE,
    .data$CB_METHOD,
    .data$CB_TRANSECTID,
    .data$OCC_SITEID_TRANSECT,
    .data$CORAL_GROUP,
    .data$CORAL_GROUP_NAME,
  ) %>%
    summarize(
      TRANSECT_PLANAR_LENGTH_M =
        first(.data$TRANSECT_PLANAR_LENGTH_M),
      TRANSECT_TOTAL_SUBSTRATE_COVER_M =
        first(.data$TRANSECT_TOTAL_SUBSTRATE_COVER_M),
      SUBSTRATE_COVER_CM =
        sum(.data$SUBSTRATE_COVER_CM, na.rm = TRUE),
      SUBSTRATE_COVER_PCT =
        .data$SUBSTRATE_COVER_CM / .data$TRANSECT_TOTAL_SUBSTRATE_COVER_M,
      SUBSTRATE_PLANAR_PROD_KG_CM_YR =
        sum(.data$COLONY_PROD_G_YR, na.rm = TRUE) /
        (.data$TRANSECT_PLANAR_LENGTH_M * 100),
      SUBSTRATE_PLANAR_PROD_KG_CM_YR_L95 =
        sum(.data$COLONY_PROD_G_YR_L95, na.rm = TRUE) /
        (.data$TRANSECT_PLANAR_LENGTH_M * 100),
      SUBSTRATE_PLANAR_PROD_KG_CM_YR_U95 =
        sum(.data$COLONY_PROD_G_YR_U95, na.rm = TRUE) /
        (.data$TRANSECT_PLANAR_LENGTH_M * 100),
      SUBSTRATE_CARB_PROD_KG_M2_YR =
        .data$SUBSTRATE_PLANAR_PROD_KG_CM_YR * 10,
      SUBSTRATE_CARB_PROD_KG_M2_YR_L95 =
        .data$SUBSTRATE_PLANAR_PROD_KG_CM_YR_L95 * 10,
      SUBSTRATE_CARB_PROD_KG_M2_YR_U95 =
        .data$SUBSTRATE_PLANAR_PROD_KG_CM_YR_U95 * 10
    )

  transect_coral_group <-
    subset(transect_coral_group,
           transect_coral_group$CORAL_GROUP != "NA")

  coral_group_all <- c(
    "ACBR",
    "ASMA",
    "BR",
    "ENC",
    "FOL",
    "FREE",
    "GOMA",
    "MASS",
    "MOBR",
    "MOEN",
    "PAEN",
    "PAMA",
    "POBR",
    "POCS",
    "POEN",
    "POLC",
    "POMA"
  )

  coral_full_table <- expand.grid(
    OCC_SITEID = unique(transect_summary$OCC_SITEID),
    CB_TRANSECTID = unique(transect_summary$CB_TRANSECTID),
    CORAL_GROUP = coral_group_all
  )

  coral_full_table$OCC_SITEID_TRANSECT <-
    paste(coral_full_table$OCC_SITEID,
          coral_full_table$CB_TRANSECTID,
          sep = "-")

  # Populate full table with transect-level data
  summary_transect_coral <-  merge(
    coral_full_table,
    transect_coral_group,
    by = c(
      "OCC_SITEID",
      "CB_TRANSECTID",
      "OCC_SITEID_TRANSECT",
      "CORAL_GROUP"
    ),
    all.x = TRUE
  )

  summary_transect_coral <- summary_transect_coral[c(
    "REGION",
    "REGIONCODE",
    "YEAR",
    "CRUISE_ID",
    "LOCATION",
    "LOCATIONCODE",
    "OCC_SITEID",
    "OCC_SITENAME",
    "LATITUTE",
    "LONGITUDE",
    "DEPTH_M",
    "LOCALDATE",
    "CB_METHOD",
    "CB_TRANSECTID",
    "OCC_SITEID_TRANSECT",
    "CORAL_GROUP",
    "CORAL_GROUP_NAME",
    "TRANSECT_PLANAR_LENGTH_M",
    "TRANSECT_TOTAL_SUBSTRATE_COVER_M",
    "SUBSTRATE_COVER_CM",
    "SUBSTRATE_COVER_PCT",
    "SUBSTRATE_PLANAR_PROD_KG_CM_YR",
    "SUBSTRATE_PLANAR_PROD_KG_CM_YR_L95",
    "SUBSTRATE_PLANAR_PROD_KG_CM_YR_U95",
    "SUBSTRATE_CARB_PROD_KG_M2_YR",
    "SUBSTRATE_CARB_PROD_KG_M2_YR_L95",
    "SUBSTRATE_CARB_PROD_KG_M2_YR_U95"
  )]

  summary_transect_coral <-
    summary_transect_coral %>% group_by(OCC_SITEID_TRANSECT) %>%
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
        "LATITUTE",
        "LONGITUDE",
        "DEPTH_M",
        "LOCALDATE",
        "CB_METHOD"
      ),
      .direction = "downup"
    )

  summary_transect_coral <-
    summary_transect_coral %>% group_by(CORAL_GROUP ) %>%
    fill(
      c(
        "CORAL_GROUP_NAME"
      ),
      .direction = "downup"
    )

  summary_transect_coral[, 18:27][is.na(summary_transect_coral[, 18:27])] <- 0



  # Summarize production by TRANSECT ----------------------------------------
  if (dbase_type == "IPRB") {
    summary_transect <- summary_transect_substratecode  %>%
      dplyr::group_by(
        .data$REGION,
        .data$REGIONCODE,
        .data$YEAR,
        .data$CRUISE_ID,
        .data$LOCATION,
        .data$LOCATIONCODE,
        .data$OCC_SITEID,
        .data$OCC_SITENAME,
        .data$LATITUTE,
        .data$LONGITUDE,
        .data$DEPTH_M,
        .data$LOCALDATE,
        .data$CB_METHOD,
        .data$CB_TRANSECTID,
        .data$OCC_SITEID_TRANSECT,
      ) %>%
      summarize(
        TRANSECT_PLANAR_LENGTH_M =
          (.data$TRANSECT_PLANAR_LENGTH_M),
        TRANSECT_TOTAL_SUBSTRATE_COVER_M =
          mean(.data$TRANSECT_TOTAL_SUBSTRATE_COVER_M),
        GROSS_CARB_PROD_KG_M2_YR =
          sum(.data$SUBSTRATE_CARB_PROD_KG_M2_YR, na.rm = TRUE),
        GROSS_CARB_PROD_KG_M2_YR_L95 =
          sum(.data$SUBSTRATE_CARB_PROD_KG_M2_YR_L95, na.rm = TRUE),
        GROSS_CARB_PROD_KG_M2_YR_U95 =
          sum(.data$SUBSTRATE_CARB_PROD_KG_M2_YR_U95, na.rm = TRUE),
        RUGOSITY =
          .data$TRANSECT_TOTAL_SUBSTRATE_COVER_M /
          .data$TRANSECT_PLANAR_LENGTH_M,
        SUBSTRATE_AVAILABLE_MACRO_PCT =
          sum(.data$SUBSTRATE_COVER_PCT[.data$SUBSTRATE_CODE %in% c("CCA",
                                                               "MCA",
                                                               "RUBC",
                                                               "SCA",
                                                               "BOR",
                                                               "DC",
                                                               "MAC",
                                                               "OCE",
                                                               "OTH",
                                                               "RUB",
                                                               "RUBT",
                                                               "TF")]),
        SUBSTRATE_AVAILABLE_MACRO_INDEX =
          .data$RUGOSITY * (.data$SUBSTRATE_AVAILABLE_MACRO_PCT / 100),
        SUBSTRATE_AVAILABLE_MICRO_PCT =
          100 - sum(.data$SUBSTRATE_COVER_PCT[.data$SUBSTRATE_CODE %in%
                                                 c("SG", "RCK", "SAND")]),
        SUBSTRATE_AVAILABLE_MICRO_INDEX =
          .data$RUGOSITY * (.data$SUBSTRATE_AVAILABLE_MICRO_PCT / 100),
        HARD_CORAL_COVER_PCT =
          sum(.data$SUBSTRATE_COVER_PCT[.data$SUBSTRATE_CLASS == "CORAL"]),
        CCA_COVER_PCT =
          sum(.data$SUBSTRATE_COVER_PCT[.data$SUBSTRATE_CLASS == "CCA"]),
        HARD_CORAL_CARB_PROD_KG_M2_YR =
          sum(.data$SUBSTRATE_CARB_PROD_KG_M2_YR[
            .data$SUBSTRATE_CLASS == "CORAL"],
              na.rm = TRUE),
        CCA_CARB_PROD_KG_M2_YR =
          sum(.data$SUBSTRATE_CARB_PROD_KG_M2_YR[
            .data$SUBSTRATE_CLASS == "CCA"],
              na.rm = TRUE),
        MACROBIOEROSION_KG_M2_YR =
          .data$SUBSTRATE_AVAILABLE_MACRO_INDEX * macro_rate,
        MACROBIOEROSION_KG_M2_YR_L95 =
          .data$SUBSTRATE_AVAILABLE_MACRO_INDEX * (macro_rate - macro_rate_ci),
        MACROBIOEROSION_KG_M2_YR_U95 =
          .data$SUBSTRATE_AVAILABLE_MACRO_INDEX * (macro_rate + macro_rate_ci),
        MICROBIOEROSION_KG_M2_YR =
          .data$SUBSTRATE_AVAILABLE_MICRO_INDEX * micro_rate,
        MICROBIOEROSION_KG_M2_YR_L95 =
          .data$SUBSTRATE_AVAILABLE_MICRO_INDEX * (micro_rate - micro_rate_ci),
        MICROBIOEROSION_KG_M2_YR_U95 =
          .data$SUBSTRATE_AVAILABLE_MICRO_INDEX * (micro_rate + micro_rate_ci),
        BIOEROSION_KG_M2_YR =
          .data$MACROBIOEROSION_KG_M2_YR +
          .data$MICROBIOEROSION_KG_M2_YR,
        BIOEROSION_KG_M2_YR_L95 =
          .data$MACROBIOEROSION_KG_M2_YR_L95 +
          .data$MICROBIOEROSION_KG_M2_YR_L95,
        BIOEROSION_KG_M2_YR_U95 =
          .data$MACROBIOEROSION_KG_M2_YR_U95 +
          .data$MICROBIOEROSION_KG_M2_YR_U95
      )
  }
  if (dbase_type == "NCRMP") {
    summary_transect <-
      summary_transect_substratecode  %>% dplyr::group_by(
        .data$REGION,
        .data$REGIONCODE,
        .data$YEAR,
        .data$CRUISE_ID,
        .data$LOCATION,
        .data$LOCATIONCODE,
        .data$OCC_SITEID,
        .data$OCC_SITENAME,
        .data$LATITUTE,
        .data$LONGITUDE,
        .data$DEPTH_M,
        .data$LOCALDATE,
        .data$CB_METHOD,
        .data$CB_TRANSECTID,
        .data$OCC_SITEID_TRANSECT,
      ) %>%
      summarize(
        TRANSECT_PLANAR_LENGTH_M =
          mean(.data$TRANSECT_PLANAR_LENGTH_M),
        TRANSECT_TOTAL_SUBSTRATE_COVER_M =
          mean(.data$TRANSECT_TOTAL_SUBSTRATE_COVER_M),
        GROSS_CARB_PROD_KG_M2_YR =
          sum(.data$SUBSTRATE_CARB_PROD_KG_M2_YR, na.rm = TRUE),
        GROSS_CARB_PROD_KG_M2_YR_L95 =
          sum(.data$SUBSTRATE_CARB_PROD_KG_M2_YR_L95, na.rm = TRUE),
        GROSS_CARB_PROD_KG_M2_YR_U95 =
          sum(.data$SUBSTRATE_CARB_PROD_KG_M2_YR_U95, na.rm = TRUE),
        RUGOSITY =
          .data$TRANSECT_TOTAL_SUBSTRATE_COVER_M /
          .data$TRANSECT_PLANAR_LENGTH_M,
        SUBSTRATE_AVAILABLE_MACRO_PCT =
          sum(.data$SUBSTRATE_COVER_PCT[.data$SUBSTRATE_CODE %in% c("HARD",
                                                                     "RUB",
                                                                     "DC",
                                                                     "CCA",
                                                                     "MA",
                                                                     "PESP",
                                                                     "TURF",
                                                                     "SP")]),
        SUBSTRATE_AVAILABLE_MACRO_INDEX =
          .data$RUGOSITY * (.data$SUBSTRATE_AVAILABLE_MACRO_PCT / 100),
        SUBSTRATE_AVAILABLE_MICRO_PCT =
          100 - sum(.data$SUBSTRATE_COVER_PCT[.data$SUBSTRATE_CODE %in%
                                                 c("SG", "RCK", "SAND")]),
        SUBSTRATE_AVAILABLE_MICRO_INDEX =
          .data$RUGOSITY * (.data$SUBSTRATE_AVAILABLE_MICRO_PCT / 100),
        HARD_CORAL_COVER_PCT =
          sum(.data$SUBSTRATE_COVER_PCT[.data$SUBSTRATE_CLASS == "CORAL"]),
        CCA_COVER_PCT =
          sum(.data$SUBSTRATE_COVER_PCT[.data$SUBSTRATE_CLASS == "CCA"]),
        HARD_CORAL_CARB_PROD_KG_M2_YR =
          sum(.data$SUBSTRATE_CARB_PROD_KG_M2_YR[
            .data$SUBSTRATE_CLASS == "CORAL"],
              na.rm = TRUE),
        CCA_CARB_PROD_KG_M2_YR =
          sum(.data$SUBSTRATE_CARB_PROD_KG_M2_YR[
            .data$SUBSTRATE_CLASS == "CCA"],
              na.rm = TRUE),
        MACROBIOEROSION_KG_M2_YR =
          .data$SUBSTRATE_AVAILABLE_MACRO_INDEX * macro_rate,
        MACROBIOEROSION_KG_M2_YR_L95 =
          .data$SUBSTRATE_AVAILABLE_MACRO_INDEX * (macro_rate - macro_rate_ci),
        MACROBIOEROSION_KG_M2_YR_U95 =
          .data$SUBSTRATE_AVAILABLE_MACRO_INDEX * (macro_rate + macro_rate_ci),
        MICROBIOEROSION_KG_M2_YR =
          .data$SUBSTRATE_AVAILABLE_MICRO_INDEX * micro_rate,
        MICROBIOEROSION_KG_M2_YR_L95 =
          .data$SUBSTRATE_AVAILABLE_MICRO_INDEX * (micro_rate - micro_rate_ci),
        MICROBIOEROSION_KG_M2_YR_U95 =
          .data$SUBSTRATE_AVAILABLE_MICRO_INDEX * (micro_rate + micro_rate_ci),
        BIOEROSION_KG_M2_YR =
          .data$MACROBIOEROSION_KG_M2_YR +
          .data$MICROBIOEROSION_KG_M2_YR,
        BIOEROSION_KG_M2_YR_L95 =
          .data$MACROBIOEROSION_KG_M2_YR_L95 +
          .data$MICROBIOEROSION_KG_M2_YR_L95,
        BIOEROSION_KG_M2_YR_U95 =
          .data$MACROBIOEROSION_KG_M2_YR_U95 +
          .data$MICROBIOEROSION_KG_M2_YR_U95
      )
  }

  # Summarize production by SUBSTRATE_CODE at SITE level -------------------

  se <- function(x) {
    sd(x, na.rm = TRUE) / sqrt(length(x))
  }

  ci95 <- function(x) {
    sd(x, na.rm = TRUE) / sqrt(length(x)) * 1.97
  }

  summary_site_substratecode <-
    summary_transect_substratecode %>% dplyr::group_by(
      .data$REGION,
      .data$REGIONCODE,
      .data$YEAR,
      .data$CRUISE_ID,
      .data$LOCATION,
      .data$LOCATIONCODE,
      .data$OCC_SITEID,
      .data$OCC_SITENAME,
      .data$LATITUTE,
      .data$LONGITUDE,
      .data$DEPTH_M,
      .data$LOCALDATE,
      .data$CB_METHOD,
      .data$SUBSTRATE_CLASS,
      .data$SUBSTRATE_CODE,
      .data$MORPHOLOGY,
      .data$MORPHOLOGYCODE
    ) %>%
    dplyr::summarize(across(
      c(SUBSTRATE_COVER_PCT,
        SUBSTRATE_CARB_PROD_KG_M2_YR),
      list(
        MEAN = ~mean(.x, na.rm = TRUE),
        SD = ~sd(.x, na.rm = TRUE),
        SE = se,
        CI95 = ci95,
        N = length
      )
    ))

  # Summarize production by SUBSTRATE_CLASS at SITE level --------------------
  summary_site_substrateclass <-
    summary_transect_substrateclass %>%
    dplyr::group_by(
      .data$REGION,
      .data$REGIONCODE,
      .data$YEAR,
      .data$CRUISE_ID,
      .data$LOCATION,
      .data$LOCATIONCODE,
      .data$OCC_SITEID,
      .data$OCC_SITENAME,
      .data$LATITUTE,
      .data$LONGITUDE,
      .data$DEPTH_M,
      .data$LOCALDATE,
      .data$CB_METHOD,
      .data$SUBSTRATE_CLASS
    ) %>%
  dplyr::summarize(across(
    c(SUBSTRATE_COVER_PCT,
      SUBSTRATE_CARB_PROD_KG_M2_YR),
    list(
      MEAN = ~mean(.x, na.rm = TRUE),
      SD = ~sd(.x, na.rm = TRUE),
      SE = se,
      CI95 = ci95,
      N = length
    )
  ))

  # Summarize production by CORAL_GROUP at SITE level -----------------------
  summary_site_coral <- summary_transect_coral %>%
    dplyr::group_by(
      .data$REGION,
      .data$REGIONCODE,
      .data$YEAR,
      .data$CRUISE_ID,
      .data$LOCATION,
      .data$LOCATIONCODE,
      .data$OCC_SITEID,
      .data$OCC_SITENAME,
      .data$LATITUTE,
      .data$LONGITUDE,
      .data$DEPTH_M,
      .data$LOCALDATE,
      .data$CB_METHOD,
      .data$CORAL_GROUP,
      .data$CORAL_GROUP_NAME
    ) %>%
    dplyr::summarize(across(
      c(SUBSTRATE_COVER_PCT,
        SUBSTRATE_CARB_PROD_KG_M2_YR),
      list(
        MEAN = ~mean(.x, na.rm = TRUE),
        SD = ~sd(.x, na.rm = TRUE),
        SE = se,
        CI95 = ci95,
        N = length
      )
    ))


  # Summarize production by SITE --------------------------------------------
  summary_site <- summary_transect %>%
    dplyr::group_by(
      .data$REGION,
      .data$REGIONCODE,
      .data$YEAR,
      .data$CRUISE_ID,
      .data$LOCATION,
      .data$LOCATIONCODE,
      .data$OCC_SITEID,
      .data$OCC_SITENAME,
      .data$LATITUTE,
      .data$LONGITUDE,
      .data$DEPTH_M,
      .data$LOCALDATE,
      .data$CB_METHOD,
    ) %>%
    dplyr::summarize(across(
      c(RUGOSITY,
        HARD_CORAL_COVER_PCT,
        CCA_COVER_PCT,
        GROSS_CARB_PROD_KG_M2_YR,
        MACROBIOEROSION_KG_M2_YR,
        MICROBIOEROSION_KG_M2_YR,
        BIOEROSION_KG_M2_YR,
        HARD_CORAL_CARB_PROD_KG_M2_YR,
        CCA_CARB_PROD_KG_M2_YR
      ),
      list(
        MEAN = ~ mean(.x, na.rm = TRUE),
        SD = ~ sd(.x, na.rm = TRUE),
        SE = se,
        CI95 = ci95,
        N = length
      )
    ))


  # Export data -------------------------------------------------------------
  if (summarize_by == "substrate code" && level == "transect") {
    return(summary_transect_substratecode)
  }
  if (summarize_by == "substrate class" && level == "transect") {
    return(summary_transect_substrateclass)
  }
  if (summarize_by == "coral group" && level == "transect") {
    return(summary_transect_coral)
  }
  if (summarize_by == "overall" && level == "transect") {
    return(summary_transect)
  }
  if (summarize_by == "substrate code" && level == "site") {
    return(summary_site_substratecode)
  }
  if (summarize_by == "substrate class" && level == "site") {
    return(summary_site_substrateclass)
  }
  if (summarize_by == "coral group" && level == "site") {
    return(summary_site_coral)
  }
  if (summarize_by == "overall" && level == "site") {
    return(summary_site)
  }
}

