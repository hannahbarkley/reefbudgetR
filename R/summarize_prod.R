#' Summarize carbonate production rates at transect and site level
#'
#'@author Hannah Barkley
#'
#'@param data Data set to summarize; product of `run_calc_prod`.
#'@param transect_summary Transect sumamry; product of `run_calc_prod`.
#'@param dbase_type Production database to use, either Indo-Pacific ReefBudget ("IPRB")
#'or U.S. Pacific Islands NCRMP-specific database ("NCRMP"). The Indo-Pacific ReefBudget
#'database is derived from "IP Calcification and bioerosion rates database v.1.3",
#'downloaded from https://geography.exeter.ac.uk/reefbudget/indopacific/.
#'@param summarize_by Grouping factor to summarize by ("substrate code",
#'"substrate class", "coral group",or "overall")
#'@param level Summarize at "transect" or "site" level.
#'@param macro_rate Rate of macrobioerosion in kg/cm2/yr. Default is 0.209.
#'@param macro_rate_ci Confidence interval for rate of macrobioerosion in kg/cm2/yr, Default is 0.129.
#'@param micro_rate Rate of microbioerosion in kg/cm2/yr. Default is 0.262.
#'@param micro_rate_ci Confidence interval for rate of microbioerosion in kg/cm2/yr. Default is 0.180.
#'
#'@import dplyr
#'@importFrom sjmisc seq_row
#'@importFrom rlang .data
#'
#'@export summarize_prod
#'
#'@examples
#' calc_prod_output <- run_calc_prod(
#'     data = data,
#'     method_name = "IPRB",
#'     dbase_type = "NCRMP")
#'
#' data <- calc_prod_output$data
#' transect_summary <- calc_prod_output$transect_summary
#'
#' prod_transect_substratecode <-
#'   summarize_prod(data,
#'                  transect_summary,
#'                  dbase_type,
#'                  summarize_by = "substrate code",
#'                  level = "transect")
#'
#' prod_transect_substrateclass <-
#'   summarize_prod(data,
#'                  transect_summary,
#'                  dbase_type,
#'                  summarize_by = "substrate class",
#'                  level = "transect")
#'
#' prod_transect_coral <-
#'   summarize_prod(data,
#'                  transect_summary,
#'                  dbase_type,
#'                  summarize_by = "coral group",
#'                  level = "transect")
#'
#' prod_transect <-
#'   summarize_prod(data,
#'                  transect_summary,
#'                  dbase_type,
#'                  summarize_by = "overall",
#'                  level = "transect")

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

  transect_substratecode <- suppressMessages(
    data %>% dplyr::group_by(
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      SITEVISITID,
      LATITUDE,
      LONGITUDE,
      SITE_DEPTH_M,
      LOCALDATE,
      CB_TRANSECTID,
      OCC_SITEID_TRANSECT,
      SUBSTRATE_CLASS,
      SUBSTRATE_NAME,
      SUBSTRATE_CODE,
      MORPHOLOGY,
      MORPHOLOGYCODE
    ) %>%
      reframe(
        TRANSECT_PLANAR_LENGTH_M =
          mean(TRANSECT_PLANAR_LENGTH_M),
        TRANSECT_TOTAL_SUBSTRATE_COVER_M =
          mean(TRANSECT_TOTAL_SUBSTRATE_COVER_M),
        SUBSTRATE_COVER_CM =
          sum(SUBSTRATE_COVER_CM, na.rm = TRUE),
        SUBSTRATE_COVER_PCT =
          SUBSTRATE_COVER_CM /
          TRANSECT_TOTAL_SUBSTRATE_COVER_M,
        SUBSTRATE_PLANAR_PROD_KG_YR =
          sum(COLONY_PROD_G_YR, na.rm = TRUE) / 1000,
        SUBSTRATE_PLANAR_PROD_KG_YR_L95 =
          sum(COLONY_PROD_G_YR_L95, na.rm = TRUE) / 1000,
        SUBSTRATE_PLANAR_PROD_KG_YR_U95 =
          sum(COLONY_PROD_G_YR_U95, na.rm = TRUE) / 1000,
        SUBSTRATE_CARB_PROD_KG_M2_YR =
          SUBSTRATE_PLANAR_PROD_KG_YR *
          (10000 / (TRANSECT_PLANAR_LENGTH_M * 100)),
        SUBSTRATE_CARB_PROD_KG_M2_YR_L95 =
          SUBSTRATE_PLANAR_PROD_KG_YR_L95 *
          (10000 / (TRANSECT_PLANAR_LENGTH_M * 100)),
        SUBSTRATE_CARB_PROD_KG_M2_YR_U95 =
          SUBSTRATE_PLANAR_PROD_KG_YR_U95 *
          (10000 / (TRANSECT_PLANAR_LENGTH_M * 100))
      )
  )
  

  substrate_code_morphology_all <-
    unique(data[c("SUBSTRATE_CODE", "MORPHOLOGYCODE")])

  substrate_code_morphology_all$SUBSTRATE_CODE_MORPHOLOGYCODE <-
    paste0(substrate_code_morphology_all$SUBSTRATE_CODE, "-", substrate_code_morphology_all$MORPHOLOGYCODE)

  substrate_code_full_table <- expand.grid(
    OCC_SITEID = unique(transect_summary$OCC_SITEID),
    CB_TRANSECTID = unique(transect_summary$CB_TRANSECTID),
    SUBSTRATE_CODE_MORPHOLOGYCODE = substrate_code_morphology_all$SUBSTRATE_CODE_MORPHOLOGYCODE
  )
  
  substrate_code_full_table$SUBSTRATE_CODE <- substrate_code_morphology_all$SUBSTRATE_CODE[
    match(substrate_code_full_table$SUBSTRATE_CODE_MORPHOLOGYCODE,
          substrate_code_morphology_all$SUBSTRATE_CODE_MORPHOLOGYCODE)]

  substrate_code_full_table$MORPHOLOGYCODE <- substrate_code_morphology_all$MORPHOLOGYCODE[
    match(substrate_code_full_table$SUBSTRATE_CODE_MORPHOLOGYCODE,
          substrate_code_morphology_all$SUBSTRATE_CODE_MORPHOLOGYCODE)]

  substrate_code_full_table$OCC_SITEID_TRANSECT <-
    paste(substrate_code_full_table$OCC_SITEID,
          substrate_code_full_table$CB_TRANSECTID,
          sep = "-")
  
   substrate_code_full_table <- substrate_code_full_table[substrate_code_full_table$OCC_SITEID_TRANSECT %in% transect_summary$OCC_SITEID_TRANSECT ,]

  # Populate full table with transect-level data
  summary_transect_substratecode <-  merge(
    substrate_code_full_table,
    transect_substratecode,
    by = c(
      "OCC_SITEID",
      "CB_TRANSECTID",
      "OCC_SITEID_TRANSECT",
      "SUBSTRATE_CODE",
      "MORPHOLOGYCODE"
    ),
    all.x = TRUE
  )

  prod_dbase$SUBSTRATE_CODE_MORPHOLOGYCODE <-
    paste0(prod_dbase$SUBSTRATE_CODE, "-", prod_dbase$MORPHOLOGYCODE)

  summary_transect_substratecode$MORPHOLOGY <- prod_dbase$MORPHOLOGY[match(summary_transect_substratecode$SUBSTRATE_CODE_MORPHOLOGYCODE, as.factor(prod_dbase$SUBSTRATE_CODE_MORPHOLOGYCODE))]

  summary_transect_substratecode$MORPHOLOGY[summary_transect_substratecode$SUBSTRATE_CODE == "PRUS"] <- "Laminar Columnar"
  summary_transect_substratecode$MORPHOLOGY[summary_transect_substratecode$SUBSTRATE_CODE == "PMRC"] <- "Laminar Columnar"
  
  summary_transect_substratecode$SUBSTRATE_CLASS <- prod_dbase$SUBSTRATE_CLASS[match(summary_transect_substratecode$SUBSTRATE_CODE_MORPHOLOGYCODE, as.factor(prod_dbase$SUBSTRATE_CODE_MORPHOLOGYCODE))]

  summary_transect_substratecode$SUBSTRATE_CLASS[summary_transect_substratecode$SUBSTRATE_CODE == "PRUS"] <- "CORAL"
  summary_transect_substratecode$SUBSTRATE_CLASS[summary_transect_substratecode$SUBSTRATE_CODE == "PMRC"] <- "CORAL"

  summary_transect_substratecode$SUBSTRATE_NAME <- prod_dbase$SUBSTRATE_NAME[match(summary_transect_substratecode$SUBSTRATE_CODE_MORPHOLOGYCODE, as.factor(prod_dbase$SUBSTRATE_CODE_MORPHOLOGYCODE))]

  summary_transect_substratecode$SUBSTRATE_NAME[summary_transect_substratecode$SUBSTRATE_CODE == "PRUS"] <- "Porites rus"
  summary_transect_substratecode$SUBSTRATE_NAME[summary_transect_substratecode$SUBSTRATE_CODE == "PMRC"] <- "Porites monticulosa/rus complex"
  

  
  summary_transect_substratecode <- summary_transect_substratecode[c(
    "REGION",
    "REGIONCODE",
    "YEAR",
    "CRUISE_ID",
    "LOCATION",
    "LOCATIONCODE",
    "OCC_SITEID",
    "SITEVISITID",
    "LATITUDE",
    "LONGITUDE",
    "SITE_DEPTH_M",
    "LOCALDATE",
    "CB_TRANSECTID",
    "OCC_SITEID_TRANSECT",
    "SUBSTRATE_CLASS",
    "SUBSTRATE_NAME",
    "SUBSTRATE_CODE",
    "MORPHOLOGY",
    "MORPHOLOGYCODE",
    "TRANSECT_PLANAR_LENGTH_M",
    "TRANSECT_TOTAL_SUBSTRATE_COVER_M",
    "SUBSTRATE_COVER_CM",
    "SUBSTRATE_COVER_PCT",
    "SUBSTRATE_PLANAR_PROD_KG_YR",
    "SUBSTRATE_PLANAR_PROD_KG_YR_L95",
    "SUBSTRATE_PLANAR_PROD_KG_YR_U95",
    "SUBSTRATE_CARB_PROD_KG_M2_YR",
    "SUBSTRATE_CARB_PROD_KG_M2_YR_L95",
    "SUBSTRATE_CARB_PROD_KG_M2_YR_U95"
  )]


  summary_transect_substratecode <-
    summary_transect_substratecode %>%
    group_by(OCC_SITEID) %>%
    fill(
      c(
        "REGION",
        "REGIONCODE",
        "YEAR",
        "CRUISE_ID",
        "LOCATION",
        "LOCATIONCODE",
        "OCC_SITEID",
        "SITEVISITID",
        "LATITUDE",
        "LONGITUDE",
        "SITE_DEPTH_M",
        "LOCALDATE",
        "TRANSECT_PLANAR_LENGTH_M",
        "TRANSECT_TOTAL_SUBSTRATE_COVER_M",
      ),
      .direction = 'downup'
    )
  
  summary_transect_substratecode[, c((which(
    colnames(summary_transect_substratecode) == "SUBSTRATE_COVER_CM"
  )):(which(
    colnames(summary_transect_substratecode) == "SUBSTRATE_CARB_PROD_KG_M2_YR_U95"
  )))][is.na(summary_transect_substratecode[, c((which(
    colnames(summary_transect_substratecode) == "SUBSTRATE_COVER_CM"
  )):(which(
    colnames(summary_transect_substratecode) == "SUBSTRATE_CARB_PROD_KG_M2_YR_U95"
  )))])] <- 0
  
  # Summarize production by SUBSTRATE_CLASS at TRANSECT level  ---------------
  transect_substrateclass <- suppressMessages(
    data %>% dplyr::group_by(
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      SITEVISITID,
      LATITUDE,
      LONGITUDE,
      SITE_DEPTH_M,
      LOCALDATE,
      CB_TRANSECTID,
      OCC_SITEID_TRANSECT,
      SUBSTRATE_CLASS
    ) %>%
      reframe(
        TRANSECT_PLANAR_LENGTH_M = first(TRANSECT_PLANAR_LENGTH_M),
        TRANSECT_TOTAL_SUBSTRATE_COVER_M =
          first(TRANSECT_TOTAL_SUBSTRATE_COVER_M),
        SUBSTRATE_COVER_CM =
          sum(SUBSTRATE_COVER_CM, na.rm = TRUE),
        SUBSTRATE_COVER_PCT =
          SUBSTRATE_COVER_CM / TRANSECT_TOTAL_SUBSTRATE_COVER_M,
        SUBSTRATE_PLANAR_PROD_KG_YR =
          sum(COLONY_PROD_G_YR, na.rm = TRUE) / 1000,
        SUBSTRATE_PLANAR_PROD_KG_YR_L95 =
          sum(COLONY_PROD_G_YR_L95, na.rm = TRUE) / 1000,
        SUBSTRATE_PLANAR_PROD_KG_YR_U95 =
          sum(COLONY_PROD_G_YR_U95, na.rm = TRUE) / 1000,
        SUBSTRATE_CARB_PROD_KG_M2_YR =
          SUBSTRATE_PLANAR_PROD_KG_YR *
          (10000 / (TRANSECT_PLANAR_LENGTH_M * 100)),
        SUBSTRATE_CARB_PROD_KG_M2_YR_L95 =
          SUBSTRATE_PLANAR_PROD_KG_YR_L95 *
          (10000 / (TRANSECT_PLANAR_LENGTH_M * 100)),
        SUBSTRATE_CARB_PROD_KG_M2_YR_U95 =
          SUBSTRATE_PLANAR_PROD_KG_YR_U95 *
          (10000 / (TRANSECT_PLANAR_LENGTH_M * 100))
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
  
 substrate_full_table <- substrate_full_table[substrate_full_table$OCC_SITEID_TRANSECT %in% transect_summary$OCC_SITEID_TRANSECT ,]

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
    "SITEVISITID",
    "LATITUDE",
    "LONGITUDE",
    "SITE_DEPTH_M",
    "LOCALDATE",
    "CB_TRANSECTID",
    "OCC_SITEID_TRANSECT",
    "SUBSTRATE_CLASS",
    "TRANSECT_PLANAR_LENGTH_M",
    "TRANSECT_TOTAL_SUBSTRATE_COVER_M",
    "SUBSTRATE_COVER_CM",
    "SUBSTRATE_COVER_PCT",
    "SUBSTRATE_PLANAR_PROD_KG_YR",
    "SUBSTRATE_PLANAR_PROD_KG_YR_L95",
    "SUBSTRATE_PLANAR_PROD_KG_YR_U95",
    "SUBSTRATE_CARB_PROD_KG_M2_YR",
    "SUBSTRATE_CARB_PROD_KG_M2_YR_L95",
    "SUBSTRATE_CARB_PROD_KG_M2_YR_U95"
  )]


  summary_transect_substrateclass <-
    summary_transect_substrateclass %>%
    group_by(OCC_SITEID) %>%
    fill(
      c(
        "REGION",
        "REGIONCODE",
        "YEAR",
        "CRUISE_ID",
        "LOCATION",
        "LOCATIONCODE",
        "OCC_SITEID",
        "SITEVISITID",
        "LATITUDE",
        "LONGITUDE",
        "SITE_DEPTH_M",
        "LOCALDATE",
        "TRANSECT_PLANAR_LENGTH_M",
        "TRANSECT_TOTAL_SUBSTRATE_COVER_M",
      ),
      .direction = 'downup'
    )

  summary_transect_substrateclass[, c((which(
    colnames(summary_transect_substrateclass) == "SUBSTRATE_COVER_CM"
  )):(which(
    colnames(summary_transect_substrateclass) == "SUBSTRATE_CARB_PROD_KG_M2_YR_U95"
  )))][is.na(summary_transect_substrateclass[, c((which(
    colnames(summary_transect_substrateclass) == "SUBSTRATE_COVER_CM"
  )):(which(
    colnames(summary_transect_substrateclass) == "SUBSTRATE_CARB_PROD_KG_M2_YR_U95"
  )))])] <- 0

  # Summarize production by CORAL_GROUP at TRANSECT level  -------------------
  transect_coral_group <- data %>% dplyr::group_by(
    REGION,
    REGIONCODE,
    YEAR,
    CRUISE_ID,
    LOCATION,
    LOCATIONCODE,
    OCC_SITEID,
    SITEVISITID,
    LATITUDE,
    LONGITUDE,
    SITE_DEPTH_M,
    LOCALDATE,
    CB_TRANSECTID,
    OCC_SITEID_TRANSECT,
    CORAL_GROUP,
    CORAL_GROUP_NAME,
  ) %>%
    reframe(
      TRANSECT_PLANAR_LENGTH_M =
        first(TRANSECT_PLANAR_LENGTH_M),
      TRANSECT_TOTAL_SUBSTRATE_COVER_M =
        first(TRANSECT_TOTAL_SUBSTRATE_COVER_M),
      SUBSTRATE_COVER_CM =
        sum(SUBSTRATE_COVER_CM, na.rm = TRUE),
      SUBSTRATE_COVER_PCT =
        SUBSTRATE_COVER_CM / TRANSECT_TOTAL_SUBSTRATE_COVER_M,
      SUBSTRATE_PLANAR_PROD_KG_YR =
        sum(COLONY_PROD_G_YR, na.rm = TRUE) / 1000,
      SUBSTRATE_PLANAR_PROD_KG_YR_L95 =
        sum(COLONY_PROD_G_YR_L95, na.rm = TRUE) / 1000,
      SUBSTRATE_PLANAR_PROD_KG_YR_U95 =
        sum(COLONY_PROD_G_YR_U95, na.rm = TRUE) / 1000,
      SUBSTRATE_CARB_PROD_KG_M2_YR =
        SUBSTRATE_PLANAR_PROD_KG_YR *
        (10000 / (TRANSECT_PLANAR_LENGTH_M * 100)),
      SUBSTRATE_CARB_PROD_KG_M2_YR_L95 =
        SUBSTRATE_PLANAR_PROD_KG_YR_L95 *
        (10000 / (TRANSECT_PLANAR_LENGTH_M * 100)),
      SUBSTRATE_CARB_PROD_KG_M2_YR_U95 =
        SUBSTRATE_PLANAR_PROD_KG_YR_U95 *
        (10000 / (TRANSECT_PLANAR_LENGTH_M * 100))
    )

  # transect_coral_group <-
  #   subset(transect_coral_group,
  #          transect_coral_group$CORAL_GROUP != "NA")

  coral_group_all <- c(
    "ACBR",
    "ACTA",
    "ASMA",
    "BR",
    "ENC",
    "FOL",
    "FREE",
    "GOMA",
    "MASS",
    "MOBR",
    "MOEN",
    "MOFO",
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
  
  coral_full_table <- coral_full_table[coral_full_table$OCC_SITEID_TRANSECT %in% transect_summary$OCC_SITEID_TRANSECT ,]

  # Populate full table with transect-level data
  
  if (all(is.na(transect_coral_group$CORAL_GROUP)) == TRUE) {
    
    transect_coral_group <- transect_coral_group %>% select(-c(CORAL_GROUP, CORAL_GROUP_NAME))
    
    summary_transect_coral <- merge(
      coral_full_table,
      transect_coral_group,
      by = c("OCC_SITEID", "CB_TRANSECTID", "OCC_SITEID_TRANSECT"),
      all.x = TRUE
      )
    
      summary_transect_coral$CORAL_GROUP_NAME <- prod_dbase$CORAL_GROUP_NAME[match(summary_transect_coral$CORAL_GROUP, prod_dbase$CORAL_GROUP)] 
      
    
  } else {
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
  } 
  
  summary_transect_coral <- summary_transect_coral[c(
    "REGION",
    "REGIONCODE",
    "YEAR",
    "CRUISE_ID",
    "LOCATION",
    "LOCATIONCODE",
    "OCC_SITEID",
    "SITEVISITID",
    "LATITUDE",
    "LONGITUDE",
    "SITE_DEPTH_M",
    "LOCALDATE",
    "CB_TRANSECTID",
    "OCC_SITEID_TRANSECT",
    "CORAL_GROUP",
    "CORAL_GROUP_NAME",
    "TRANSECT_PLANAR_LENGTH_M",
    "TRANSECT_TOTAL_SUBSTRATE_COVER_M",
    "SUBSTRATE_COVER_CM",
    "SUBSTRATE_COVER_PCT",
    "SUBSTRATE_PLANAR_PROD_KG_YR",
    "SUBSTRATE_PLANAR_PROD_KG_YR_L95",
    "SUBSTRATE_PLANAR_PROD_KG_YR_U95",
    "SUBSTRATE_CARB_PROD_KG_M2_YR",
    "SUBSTRATE_CARB_PROD_KG_M2_YR_L95",
    "SUBSTRATE_CARB_PROD_KG_M2_YR_U95"
  )]


  summary_transect_coral <-
    summary_transect_coral %>% group_by(OCC_SITEID) %>%
    fill(
      c(
        "REGION",
        "REGIONCODE",
        "YEAR",
        "CRUISE_ID",
        "LOCATION",
        "LOCATIONCODE",
        "OCC_SITEID",
        "SITEVISITID",
        "LATITUDE",
        "LONGITUDE",
        "SITE_DEPTH_M",
        "LOCALDATE",
        "TRANSECT_PLANAR_LENGTH_M",
        "TRANSECT_TOTAL_SUBSTRATE_COVER_M",
      ),
      .direction = "downup"
    )
  
  
  for (i in 1:nrow(summary_transect_coral))
  {
    if (is.na(unique(summary_transect_coral$LATITUDE[i])) == TRUE) {
      summary_transect_coral$REGION[i]  <- transect_summary$REGION[match(summary_transect_coral$OCC_SITEID[i],
                                                                      transect_summary$OCC_SITEID[i])]
      summary_transect_coral$REGIONCODE[i]  <- transect_summary$REGIONCODE[match(summary_transect_coral$OCC_SITEID[i],
                                                                              transect_summary$OCC_SITEID[i])]
      summary_transect_coral$YEAR[i]  <- transect_summary$YEAR[match(summary_transect_coral$OCC_SITEID[i],
                                                                  transect_summary$OCC_SITEID[i])]
      summary_transect_coral$CRUISE_ID[i]  <- transect_summary$CRUISE_ID[match(summary_transect_coral$OCC_SITEID[i],
                                                                            transect_summary$OCC_SITEID[i])]
      summary_transect_coral$LOCATIONCODE[i]  <- transect_summary$LOCATIONCODE[match(summary_transect_coral$OCC_SITEID[i],
                                                                                  transect_summary$OCC_SITEID[i])]
      summary_transect_coral$LOCATION[i]  <- transect_summary$LOCATION[match(summary_transect_coral$OCC_SITEID[i],
                                                                          transect_summary$OCC_SITEID[i])]
      summary_transect_coral$LATITUDE[i]  <- transect_summary$LATITUDE[match(summary_transect_coral$OCC_SITEID[i],
                                                                          transect_summary$OCC_SITEID[i])]
      summary_transect_coral$LONGITUDE[i]  <- transect_summary$LONGITUDE[match(summary_transect_coral$OCC_SITEID[i],
                                                                            transect_summary$OCC_SITEID[i])]
      summary_transect_coral$SITE_DEPTH_M[i]  <- transect_summary$SITE_DEPTH_M[match(summary_transect_coral$OCC_SITEID[i],
                                                                                  transect_summary$OCC_SITEID[i])]
      summary_transect_coral$LOCALDATE[i]  <- transect_summary$LOCALDATE[match(summary_transect_coral$OCC_SITEID[i],
                                                                            transect_summary$OCC_SITEID[i])]
      summary_transect_coral$TRANSECT_PLANAR_LENGTH_M[i]  <- transect_summary$TRANSECT_PLANAR_LENGTH_M[match(summary_transect_coral$OCC_SITEID[i],
                                                                                                          transect_summary$OCC_SITEID[i])]
      summary_transect_coral$TRANSECT_TOTAL_SUBSTRATE_COVER_M[i]  <- transect_summary$TRANSECT_TOTAL_SUBSTRATE_COVER_M[match(summary_transect_coral$OCC_SITEID[i],
                                                                                                                          transect_summary$OCC_SITEID[i])]
    }
  }
  
  summary_transect_coral <-
    summary_transect_coral %>% group_by(CORAL_GROUP) %>%
    fill(
      c(
        "CORAL_GROUP_NAME"
      ),
      .direction = "downup"
    )

  summary_transect_coral$CORAL_GROUP_NAME <- prod_dbase$CORAL_GROUP_NAME[match(summary_transect_coral$CORAL_GROUP, as.factor(prod_dbase$CORAL_GROUP))]

  summary_transect_coral[, c((which(
    colnames(summary_transect_coral) == "SUBSTRATE_COVER_CM"
  )):(which(
    colnames(summary_transect_coral) == "SUBSTRATE_CARB_PROD_KG_M2_YR_U95"
  )))][is.na(summary_transect_coral[, c((which(
    colnames(summary_transect_coral) == "SUBSTRATE_COVER_CM"
  )):(which(
    colnames(summary_transect_coral) == "SUBSTRATE_CARB_PROD_KG_M2_YR_U95"
  )))])] <- 0
  
  # Summarize production by TRANSECT ----------------------------------------
  if (dbase_type == "IPRB") {

    summary_transect <- summary_transect_substratecode  %>%
      dplyr::group_by(
        REGION,
        REGIONCODE,
        YEAR,
        CRUISE_ID,
        LOCATION,
        LOCATIONCODE,
        OCC_SITEID,
        SITEVISITID,
        LATITUDE,
        LONGITUDE,
        SITE_DEPTH_M,
        LOCALDATE,
        CB_TRANSECTID,
        OCC_SITEID_TRANSECT
      ) %>%
      reframe(
        TRANSECT_PLANAR_LENGTH_M =
          mean(TRANSECT_PLANAR_LENGTH_M),
        TRANSECT_TOTAL_SUBSTRATE_COVER_M =
          mean(TRANSECT_TOTAL_SUBSTRATE_COVER_M),
        GROSS_CARB_PROD_KG_M2_YR =
          sum(SUBSTRATE_CARB_PROD_KG_M2_YR, na.rm = TRUE),
        GROSS_CARB_PROD_KG_M2_YR_L95 =
          sum(SUBSTRATE_CARB_PROD_KG_M2_YR_L95, na.rm = TRUE),
        GROSS_CARB_PROD_KG_M2_YR_U95 =
          sum(SUBSTRATE_CARB_PROD_KG_M2_YR_U95, na.rm = TRUE),
        RUGOSITY =
          TRANSECT_TOTAL_SUBSTRATE_COVER_M /
          TRANSECT_PLANAR_LENGTH_M,
        SUBSTRATE_AVAILABLE_MACRO_PCT =
          sum(SUBSTRATE_COVER_PCT[SUBSTRATE_CODE %in% c("ART",
                                                        "BOR",
                                                        "CCA",
                                                        "CYA",
                                                        "DC",
                                                        "HA",
                                                        "LSP",
                                                        "MAC",
                                                        "MCA",
                                                        "OCE",
                                                        "OTH",
                                                        "RUB",
                                                        "RUBT",
                                                        "RUBC",
                                                        "SCA",
                                                        "TF"
                                                        )]),
        SUBSTRATE_AVAILABLE_MACRO_INDEX =
          RUGOSITY * (SUBSTRATE_AVAILABLE_MACRO_PCT / 100),
        SUBSTRATE_AVAILABLE_MICRO_PCT =
          sum(SUBSTRATE_COVER_PCT[SUBSTRATE_CODE %in% c("ART",
                                                        "CYA",
                                                        "DC",
                                                        "HA",
                                                        "LSP",
                                                        "MAC",
                                                        "RUB",
                                                        "RUBT",
                                                        "TF"
                                                        )]),
        SUBSTRATE_AVAILABLE_MICRO_INDEX =
          RUGOSITY * (SUBSTRATE_AVAILABLE_MICRO_PCT / 100),
        HARD_CORAL_COVER_PCT =
          sum(SUBSTRATE_COVER_PCT[SUBSTRATE_CLASS == "CORAL"]),
        CCA_COVER_PCT =
          sum(SUBSTRATE_COVER_PCT[SUBSTRATE_CLASS == "CCA"]),
        HARD_CORAL_CARB_PROD_KG_M2_YR =
          sum(SUBSTRATE_CARB_PROD_KG_M2_YR[
            SUBSTRATE_CLASS == "CORAL"],
              na.rm = TRUE),
        CCA_CARB_PROD_KG_M2_YR =
          sum(SUBSTRATE_CARB_PROD_KG_M2_YR[
            SUBSTRATE_CLASS == "CCA"],
              na.rm = TRUE),
        MACROBIOEROSION_KG_M2_YR =
          SUBSTRATE_AVAILABLE_MACRO_INDEX * macro_rate,
        MACROBIOEROSION_KG_M2_YR_L95 =
          SUBSTRATE_AVAILABLE_MACRO_INDEX * (macro_rate - macro_rate_ci),
        MACROBIOEROSION_KG_M2_YR_U95 =
          SUBSTRATE_AVAILABLE_MACRO_INDEX * (macro_rate + macro_rate_ci),
        MICROBIOEROSION_KG_M2_YR =
          SUBSTRATE_AVAILABLE_MICRO_INDEX * micro_rate,
        MICROBIOEROSION_KG_M2_YR_L95 =
          SUBSTRATE_AVAILABLE_MICRO_INDEX * (micro_rate - micro_rate_ci),
        MICROBIOEROSION_KG_M2_YR_U95 =
          SUBSTRATE_AVAILABLE_MICRO_INDEX * (micro_rate + micro_rate_ci),
        BIOEROSION_KG_M2_YR =
          MACROBIOEROSION_KG_M2_YR +
          MICROBIOEROSION_KG_M2_YR,
        BIOEROSION_KG_M2_YR_L95 =
          MACROBIOEROSION_KG_M2_YR_L95 +
          MICROBIOEROSION_KG_M2_YR_L95,
        BIOEROSION_KG_M2_YR_U95 =
          MACROBIOEROSION_KG_M2_YR_U95 +
          MICROBIOEROSION_KG_M2_YR_U95
      )
  }
  if (dbase_type == "NCRMP") {
    summary_transect <-
      summary_transect_substratecode  %>% dplyr::group_by(
        REGION,
        REGIONCODE,
        YEAR,
        CRUISE_ID,
        LOCATION,
        LOCATIONCODE,
        OCC_SITEID,
        SITEVISITID,
        LATITUDE,
        LONGITUDE,
        SITE_DEPTH_M,
        LOCALDATE,
        CB_TRANSECTID,
        OCC_SITEID_TRANSECT,
      ) %>%
      reframe(
        TRANSECT_PLANAR_LENGTH_M =
          mean(TRANSECT_PLANAR_LENGTH_M),
        TRANSECT_TOTAL_SUBSTRATE_COVER_M =
          mean(TRANSECT_TOTAL_SUBSTRATE_COVER_M),
        GROSS_CARB_PROD_KG_M2_YR =
          sum(SUBSTRATE_CARB_PROD_KG_M2_YR, na.rm = TRUE),
        GROSS_CARB_PROD_KG_M2_YR_L95 =
          sum(SUBSTRATE_CARB_PROD_KG_M2_YR_L95, na.rm = TRUE),
        GROSS_CARB_PROD_KG_M2_YR_U95 =
          sum(SUBSTRATE_CARB_PROD_KG_M2_YR_U95, na.rm = TRUE),
        RUGOSITY =
          TRANSECT_TOTAL_SUBSTRATE_COVER_M /
          TRANSECT_PLANAR_LENGTH_M,
        SUBSTRATE_AVAILABLE_MACRO_PCT =
          sum(SUBSTRATE_COVER_PCT[SUBSTRATE_CODE %in% c("HARD",
                                                        "DC",
                                                        "RUB",
                                                        "TURF",
                                                        "CYANO",
                                                        "CCA",
                                                        "MA",
                                                        "PESP",
                                                        "HALI",
                                                        "SP"
                                                        )]),
        SUBSTRATE_AVAILABLE_MACRO_INDEX =
          RUGOSITY * (SUBSTRATE_AVAILABLE_MACRO_PCT / 100),
        SUBSTRATE_AVAILABLE_MICRO_PCT =
          sum(SUBSTRATE_COVER_PCT[SUBSTRATE_CODE %in% c("HARD",
                                                        "DC",
                                                        "RUB",
                                                        "TURF",
                                                        "CYANO",
                                                        "MA",
                                                        "HALI"
                                                        )]),
        SUBSTRATE_AVAILABLE_MICRO_INDEX =
          RUGOSITY * (SUBSTRATE_AVAILABLE_MICRO_PCT / 100),
        HARD_CORAL_COVER_PCT =
          sum(SUBSTRATE_COVER_PCT[SUBSTRATE_CLASS == "CORAL"]),
        CCA_COVER_PCT =
          sum(SUBSTRATE_COVER_PCT[SUBSTRATE_CLASS == "CCA"]),
        HARD_CORAL_CARB_PROD_KG_M2_YR =
          sum(SUBSTRATE_CARB_PROD_KG_M2_YR[
            SUBSTRATE_CLASS == "CORAL"],
              na.rm = TRUE),
        CCA_CARB_PROD_KG_M2_YR =
          sum(SUBSTRATE_CARB_PROD_KG_M2_YR[
            SUBSTRATE_CLASS == "CCA"],
              na.rm = TRUE),
        MACROBIOEROSION_KG_M2_YR =
          SUBSTRATE_AVAILABLE_MACRO_INDEX * macro_rate,
        MACROBIOEROSION_KG_M2_YR_L95 =
          SUBSTRATE_AVAILABLE_MACRO_INDEX * (macro_rate - macro_rate_ci),
        MACROBIOEROSION_KG_M2_YR_U95 =
          SUBSTRATE_AVAILABLE_MACRO_INDEX * (macro_rate + macro_rate_ci),
        MICROBIOEROSION_KG_M2_YR =
          SUBSTRATE_AVAILABLE_MICRO_INDEX * micro_rate,
        MICROBIOEROSION_KG_M2_YR_L95 =
          SUBSTRATE_AVAILABLE_MICRO_INDEX * (micro_rate - micro_rate_ci),
        MICROBIOEROSION_KG_M2_YR_U95 =
          SUBSTRATE_AVAILABLE_MICRO_INDEX * (micro_rate + micro_rate_ci),
        BIOEROSION_KG_M2_YR =
          MACROBIOEROSION_KG_M2_YR +
          MICROBIOEROSION_KG_M2_YR,
        BIOEROSION_KG_M2_YR_L95 =
          MACROBIOEROSION_KG_M2_YR_L95 +
          MICROBIOEROSION_KG_M2_YR_L95,
        BIOEROSION_KG_M2_YR_U95 =
          MACROBIOEROSION_KG_M2_YR_U95 +
          MICROBIOEROSION_KG_M2_YR_U95
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
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      SITEVISITID,
      LATITUDE,
      LONGITUDE,
      SITE_DEPTH_M,
      LOCALDATE,
      SUBSTRATE_CLASS,
      SUBSTRATE_CODE,
      MORPHOLOGY,
      MORPHOLOGYCODE
    ) %>%
    dplyr::reframe(across(
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
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      SITEVISITID,
      LATITUDE,
      LONGITUDE,
      SITE_DEPTH_M,
      LOCALDATE,
      SUBSTRATE_CLASS
    ) %>%
  dplyr::reframe(across(
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
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      SITEVISITID,
      LATITUDE,
      LONGITUDE,
      SITE_DEPTH_M,
      LOCALDATE,
      CORAL_GROUP,
      CORAL_GROUP_NAME
    ) %>%
    dplyr::reframe(across(
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
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      SITEVISITID,
      LATITUDE,
      LONGITUDE,
      SITE_DEPTH_M,
      LOCALDATE
    ) %>%
    dplyr::reframe(across(
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


