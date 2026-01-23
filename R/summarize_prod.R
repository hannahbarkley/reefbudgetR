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
                           macro_rates = c("IPRB", "NCRMP"),
                           micro_rates = c("IPRB", "NCRMP"),
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

  transect_substratecode <-
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
          first(TRANSECT_PLANAR_LENGTH_M),
        TRANSECT_TOTAL_SUBSTRATE_COVER_M =
          first(TRANSECT_TOTAL_SUBSTRATE_COVER_M),
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
  
  # Create unique pairs and lookup ID
  substrate_code_morphology_all <- data %>%
    distinct(SUBSTRATE_CODE, MORPHOLOGYCODE) %>%
    unite(
      "SUBSTRATE_CODE_MORPHOLOGYCODE", 
      SUBSTRATE_CODE, MORPHOLOGYCODE, 
      sep = "-", 
      remove = FALSE, # Keeps original columns
      na.rm = TRUE    # If MORPHOLOGYCODE is NA, it won't add "-NA"
    )

  # Build the full substrate code table for each OCC SITE and substrate/morphology pair
  
  substrate_code_full_table <- expand_grid(
    OCC_SITEID = unique(transect_summary$OCC_SITEID),
    CB_TRANSECTID = unique(transect_summary$CB_TRANSECTID),
    substrate_code_morphology_all %>%
      select(
        SUBSTRATE_CODE,
        MORPHOLOGYCODE,
        SUBSTRATE_CODE_MORPHOLOGYCODE
      )
  ) %>%
    mutate(OCC_SITEID_TRANSECT = str_c(OCC_SITEID, CB_TRANSECTID, sep = "-"))
  
  substrate_code_full_table <- substrate_code_full_table[substrate_code_full_table$OCC_SITEID_TRANSECT %in% transect_summary$OCC_SITEID_TRANSECT , ]
  
  
  # Populate full table with transect-level data
  
  summary_transect_substratecode <- substrate_code_full_table %>%
    left_join(
      transect_substratecode, 
      by = c("OCC_SITEID", "CB_TRANSECTID", "OCC_SITEID_TRANSECT", "SUBSTRATE_CODE", "MORPHOLOGYCODE")
    ) 
  
  summary_transect_substratecode <- summary_transect_substratecode %>%
    left_join(
      prod_dbase %>% 
        select(SUBSTRATE_CODE, MORPHOLOGYCODE, MORPHOLOGY, SUBSTRATE_CLASS, SUBSTRATE_NAME) %>%
        distinct(), 
      by = c("SUBSTRATE_CODE", "MORPHOLOGYCODE"),
      suffix = c("", ".db") # Adds .db to columns from prod_dbase if names overlap
    ) %>%
    mutate(
      SUBSTRATE_CLASS = coalesce(SUBSTRATE_CLASS, SUBSTRATE_CLASS.db),
      SUBSTRATE_NAME = coalesce(SUBSTRATE_NAME, SUBSTRATE_NAME.db),
      MORPHOLOGY = coalesce(MORPHOLOGY,MORPHOLOGY.db),
    ) %>%
    select(-ends_with(".db")) %>%
    mutate(
      MORPHOLOGY = case_when(
        SUBSTRATE_CODE %in% c("PRUS", "PMRC") ~ "Laminar Columnar",
        TRUE ~ MORPHOLOGY # Keeps the joined value if not PRUS or PMRC
      ),
      SUBSTRATE_CLASS = case_when(
        SUBSTRATE_CODE %in% c("PRUS", "PMRC") ~ "CORAL",
        TRUE ~ SUBSTRATE_CLASS
      ),
      SUBSTRATE_NAME = case_when(
        SUBSTRATE_CODE == "PRUS" ~ "Porites rus",
        SUBSTRATE_CODE == "PMRC" ~ "Porites monticulosa/rus complex",
        TRUE ~ SUBSTRATE_NAME
      )
    ) %>%
    group_by(OCC_SITEID) %>%
    fill(
      REGION, REGIONCODE, YEAR, CRUISE_ID, LOCATION, LOCATIONCODE, 
      SITEVISITID, LATITUDE, LONGITUDE, SITE_DEPTH_M, LOCALDATE, 
      TRANSECT_PLANAR_LENGTH_M, TRANSECT_TOTAL_SUBSTRATE_COVER_M,
      .direction = 'downup'
    ) %>%
    ungroup() %>%
    mutate(across(
      c(SUBSTRATE_COVER_CM:SUBSTRATE_CARB_PROD_KG_M2_YR_U95),
      ~ replace_na(.x, 0)
    ))
  
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

  
  # Summarize production by SUBSTRATE_CLASS at TRANSECT level  ---------------
  transect_substrateclass <- 
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
  

  # Create a table of all substrate classes for each transect
  substrate_class_all <- c("CORAL",
                           "CCA",
                           "CARB",
                           "MA",
                           "RCK",
                           "SAND",
                           "TURF",
                           "OTHER")
  

  substrate_full_table <- expand_grid(
    OCC_SITEID = unique(transect_summary$OCC_SITEID),
    CB_TRANSECTID = unique(transect_summary$CB_TRANSECTID),
    SUBSTRATE_CLASS = substrate_class_all
  ) %>%
    mutate(OCC_SITEID_TRANSECT = str_c(OCC_SITEID, CB_TRANSECTID, sep = "-"))
  
  
  substrate_full_table <- substrate_full_table[substrate_full_table$OCC_SITEID_TRANSECT %in% transect_summary$OCC_SITEID_TRANSECT , ]

  # Populate full table with transect-level data
  summary_transect_substrateclass <-  left_join(
    substrate_full_table,
    transect_substrateclass,
    by = c(
      "OCC_SITEID",
      "CB_TRANSECTID",
      "OCC_SITEID_TRANSECT",
      "SUBSTRATE_CLASS"
    )
  )
  
  summary_transect_substrateclass <- summary_transect_substrateclass %>%
    group_by(OCC_SITEID) %>%
    fill(
      REGION:TRANSECT_TOTAL_SUBSTRATE_COVER_M, 
      .direction = 'downup'
    ) %>%
    ungroup() %>%
    mutate(across(
      SUBSTRATE_COVER_CM:SUBSTRATE_CARB_PROD_KG_M2_YR_U95, 
      ~ replace_na(.x, 0)
    ))

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
  ) %>%
    mutate(OCC_SITEID_TRANSECT = str_c(OCC_SITEID, CB_TRANSECTID, sep = "-")) %>%
    left_join(
      prod_dbase %>% select(CORAL_GROUP, CORAL_GROUP_NAME) %>% distinct (CORAL_GROUP, CORAL_GROUP_NAME), by = "CORAL_GROUP"
    )
  
  coral_full_table <- coral_full_table[coral_full_table$OCC_SITEID_TRANSECT %in% transect_summary$OCC_SITEID_TRANSECT ,]

  
  # Populate full table with transect-level data
  
  if (all(is.na(transect_coral_group$CORAL_GROUP)) == TRUE) {
    
    transect_coral_group <- transect_coral_group %>% select(-c(CORAL_GROUP, CORAL_GROUP_NAME)) %>%
      mutate(across(
        SUBSTRATE_COVER_CM:SUBSTRATE_CARB_PROD_KG_M2_YR_U95, 
        ~ 0 ))
    
    summary_transect_coral <- left_join(
      coral_full_table,
      transect_coral_group,
      by = c("OCC_SITEID", "CB_TRANSECTID", "OCC_SITEID_TRANSECT")
      )
    
  } else {
    summary_transect_coral <-  left_join(
      coral_full_table,
      transect_coral_group,
      by = c(
        "OCC_SITEID",
        "CB_TRANSECTID",
        "OCC_SITEID_TRANSECT",
        "CORAL_GROUP",
        "CORAL_GROUP_NAME"
      )
    )
  } 
  
  summary_transect_coral <- summary_transect_coral %>%
    group_by(OCC_SITEID) %>%
    fill(REGION:TRANSECT_TOTAL_SUBSTRATE_COVER_M, .direction = "downup") %>%
    ungroup() %>%
    mutate(across(
      SUBSTRATE_COVER_CM:SUBSTRATE_CARB_PROD_KG_M2_YR_U95,
      ~ replace_na(.x, 0)
    )) 
  
  
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
  
  # Summarize production by TRANSECT ----------------------------------------
 
  if (macro_rates == "IPRB"){
    
    summary_transect_substratecode$MACROBIOEROSION_RATE_KG_M2_YR <- 0.209
    summary_transect_substratecode$MACROBIOEROSION_RATE_KG_M2_YR_CI <- 0.129
  }
  
  if (micro_rates == "IPRB"){
    
    summary_transect_substratecode$MICROBIOEROSION_RATE_KG_M2_YR <- 0.262
    summary_transect_substratecode$MICROBIOEROSION_RATE_KG_M2_YR_CI <- 0.180
  }
  
  if (macro_rates == "NCRMP"){
    
    summary_transect_substratecode <- summary_transect_substratecode %>%
      left_join(macro_dbase_ncrmp %>% select(LOCATIONCODE, MACROBIOEROSION_RATE_KG_M2_YR, MACROBIOEROSION_RATE_KG_M2_YR_CI),
                by = "LOCATIONCODE"
      )
  }
  
  if (micro_rates == "NCRMP"){
    
    summary_transect_substratecode$MICROBIOEROSION_RATE_KG_M2_YR <- 0.184
    summary_transect_substratecode$MICROBIOEROSION_RATE_KG_M2_YR_CI <- 0.121
  }
  
  
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
          SUBSTRATE_AVAILABLE_MACRO_INDEX * first(MACROBIOEROSION_RATE_KG_M2_YR),
        MACROBIOEROSION_KG_M2_YR_L95 =
          SUBSTRATE_AVAILABLE_MACRO_INDEX * (first(MACROBIOEROSION_RATE_KG_M2_YR) - first(MACROBIOEROSION_RATE_KG_M2_YR_CI)),
        MACROBIOEROSION_KG_M2_YR_U95 =
          SUBSTRATE_AVAILABLE_MACRO_INDEX * (first(MACROBIOEROSION_RATE_KG_M2_YR) + first(MACROBIOEROSION_RATE_KG_M2_YR_CI)),
        MICROBIOEROSION_KG_M2_YR =
          SUBSTRATE_AVAILABLE_MICRO_INDEX * first(MICROBIOEROSION_RATE_KG_M2_YR),
        MICROBIOEROSION_KG_M2_YR_L95 =
          SUBSTRATE_AVAILABLE_MICRO_INDEX * (first(MICROBIOEROSION_RATE_KG_M2_YR) - first(MICROBIOEROSION_RATE_KG_M2_YR_CI)),
        MICROBIOEROSION_KG_M2_YR_U95 =
          SUBSTRATE_AVAILABLE_MICRO_INDEX * (first(MICROBIOEROSION_RATE_KG_M2_YR) + first(MICROBIOEROSION_RATE_KG_M2_YR_CI)),
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
          SUBSTRATE_AVAILABLE_MACRO_INDEX * first(MACROBIOEROSION_RATE_KG_M2_YR),
        MACROBIOEROSION_KG_M2_YR_L95 =
          SUBSTRATE_AVAILABLE_MACRO_INDEX * (first(MACROBIOEROSION_RATE_KG_M2_YR) - first(MACROBIOEROSION_RATE_KG_M2_YR_CI)),
        MACROBIOEROSION_KG_M2_YR_U95 =
          SUBSTRATE_AVAILABLE_MACRO_INDEX * (first(MACROBIOEROSION_RATE_KG_M2_YR) + first(MACROBIOEROSION_RATE_KG_M2_YR_CI)),
        MICROBIOEROSION_KG_M2_YR =
          SUBSTRATE_AVAILABLE_MICRO_INDEX * first(MICROBIOEROSION_RATE_KG_M2_YR),
        MICROBIOEROSION_KG_M2_YR_L95 =
          SUBSTRATE_AVAILABLE_MICRO_INDEX * (first(MICROBIOEROSION_RATE_KG_M2_YR) - first(MICROBIOEROSION_RATE_KG_M2_YR_CI)),
        MICROBIOEROSION_KG_M2_YR_U95 =
          SUBSTRATE_AVAILABLE_MICRO_INDEX * (first(MICROBIOEROSION_RATE_KG_M2_YR) + first(MICROBIOEROSION_RATE_KG_M2_YR_CI)),
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


