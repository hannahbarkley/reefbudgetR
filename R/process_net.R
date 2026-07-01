#' Calculate site-level net production rates
#'
#'@author Hannah Barkley
#'
#'@param prod Production data, output of `process_prod`.
#'@param urch Urchin data, output of `process_urchins`.
#'@param parrotfish Parrotfish data, output of `process_fish`.
#'@param sum_by Calculate net production from transect or site level data.
#'@param format Output data frame format ("wide" or "long"). Default is "wide".
#'
#'
#'@import dplyr
#'@import tools
#'@import tidyr
#'@import reshape2
#'@importFrom sjmisc seq_row
#'
#'@export process_net

process_net <- function(prod,
                        urch,
                        fish,
                        sum_by = c("site", "transect"),
                        format = "wide") {
  
  sum_by <- match.arg(sum_by)
  
  # Ensure Year is factor across all inputs 
  prod <- prod %>% mutate(YEAR = as.factor(YEAR))
  urch <- urch %>% mutate(YEAR = as.factor(YEAR))
  fish <- fish %>% mutate(YEAR = as.factor(YEAR))
  
  # Helper functions for standard error and confidence intervals
  se_fn <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  ci_fn <- function(x) se_fn(x) * 1.97
  
  # ----------------------------------------------------------------------------
  # TRANSECT LEVEL SUMMARY (Corrected Error Propagation & N_TRANSECTS)
  # ----------------------------------------------------------------------------
  if (sum_by == "transect") {
    
    # 1. Join ONLY Transect-level data (Prod + Urchin) first
    net_transect <- prod %>%
      mutate(YEAR = as.character(YEAR)) %>%
      left_join(
        urch %>% mutate(YEAR = as.character(YEAR)), 
        by = c("OCC_SITEID_TRANSECT", "YEAR"), 
        suffix = c("", ".y")
      ) %>%
      select(-ends_with(".y")) %>% 
      mutate(across(c(REGION, LOCALDATE, YEAR), as.factor)) %>%
      mutate(across(where(is.numeric), ~tidyr::replace_na(., 0)))
    
    # 2. Aggregate Transect data to Site level 
    net_site <- net_transect %>%
      group_by(REGION, REGIONCODE, YEAR, CRUISE_ID, LOCATION, LOCATIONCODE, 
               OCC_SITEID, SITEVISITID, LATITUDE, LONGITUDE, SITE_DEPTH_M, 
               LOCALDATE, METHOD) %>%
      summarise(
        N_TRANSECTS = n(), # <--- ADDED: Count number of transects per site
        across(
          c(
            RUGOSITY, HARD_CORAL_COVER_PCT, CCA_COVER_PCT,
            GROSS_CARB_PROD_KG_M2_YR, MACROBIOEROSION_KG_M2_YR,
            MICROBIOEROSION_KG_M2_YR, BIOEROSION_KG_M2_YR,
            HARD_CORAL_CARB_PROD_KG_M2_YR, CCA_CARB_PROD_KG_M2_YR,
            URCHIN_DENSITY_NO_M2, URCHIN_EROSION_KG_M2_YR
          ),
          list(
            MEAN = ~ mean(.x, na.rm = TRUE),
            SD   = ~ sd(.x, na.rm = TRUE),
            SE   = se_fn,
            CI95 = ci_fn,
            N    = ~ sum(!is.na(.x))
          )
        ),
        .groups = "drop"
      )
    
    # 3. NOW Join Fish data
    net_site <- net_site %>%
      left_join(
        fish %>% select(OCC_SITEID, starts_with("FISH_")), 
        by = "OCC_SITEID"
      ) %>%
      mutate(across(starts_with("FISH_"), ~tidyr::replace_na(., 0)))
    
    # 4. Final Math & Error Propagation (Quadrature)
    net_site <- net_site %>%
      mutate(
        # Means
        GROSS_EROSION_KG_M2_YR_MEAN = MACROBIOEROSION_KG_M2_YR_MEAN + 
          MICROBIOEROSION_KG_M2_YR_MEAN + 
          URCHIN_EROSION_KG_M2_YR_MEAN + 
          FISH_EROSION_KG_M2_YR_ALL_MEAN,
        
        NET_CARB_PROD_KG_M2_YR_MEAN = GROSS_CARB_PROD_KG_M2_YR_MEAN - GROSS_EROSION_KG_M2_YR_MEAN,
        
        # Standard Errors via Quadrature 
        GROSS_EROSION_KG_M2_YR_SE = sqrt(
          MACROBIOEROSION_KG_M2_YR_SE^2 + 
            MICROBIOEROSION_KG_M2_YR_SE^2 + 
            URCHIN_EROSION_KG_M2_YR_SE^2 + 
            FISH_EROSION_KG_M2_YR_ALL_SE^2
        ),
        
        NET_CARB_PROD_KG_M2_YR_SE = sqrt(
          GROSS_CARB_PROD_KG_M2_YR_SE^2 + 
            GROSS_EROSION_KG_M2_YR_SE^2
        )
      )
  }
  
  # ----------------------------------------------------------------------------
  # SITE LEVEL SUMMARY 
  # ----------------------------------------------------------------------------
  if (sum_by == "site") {
    
    net_site <- prod %>%
      full_join(urch, by = "OCC_SITEID", suffix = c("", ".y")) %>%
      full_join(fish, by = "OCC_SITEID", suffix = c("", ".y")) %>%
      select(-ends_with(".y")) %>%
      mutate(across(REGION:LOCALDATE, as.factor)) %>%
      mutate(across(where(is.numeric), ~tidyr::replace_na(., 0)))
    
    net_site <- net_site %>%
      mutate(
        N_TRANSECTS = GROSS_CARB_PROD_KG_M2_YR_N, 
        
        # Means
        GROSS_EROSION_KG_M2_YR_MEAN = MACROBIOEROSION_KG_M2_YR_MEAN + 
          MICROBIOEROSION_KG_M2_YR_MEAN + 
          URCHIN_EROSION_KG_M2_YR_MEAN + 
          FISH_EROSION_KG_M2_YR_ALL_MEAN,
        
        NET_CARB_PROD_KG_M2_YR_MEAN = GROSS_CARB_PROD_KG_M2_YR_MEAN - GROSS_EROSION_KG_M2_YR_MEAN,
        
        # Error Propagation via Quadrature 
        GROSS_EROSION_KG_M2_YR_SE = sqrt(
          MACROBIOEROSION_KG_M2_YR_SE^2 + 
            MICROBIOEROSION_KG_M2_YR_SE^2 + 
            URCHIN_EROSION_KG_M2_YR_SE^2 + 
            FISH_EROSION_KG_M2_YR_ALL_SE^2
        ),
        
        NET_CARB_PROD_KG_M2_YR_SE = sqrt(
          GROSS_CARB_PROD_KG_M2_YR_SE^2 + 
            GROSS_EROSION_KG_M2_YR_SE^2
        )
      )
  }
  
  # ----------------------------------------------------------------------------
  # FORMATTING
  # ----------------------------------------------------------------------------
  
  if (format == "wide") {
    return(format_4ncei(data = net_site))
  }
  
  if (format == "long") {
    
    net_long <- net_site %>%
      select(
        REGIONCODE, LOCATIONCODE, OCC_SITEID, N_TRANSECTS, 
        matches("GROSS_CARB_PROD_KG_M2_YR_(MEAN|SE)"),
        matches("MACROBIOEROSION_KG_M2_YR_(MEAN|SE)"),
        matches("MICROBIOEROSION_KG_M2_YR_(MEAN|SE)"),
        matches("URCHIN_EROSION_KG_M2_YR_(MEAN|SE)"),
        matches("FISH_EROSION_KG_M2_YR_ALL_(MEAN|SE)")
      ) %>%
      tidyr::pivot_longer(
        cols = -c(REGIONCODE, LOCATIONCODE, OCC_SITEID, N_TRANSECTS),
        names_to = c("PARAMETER_RAW", ".value"), 
        names_pattern = "(.*)_(MEAN|SE)"
      ) %>%
      mutate(
        MEAN = if_else(grepl("EROSION", PARAMETER_RAW), MEAN * -1, MEAN),
        PARAMETER = dplyr::case_match(PARAMETER_RAW,
                                      "GROSS_CARB_PROD_KG_M2_YR"    ~ "Gross production",
                                      "MACROBIOEROSION_KG_M2_YR"    ~ "Macrobioerosion",
                                      "MICROBIOEROSION_KG_M2_YR"    ~ "Microbioerosion",
                                      "URCHIN_EROSION_KG_M2_YR"     ~ "Urchin erosion",
                                      "FISH_EROSION_KG_M2_YR_ALL"   ~ "Parrotfish erosion",
                                      .default = PARAMETER_RAW
        )
      ) %>%
      select(REGIONCODE, LOCATIONCODE, OCC_SITEID, N_TRANSECTS, PARAMETER, MEAN, SE)
    
    return(net_long)
  }
}