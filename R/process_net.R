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

  
  if (sum_by == "site") {
    

    prod$YEAR <- as.factor(prod$YEAR)
    urch$YEAR <- as.factor(urch$YEAR)
    
    net_site_prod_urch <- full_join(prod,
                                    urch,
                                        by = "OCC_SITEID",
                                        suffix = c("",".y"),
                                        relationship = "one-to-one") %>% select(-ends_with(".y"))
    
    net_site_prod_urch <-
      net_site_prod_urch %>% mutate_at(vars(c("REGION":"LOCALDATE")), as.factor)

    fish$LATITUDE <- net_site_prod_urch$LATITUDE[match(fish$OCC_SITEID, net_site_prod_urch$OCC_SITEID)]
    fish$LONGITUDE <- net_site_prod_urch$LONGITUDE[match(fish$OCC_SITEID, net_site_prod_urch$OCC_SITEID)]
    
    net_site <- full_join(
      net_site_prod_urch,
      fish,
      by = "OCC_SITEID",
      suffix = c("",".y"),
      relationship = "one-to-one") %>% select(-ends_with(".y"))
    
    net_site[, c((which(
      colnames(net_site) == "HARD_CORAL_CARB_PROD_KG_M2_YR_MEAN"
    )):(which(
      colnames(net_site) == "FISH_EROSION_KG_M2_YR_ALL_SE"
    )))][is.na(net_site[, c((which(
      colnames(net_site) == "HARD_CORAL_CARB_PROD_KG_M2_YR_MEAN"
    )):(which(
      colnames(net_site) == "FISH_EROSION_KG_M2_YR_ALL_SE"
    )))])] <- 0
    
    net_site$GROSS_EROSION_KG_M2_YR_MEAN <-
      net_site$MACROBIOEROSION_KG_M2_YR_MEAN +
      net_site$MICROBIOEROSION_KG_M2_YR_MEAN +
      net_site$URCHIN_EROSION_KG_M2_YR_MEAN +
      net_site$FISH_EROSION_KG_M2_YR_ALL_MEAN
    
    net_site$GROSS_EROSION_KG_M2_YR_SD <- sqrt(
      (net_site$GROSS_CARB_PROD_KG_M2_YR_SD ^ 2) +
        (net_site$BIOEROSION_KG_M2_YR_SD ^ 2) +
        (net_site$URCHIN_EROSION_KG_M2_YR_SD ^ 2) +
        (net_site$FISH_EROSION_KG_M2_YR_ALL_SD ^ 2)
    )

    net_site$NET_CARB_PROD_KG_M2_YR_MEAN <-
      net_site$GROSS_CARB_PROD_KG_M2_YR_MEAN -
      net_site$BIOEROSION_KG_M2_YR_MEAN -
      net_site$URCHIN_EROSION_KG_M2_YR_MEAN -
      net_site$FISH_EROSION_KG_M2_YR_ALL_MEAN

    net_site$NET_CARB_PROD_KG_M2_YR_SD <- sqrt(
      (net_site$GROSS_CARB_PROD_KG_M2_YR_SD ^ 2) +
        (net_site$BIOEROSION_KG_M2_YR_SD ^ 2) +
        (net_site$URCHIN_EROSION_KG_M2_YR_SD ^ 2) +
        (net_site$FISH_EROSION_KG_M2_YR_ALL_SD ^ 2)
    )

    net_site$NET_CARB_PROD_KG_M2_YR_SE <- sqrt(
      (
        net_site$GROSS_CARB_PROD_KG_M2_YR_SD ^ 2 /
          net_site$GROSS_CARB_PROD_KG_M2_YR_N
      ) +
        (
          net_site$BIOEROSION_KG_M2_YR_SD ^ 2 /
            net_site$BIOEROSION_KG_M2_YR_N
        ) +
        (
          net_site$URCHIN_EROSION_KG_M2_YR_SD ^ 2 /
            net_site$URCHIN_EROSION_KG_M2_YR_N
        ) +
        (net_site$FISH_EROSION_KG_M2_YR_ALL_SD ^ 2 /
           1)
    )
  }

  if (sum_by == "transect") {
    
    net_transect_prod_urch <- left_join(prod,
                                        urch,
                                        by = c("OCC_SITEID_TRANSECT", "YEAR"),
                                        suffix = c("",".y"),
                                        relationship = "one-to-one",
                                        unmatched = "drop") %>% select(-ends_with(".y"))
    
    
    fish$LATITUDE <-
      as.factor(net_transect_prod_urch$LATITUDE[match(fish$OCC_SITEID, net_transect_prod_urch$OCC_SITEID)])
    fish$LONGITUDE <-
      as.factor(net_transect_prod_urch$LONGITUDE[match(fish$OCC_SITEID, net_transect_prod_urch$OCC_SITEID)])
    fish$YEAR <-
      as.factor(fish$YEAR)

    net_transect_prod_urch <-
      net_transect_prod_urch %>% mutate_at(vars(c("REGION","LOCALDATE","YEAR")), as.factor)
    
    net_transect <- left_join(
      net_transect_prod_urch,
      fish,
      by = c("OCC_SITEID", "YEAR"),
      suffix = c("",".y"),
      relationship = "many-to-one",
      unmatched = "drop") %>% select(-ends_with(".y"))
    
    net_transect[, c((which(
      colnames(net_transect) == "URCHIN_ABUNDANCE_NO"
    )):(which(
      colnames(net_transect) == "URCHIN_EROSION_KG_M2_YR"
    )))][is.na(net_transect[, c((which(
      colnames(net_transect) == "URCHIN_ABUNDANCE_NO"
    )):(which(
      colnames(net_transect) == "URCHIN_EROSION_KG_M2_YR"
    )))])] <- 0
    
    net_transect[, c((which(
      colnames(net_transect) == "FISH_BIOMASS_KG_HA_ALL_L95"
    )):(which(
      colnames(net_transect) == "FISH_EROSION_KG_M2_YR_SCRAPER_U95"
    )))][is.na(net_transect[, c((which(
      colnames(net_transect) == "FISH_BIOMASS_KG_HA_ALL_L95"
    )):(which(
      colnames(net_transect) == "FISH_EROSION_KG_M2_YR_SCRAPER_U95"
    )))])] <- 0
    
    net_transect$GROSS_EROSION_KG_M2_YR <-
      net_transect$MACROBIOEROSION_KG_M2_YR +
      net_transect$MICROBIOEROSION_KG_M2_YR +
      net_transect$URCHIN_EROSION_KG_M2_YR +
      net_transect$FISH_EROSION_KG_M2_YR_ALL_MEAN
    
    net_transect$NET_CARB_PROD_KG_M2_YR <-
      net_transect$GROSS_CARB_PROD_KG_M2_YR -
      net_transect$MACROBIOEROSION_KG_M2_YR -
      net_transect$MICROBIOEROSION_KG_M2_YR -
      net_transect$URCHIN_EROSION_KG_M2_YR -
      net_transect$FISH_EROSION_KG_M2_YR_ALL_MEAN
    

    
    
    se <- function(x) {
      sd(x, na.rm = TRUE) / sqrt(length(x))
    }
    ci95 <- function(x) {
      sd(x, na.rm = TRUE) / sqrt(length(x)) * 1.97
    }
    
    net_site <- net_transect %>%
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
        METHOD
      ) %>%
      dplyr::reframe(across(
        c(
          RUGOSITY,
          HARD_CORAL_COVER_PCT,
          CCA_COVER_PCT,
          GROSS_CARB_PROD_KG_M2_YR,
          MACROBIOEROSION_KG_M2_YR,
          MICROBIOEROSION_KG_M2_YR,
          BIOEROSION_KG_M2_YR,
          HARD_CORAL_CARB_PROD_KG_M2_YR,
          CCA_CARB_PROD_KG_M2_YR,
          URCHIN_DENSITY_NO_M2,
          URCHIN_EROSION_KG_M2_YR,
          FISH_BIOMASS_KG_HA_ALL_MEAN,
          FISH_DENSITY_ABUNDANCE_HA_ALL_MEAN,
          FISH_EROSION_KG_M2_YR_ALL_MEAN,
          GROSS_EROSION_KG_M2_YR,
          NET_CARB_PROD_KG_M2_YR
        ),
        list(
          MEAN = ~ mean(.x, na.rm = TRUE),
          SD = ~ sd(.x, na.rm = TRUE),
          SE = se,
          CI95 = ci95,
          N = length
        )
      ))
    
    
    
    net_site$FISH_BIOMASS_KG_HA_ALL_MEAN <-
      fish_site$FISH_BIOMASS_KG_HA_ALL_MEAN[match(net_site$OCC_SITEID, fish_site$OCC_SITEID)]
    net_site$FISH_BIOMASS_KG_HA_ALL_SE  <-
      fish_site$FISH_BIOMASS_KG_HA_ALL_SE[match(net_site$OCC_SITEID, fish_site$OCC_SITEID)]
    net_site$FISH_DENSITY_ABUNDANCE_HA_ALL_MEAN <-
      fish_site$FISH_DENSITY_ABUNDANCE_HA_ALL_MEAN[match(net_site$OCC_SITEID, fish_site$OCC_SITEID)]
    net_site$FISH_DENSITY_ABUNDANCE_HA_ALL_SE <-
      fish_site$FISH_DENSITY_ABUNDANCE_HA_ALL_SE[match(net_site$OCC_SITEID, fish_site$OCC_SITEID)]
    net_site$FISH_EROSION_KG_M2_YR_ALL_MEAN <-
      fish_site$FISH_EROSION_KG_M2_YR_ALL_MEAN[match(net_site$OCC_SITEID, fish_site$OCC_SITEID)]
    net_site$FISH_EROSION_KG_M2_YR_ALL_SE <-
      fish_site$FISH_EROSION_KG_M2_YR_ALL_SE[match(net_site$OCC_SITEID, fish_site$OCC_SITEID)]

    net_site$EROSION_KG_M2_YR_SE <- sqrt(net_site$MACROBIOEROSION_KG_M2_YR_SE^2 +
                                           net_site$MICROBIOEROSION_KG_M2_YR_SE^2 +
                                           net_site$URCHIN_EROSION_KG_M2_YR_SE^2 +
                                           net_site$FISH_EROSION_KG_M2_YR_ALL_SE^2)

  }

  if (format == "wide") {
    net_site <- format_4ncei(data = net_site)

    return(summary_site = net_site)
  }

  if (format == "long") {
    net_sites_reshape = net_site[c(
      "REGIONCODE",
      "LOCATIONCODE",
      "OCC_SITEID",
      "GROSS_CARB_PROD_KG_M2_YR_MEAN",
      "GROSS_CARB_PROD_KG_M2_YR_SE",
      "MACROBIOEROSION_KG_M2_YR_MEAN",
      "MACROBIOEROSION_KG_M2_YR_SE",
      "MICROBIOEROSION_KG_M2_YR_MEAN",
      "MICROBIOEROSION_KG_M2_YR_SE",
      "URCHIN_EROSION_KG_M2_YR_MEAN",
      "URCHIN_EROSION_KG_M2_YR_SE",
      "FISH_EROSION_KG_M2_YR_ALL_MEAN",
      "FISH_EROSION_KG_M2_YR_ALL_SE"
    )]

    net_sites_reshape$MACROBIOEROSION_KG_M2_YR_MEAN = -1 * net_sites_reshape$MACROBIOEROSION_KG_M2_YR_MEAN
    net_sites_reshape$MICROBIOEROSION_KG_M2_YR_MEAN = -1 * net_sites_reshape$MICROBIOEROSION_KG_M2_YR_MEAN
    net_sites_reshape$URCHIN_EROSION_KG_M2_YR_MEAN  = -1 * net_sites_reshape$URCHIN_EROSION_KG_M2_YR_MEAN
    net_sites_reshape$FISH_EROSION_KG_M2_YR_ALL_MEAN = -1 * net_sites_reshape$FISH_EROSION_KG_M2_YR_ALL_MEAN



    net_sites_long_mean = reshape2::melt(
      net_sites_reshape,
      id.vars = c(
        "REGIONCODE",
        "LOCATIONCODE",
        "OCC_SITEID"
      ),
      measure.vars = c(
        "GROSS_CARB_PROD_KG_M2_YR_MEAN",
        "MACROBIOEROSION_KG_M2_YR_MEAN",
        "MICROBIOEROSION_KG_M2_YR_MEAN",
        "URCHIN_EROSION_KG_M2_YR_MEAN",
        "FISH_EROSION_KG_M2_YR_ALL_MEAN"
      ),
      variable.name = "PARAMETER",
      value.name = "MEAN"
    )

    levels(net_sites_long_mean$PARAMETER) = c(
      "Gross production",
      "Macrobioerosion",
      "Microbioerosion",
      "Urchin erosion",
      "Parrotfish erosion"
    )

    net_sites_long_se = reshape2::melt(
      net_sites_reshape,
      id.vars = c(
        "REGIONCODE",
        "LOCATIONCODE",
        "OCC_SITEID"
      ),
      measure.vars = c(
        "GROSS_CARB_PROD_KG_M2_YR_SE",
        "MACROBIOEROSION_KG_M2_YR_SE",
        "MICROBIOEROSION_KG_M2_YR_SE",
        "URCHIN_EROSION_KG_M2_YR_SE",
        "FISH_EROSION_KG_M2_YR_ALL_SE"
      ),
      variable.name = "PARAMETER",
      value.name = "SE"
    )

    levels(net_sites_long_se$PARAMETER) = c(
      "Gross production",
      "Macrobioerosion",
      "Microbioerosion",
      "Urchin erosion",
      "Parrotfish erosion"
    )

    net_sites_long = merge(
      net_sites_long_mean,
      net_sites_long_se,
      by = c(
        "REGIONCODE",
        "LOCATIONCODE",
        "OCC_SITEID",
        "PARAMETER"
      )
    )
    return(summary_site = net_sites_long)
  }

}
