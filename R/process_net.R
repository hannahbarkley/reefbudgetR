#' Calculate net production rates for each transect
#'@author Hannah Barkley
#'@param prod Production data
#'@param urch Urchin data
#'@param parrotfish Parrotfish data
#'@import dplyr
#'@import tools
#'@import tidyr
#'@import reshape2
#'@importFrom sjmisc seq_row
#'@export process_net

process_net <- function(prod,
                        urch,
                        fish,
                        format = "wide") {

    net_site_prod_urch <- left_join(prod,
                                    urch,
                                    by = colnames(prod)[colnames(prod) %in% colnames(urch)])

    for (i in sjmisc::seq_row(net_site_prod_urch)) {
      if (net_site_prod_urch$CB_METHOD[i] == "SfM") {
        net_site_prod_urch$URCHIN_EROSION_KG_M2_YR_MEAN[i] <-
          urch$URCHIN_EROSION_KG_M2_YR_MEAN[urch$OCC_SITEID ==
                                              net_site_prod_urch$OCC_SITEID[i] &
                                              urch$CB_METHOD == "Chords"]
        net_site_prod_urch$URCHIN_EROSION_KG_M2_YR_SD[i] <-
          urch$URCHIN_EROSION_KG_M2_YR_SD[urch$OCC_SITEID ==
                                            net_site_prod_urch$OCC_SITEID[i] &
                                            urch$CB_METHOD == "Chords"]
        net_site_prod_urch$URCHIN_EROSION_KG_M2_YR_N[i] <-
          urch$URCHIN_EROSION_KG_M2_YR_N[urch$OCC_SITEID ==
                                           net_site_prod_urch$OCC_SITEID[i] &
                                           urch$CB_METHOD == "Chords"]
        net_site_prod_urch$URCHIN_EROSION_KG_M2_YR_SE[i] <-
          urch$URCHIN_EROSION_KG_M2_YR_SE[urch$OCC_SITEID ==
                                            net_site_prod_urch$OCC_SITEID[i] &
                                            urch$CB_METHOD == "Chords"]
        net_site_prod_urch$URCHIN_EROSION_KG_M2_YR_CI[i] <-
          urch$URCHIN_EROSION_KG_M2_YR_CI[urch$OCC_SITEID ==
                                            net_site_prod_urch$OCC_SITEID[i] &
                                            urch$CB_METHOD == "Chords"]
      }
    }

    net_site_prod_urch <- net_site_prod_urch %>% mutate_at(vars(c("REGION":"CB_METHOD")), as.factor)


    net_site <- full_join(
      net_site_prod_urch,
      fish,
      by = c(
        "REGION",
        "REGIONCODE",
        "YEAR",
        "CRUISE_ID",
        "LOCATION",
        "LOCATIONCODE",
        "OCC_SITEID",
        "OCC_SITENAME",
        "LATITUTE",
        "LONGITUDE"
      )
    )

    colnames(net_site)[colnames(net_site) == "CB_METHOD.x"] <-
      "CB_METHOD_BENTHIC"
    colnames(net_site)[colnames(net_site) == "CB_METHOD.y"] <-
      "CB_METHOD_FISH"

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
           10)
    )

    if (format == "wide") {

      net_site <- format_4ncei(
        net_prod = net_site
      )

      return(summary_site = net_site)
    }

    if (format == "long") {

      net_sites_reshape = net_site[c(
        "REGIONCODE",
        "LOCATIONCODE",
        "OCC_SITENAME",
        "OCC_SITEID",
        "CB_METHOD_BENTHIC",
        "GROSS_CARB_PROD_KG_M2_YR_MEAN",
        "GROSS_CARB_PROD_KG_M2_YR_SE",
        "BIOEROSION_KG_M2_YR_MEAN",
        "BIOEROSION_KG_M2_YR_SE",
        "URCHIN_EROSION_KG_M2_YR_MEAN",
        "URCHIN_EROSION_KG_M2_YR_SE",
        "FISH_EROSION_KG_M2_YR_ALL_MEAN",
        "FISH_EROSION_KG_M2_YR_ALL_SE"
      )]

      net_sites_reshape$BIOEROSION_KG_M2_YR_MEAN = -1 * net_sites_reshape$BIOEROSION_KG_M2_YR_MEAN
      net_sites_reshape$URCHIN_EROSION_KG_M2_YR_MEAN  = -1 * net_sites_reshape$URCHIN_EROSION_KG_M2_YR_MEAN
      net_sites_reshape$FISH_EROSION_KG_M2_YR_ALL_MEAN = -1 * net_sites_reshape$FISH_EROSION_KG_M2_YR_ALL_MEAN



      net_sites_long_mean = reshape2::melt(
        net_sites_reshape,
        id.vars = c(
          "REGIONCODE",
          "LOCATIONCODE",
          "OCC_SITENAME",
          "OCC_SITEID",
          "CB_METHOD_BENTHIC"
        ),
        measure.vars = c(
          "GROSS_CARB_PROD_KG_M2_YR_MEAN",
          "BIOEROSION_KG_M2_YR_MEAN",
          "URCHIN_EROSION_KG_M2_YR_MEAN",
          "FISH_EROSION_KG_M2_YR_ALL_MEAN"
        ),
        variable.name = "PARAMETER",
        value.name = "MEAN"
      )

      levels(net_sites_long_mean$PARAMETER) = c(
        "Gross production",
        "Macro/microbioerosion",
        "Urchin erosion",
        "Parrotfish erosion"
      )

      net_sites_long_se = reshape2::melt(
        net_sites_reshape,
        id.vars = c(
          "REGIONCODE",
          "LOCATIONCODE",
          "OCC_SITENAME",
          "OCC_SITEID",
          "CB_METHOD_BENTHIC"
        ),
        measure.vars = c(
          "GROSS_CARB_PROD_KG_M2_YR_SE",
          "BIOEROSION_KG_M2_YR_SE",
          "URCHIN_EROSION_KG_M2_YR_SE",
          "FISH_EROSION_KG_M2_YR_ALL_SE"
        ),
        variable.name = "PARAMETER",
        value.name = "SE"
      )

      levels(net_sites_long_se$PARAMETER) = c(
        "Gross production",
        "Macro/microbioerosion",
        "Urchin erosion",
        "Parrotfish erosion"
      )

      net_sites_long = merge(
        net_sites_long_mean,
        net_sites_long_se,
        by = c(
          "REGIONCODE",
          "LOCATIONCODE",
          "OCC_SITENAME",
          "OCC_SITEID",
          "CB_METHOD_BENTHIC",
          "PARAMETER"
        )
      )
      return(summary_site = net_site_long)
    }

  }

