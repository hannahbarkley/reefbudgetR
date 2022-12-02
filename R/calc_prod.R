#' Calculate production rates for each benthic component measured on
#' a carbonate budget transect
#'
#'@author Hannah Barkley
#'
#'@param substrate_class Type of substrate observed ("CORAL", "CCA", "TURF", "MA", ...).
#'@param substrate_code Taxa code observed (e.g. "PLOB", "MCAP", "PMEA", "CCA", ...).
#'@param morphology_code Taxa morphology observed (corals only; "BR", "MD", "EM", ...). Non-corals default to NA.
#'@param substrate_cover_cm Measure surface distance of benthic component, in cm.
#'@param region_code Survey region ("MHI", "MARIAN", ...).
#'@param prod_dbase Production database to reference, either Indo-Pacific ReefBudget ("IPRB")
#'or NCRMP-specific ("NCRMP"). Defaults to "NCRMP".
#'
#'@details See included carbonate production databases `prod_dbase_iprb` and `prod_dbase_ncrmp`
#'for acceptable `substrate_class`, `substrate_code`, and `morphology_code` values.
#'
#'@export calc_prod
#'
#'@examples
#' calc_prod(
#' substrate_class = "CORAL",
#' substrate_code = "PLOB",
#' morphology_code = "MD",
#' substrate_cover_cm = 10,
#' region_code = "MHI",
#' prod_dbase = "NCRMP")
#'
#'
calc_prod <- function(substrate_class,
                      substrate_code,
                      morphology_code,
                      substrate_cover_cm,
                      region_code,
                      prod_dbase) {

  # Set production rate for all non-calcifiers to zero
  if (substrate_class %in% c("CORAL", "CCA") == FALSE) {
    cp_i <- 0
    cp_l95_i <- 0
    cp_u95_i <- 0
  }

  # Calculate production for all calcifiers (corals and CCA)
  if (substrate_class %in% c("CORAL", "CCA") == TRUE) {

    if (substrate_class == "CCA") {

      # Get extension, density, and confidence intervals from the prod database
      g <-
        prod_dbase$EXTENSION_CM_YR[
        prod_dbase$SUBSTRATE_CODE == substrate_code &
        prod_dbase$REGIONCODE == region_code
          ]

      g_ci <-
        prod_dbase$EXTENSION_CM_YR_CI[
        prod_dbase$SUBSTRATE_CODE == substrate_code &
        prod_dbase$REGIONCODE == region_code
          ]

      d <-
        prod_dbase$DENSITY_G_CM3[
        prod_dbase$SUBSTRATE_CODE == substrate_code &
        prod_dbase$REGIONCODE == region_code
                                 ]

      d_ci <-
        prod_dbase$DENSITY_G_CM3_CI[
        prod_dbase$SUBSTRATE_CODE == substrate_code &
        prod_dbase$REGIONCODE == region_code]

      c <-
        prod_dbase$CONVERSION_FACTOR[
        prod_dbase$SUBSTRATE_CODE == substrate_code &
        prod_dbase$REGIONCODE == region_code]

      x <- substrate_cover_cm

      # Calculate production rate and lower and upper 95% confidence interval
      cp_i <- x * g * d
      cp_l95_i <- x * (g - g_ci) * d
      cp_u95_i <- x * (g + g_ci) * d

    }
    if (substrate_class == "CORAL") {
      # Get extension, density, and confidence intervals from the prod database
      g <-
        prod_dbase$EXTENSION_CM_YR[
          prod_dbase$SUBSTRATE_CODE == substrate_code &
          prod_dbase$MORPHOLOGYCODE == morphology_code &
          prod_dbase$REGIONCODE == region_code
          ]

      g_ci <-
        prod_dbase$EXTENSION_CM_YR_CI[
          prod_dbase$SUBSTRATE_CODE == substrate_code &
          prod_dbase$MORPHOLOGYCODE == morphology_code &
          prod_dbase$REGIONCODE == region_code
          ]

      d <-
        prod_dbase$DENSITY_G_CM3[
          prod_dbase$SUBSTRATE_CODE == substrate_code &
          prod_dbase$MORPHOLOGYCODE == morphology_code &
          prod_dbase$REGIONCODE == region_code
          ]

      d_ci <-
        prod_dbase$DENSITY_G_CM3_CI[
          prod_dbase$SUBSTRATE_CODE == substrate_code &
          prod_dbase$MORPHOLOGYCODE == morphology_code &
          prod_dbase$REGIONCODE == region_code
          ]

      c <-
        prod_dbase$CONVERSION_FACTOR[
          prod_dbase$SUBSTRATE_CODE == substrate_code &
          prod_dbase$MORPHOLOGYCODE == morphology_code &
          prod_dbase$REGIONCODE == region_code
          ]

      x <- substrate_cover_cm

      # Calculate production rate for mounding, mounding lobate,
      # and free-living morphologies
      if ((morphology_code %in% c("MD", "ML", "FR") == TRUE)) {

        df <- data.frame(x = seq(0, 135, 1))

        df$y <-  d * (((((((((df$x * 2) / pi
        ) / 2) + g) ^ 2
        ) * pi) / 2)) - (((((((df$x * 2) / pi) / 2) ^ 2
        ) * pi) / 2)))

        df$y_l95 <-  (d - d_ci) * (((((((((df$x * 2) / pi
        ) / 2) + (g - g_ci)) ^ 2
        ) * pi) / 2)) - (((((((df$x * 2) / pi) / 2) ^ 2
        ) * pi) / 2)))

        df$y_u95 <-  (d + d_ci) * (((((((((df$x * 2) / pi
        ) / 2) + (g + g_ci)) ^ 2
        ) * pi) / 2)) - (((((((df$x * 2) / pi) / 2) ^ 2
        ) * pi) / 2)))
      }

      # Calculate production rate for encrusting, plating, foliose,
      # and table morphologies
      if ((morphology_code %in% c("EM", "PL", "FO", "TB") == TRUE)) {

        h <- 2

        df <- data.frame(x = seq(0, 135, 1))

        df$y <- h * g * d + 0.1 * g * (df$x + (h * g)) * d

        df$y_l95 <-
          h * (g - g_ci) * (d - d_ci) + 0.1 * (g - g_ci) *
          (df$x + (h * (g - g_ci))) * (d - d_ci)

        df$y_u95 <-
          h * (g + g_ci) * (d + d_ci) + 0.1 * (g + g_ci) *
          (df$x + (h * (g + g_ci))) * (d + d_ci)
      }

      # Calculate production rate for branching, columnar,
      # and knobby morphologies
      if ((morphology_code %in% c("BR", "CO", "KN", "BRFA", "BRSL")) == TRUE) {

        df <- data.frame(x = seq(0, 135, 1))

        df$y <-  ((((1 - c) * df$x) * g * 0.1) + (c * df$x * g)) * d

        df$y_l95 <-
          ((((1 - c) * df$x) * (g - g_ci) * 0.1) +
             (c * df$x * (g - g_ci))) * (d - d_ci)

        df$y_u95 <-
          ((((1 - c) * df$x) * (g + g_ci) * 0.1) +
             (c * df$x * (g + g_ci))) * (d + d_ci)
      }

      # Calculate production rate for laminar columnar morphology (
      if ((morphology_code %in% c("LC")) == TRUE) {

        df <- data.frame(x = seq(0, 135, 1))

        # Laminar portion
        h <- 2

        prop_la <- 0.55 #Proportion of surface distance that is laminar

        df_la <- data.frame(x_la = seq(0, 135, 1))

        df_la$y_la <-
          h * g * d + 0.1 * g * ((prop_la * df_la$x_la) + (h * g)) * d

        df_la$y_la_l95 <-
          h * (g - g_ci) * (d - d_ci) + 0.1 * (g - g_ci) *
          ((prop_la * df_la$x_la) + (h * (g - g_ci))) * (d - d_ci)

        df_la$y_la_u95 <-
          h * (g + g_ci) * (d + d_ci) + 0.1 * (g + g_ci) *
          ((prop_la * df_la$x_la) + (h * (g + g_ci))) * (d + d_ci)

        # Columnar portion
        df_co <- data.frame(x_co = seq(0, 135, 1))

        prop_co <- 0.45 # Proportion of surface distance that is columnar

        df_co$y_co <-
          ((((1 - c) * (prop_co * df_co$x_co)) * g * 0.1) +
             (c * (prop_co * df_co$x_co) * g)) * d

        df_co$y_co_l95 <-
          ((((1 - c) * (prop_co * df_co$x_co)) * (g - g_ci) * 0.1) +
             (c * (prop_co * df_co$x_co) * (g - g_ci))) * (d - d_ci)

        df_co$y_co_u95 <-
          ((((1 - c) * (prop_co * df_co$x_co)) * (g + g_ci) * 0.1) +
             (c * (prop_co * df_co$x_co) * (g + g_ci))) * (d + d_ci)

        # Add together laminar and columnar portions to get overall production
        df$y <- df_la$y_la + df_co$y_co

        df$y_l95 <- df_la$y_la_l95 + df_co$y_co_l95

        df$y_u95 <- df_la$y_la_u95 + df_co$y_co_u95
      }

      df[1, ] <- 0

      lm <- lm(df$y ~ df$x)
      coefficient_mean <- lm$coefficients[2]
      intercept_mean <- lm$coefficients[1]
      cp_i <- coefficient_mean * x + intercept_mean

      lm_l95 <- lm(df$y_l95 ~ df$x)
      coefficient_mean_l95 <- lm_l95$coefficients[2]
      intercept_mean_l95 <- lm_l95$coefficients[1]
      cp_l95_i <- coefficient_mean_l95 * x + intercept_mean_l95

      lm_u95 <- lm(df$y_u95 ~ df$x)
      coefficient_mean_u95 <- lm_u95$coefficients[2]
      intercept_mean_u95 <- lm_u95$coefficients[1]
      cp_u95_i <- coefficient_mean_u95 * x + intercept_mean_u95

    }
  }
  return(data.frame(cp_i, cp_l95_i, cp_u95_i))
}
