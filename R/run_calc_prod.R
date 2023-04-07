#' Run calc_prod function
#'
#'@author Hannah Barkley
#'
#'@param data Benthic data set to process
#'@param dbase_type Production database to use, either Indo-Pacific ReefBudget ("IPRB")
#'or U.S. Pacific Islands NCRMP-specific database ("NCRMP"). The Indo-Pacific ReefBudget
#'database is derived from "IP Calcification and bioerosion rates database v.1.3",
#'downloaded from https://geography.exeter.ac.uk/reefbudget/indopacific/.
#'@param method_name Transect design by which data were collected ("IPRB", "Chords", or "SfM").
#'
#'@import dplyr
#'@importFrom sjmisc seq_row
#'
#'@export run_calc_prod
#'
#' @examples
#' benthic_data <- read.csv("ESD_CarbBudget_Benthic_OAHU_2021.csv",
#'     na = "", check.names = FALSE)
#'
#' calc_prod_output <- run_calc_prod(
#'     data = benthic_data,
#'     method_name = "IPRB",
#'     dbase_type = "NCRMP"
#'     )


run_calc_prod <- function(data,
                          dbase_type,
                          method_name,
                          ...) {
  if (dbase_type == "IPRB") {
    data$SUBSTRATE_CODE <- data$SUBSTRATE_CODE_IPRB
    prod_dbase <- prod_dbase_iprb
  }

  if (dbase_type == "NCRMP") {
    prod_dbase <- prod_dbase_ncrmp
  }

  transects <-
    unique(data[c("OCC_SITEID_TRANSECT", "TRANSECT_PLANAR_LENGTH_M")])
  transect_summary <-
    data %>%
    dplyr::group_by(
      .data$REGIONCODE,
      .data$LOCATIONCODE,
      .data$OCC_SITEID,
      .data$CB_METHOD,
      .data$CB_TRANSECTID,
      .data$TRANSECT_PLANAR_LENGTH_M
    ) %>%
    summarize(
      TRANSECT_PLANAR_LENGTH_M = mean(as.numeric(.data$TRANSECT_PLANAR_LENGTH_M)),
      TRANSECT_TOTAL_SUBSTRATE_COVER_M =
        sum(.data$SUBSTRATE_COVER_CM / 100, na.rm = TRUE),
      TRANSECT_RUGOSITY =
        .data$TRANSECT_TOTAL_SUBSTRATE_COVER_M /
        .data$TRANSECT_PLANAR_LENGTH_M
    )

  transect_summary$OCC_SITEID_TRANSECT <-
    paste(transect_summary$OCC_SITEID,
          transect_summary$CB_TRANSECTID,
          sep = "-")

  data$SUBSTRATE_CLASS <- NA
  data$SUBSTRATE_NAME <- NA
  data$TAXA_LEVEL <- NA
  data$CORAL_GROUP_NAME <- NA
  data$CORAL_GROUP <- NA


  for (i in sjmisc::seq_row(data)) {
    substrate_code_i <- data$SUBSTRATE_CODE[i]

    data$SUBSTRATE_NAME[i] <-
      unique(prod_dbase$SUBSTRATE_NAME[prod_dbase$SUBSTRATE_CODE == substrate_code_i])

    data$SUBSTRATE_CLASS[i] <-
      unique(prod_dbase$SUBSTRATE_CLASS[prod_dbase$SUBSTRATE_CODE == substrate_code_i])

    data$TAXA_LEVEL[i] <-
      unique(prod_dbase$TAXA_LEVEL[prod_dbase$SUBSTRATE_CODE == substrate_code_i])

    if (data$SUBSTRATE_CODE[i] == "PRUS") {
      data$MORPHOLOGYCODE[i] <- "LC"

      data$CORAL_GROUP[i] <- "POLC"

      data$CORAL_GROUP_NAME[i] <- "Laminar columnar Porites"

      data$MORPHOLOGY[i] <- "Laminar columnar"
    }

    if (data$SUBSTRATE_CLASS[i] == "CORAL" &
        data$TAXA_LEVEL[i] == "SPECIES" &
        data$SUBSTRATE_CODE[i] != "PRUS") {
      data$MORPHOLOGYCODE[i] <-
        as.character(paste0(unique(prod_dbase$MORPHOLOGYCODE[prod_dbase$SUBSTRATE_CODE == substrate_code_i])))

      data$CORAL_GROUP[i] <-
        unique(prod_dbase$CORAL_GROUP[prod_dbase$SUBSTRATE_CODE == substrate_code_i &
                                        prod_dbase$MORPHOLOGYCODE == data$MORPHOLOGYCODE[i]])

      data$CORAL_GROUP_NAME[i] <-
        (unique(prod_dbase$CORAL_GROUP_NAME[prod_dbase$SUBSTRATE_CODE == substrate_code_i &
                                              prod_dbase$MORPHOLOGYCODE == data$MORPHOLOGYCODE[i]]))

      data$MORPHOLOGY[i] <-
        unique(prod_dbase$MORPHOLOGY[prod_dbase$SUBSTRATE_CODE == substrate_code_i &
                                       prod_dbase$MORPHOLOGYCODE == data$MORPHOLOGYCODE[i]])

    }

    if (dbase_type == "IPRB") {
      calc_i <- calc_prod(
        substrate_class = data$SUBSTRATE_CLASS[i],
        substrate_code = data$SUBSTRATE_CODE[i],
        morphology_code = data$MORPHOLOGYCODE[i],
        substrate_cover_cm = data$SUBSTRATE_COVER_CM[i],
        dbase_type = "IPRB",
        prod_dbase = prod_dbase
      )
    }

    if (dbase_type == "NCRMP") {
      calc_i <- calc_prod(
        substrate_class = data$SUBSTRATE_CLASS[i],
        substrate_code = data$SUBSTRATE_CODE[i],
        morphology_code = data$MORPHOLOGYCODE[i],
        substrate_cover_cm = data$SUBSTRATE_COVER_CM[i],
        dbase_type = "NCRMP",
        prod_dbase = prod_dbase
      )
    }

    data$COLONY_PROD_G_YR[i] <- calc_i$cp_i
    data$COLONY_PROD_G_YR_L95[i] <- calc_i$cp_l95_i
    data$COLONY_PROD_G_YR_U95[i] <- calc_i$cp_u95_i

  }


  data$OCC_SITEID_TRANSECT <-
    paste(data$OCC_SITEID, data$CB_TRANSECTID, sep = "-")

  data$TRANSECT_PLANAR_LENGTH_M <-
    transect_summary$TRANSECT_PLANAR_LENGTH_M[match(data$OCC_SITEID_TRANSECT,
                                                    transect_summary$OCC_SITEID_TRANSECT)]

  data$TRANSECT_TOTAL_SUBSTRATE_COVER_M <-
    transect_summary$TRANSECT_TOTAL_SUBSTRATE_COVER_M[match(data$OCC_SITEID_TRANSECT,
                                                            transect_summary$OCC_SITEID_TRANSECT)]

  return(list(data = data,
              transect_summary = transect_summary))
}
