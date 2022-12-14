#' Run calc_prod function
#'
#'@author Hannah Barkley
#'
#'@param data Benthic data set to process
#'@param transect_id String of transect names (e.g., ("A1", "A2", "A3", "B1", "B2", "B3")). Defaults to NULL.
#'@param transect_length String of transect lengths in meters (e.g., c(10, 10, 10, 10, 10, 10)). Defaults to NULL.
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
#'     transect_id = c("A1", "A2", "A3", "B1", "B2", "B3"),
#'     transect_length = c(10, 10, 10, 10, 10, 10),
#'     method_name = "IPRB",
#'     dbase_type = "NCRMP"
#'     )


run_calc_prod <- function(data,
                          transect_id,
                          transect_length,
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

  if(is.null(transect_id) == TRUE & is.null(transect_length) == TRUE){
    transect_summary_pairs <- unique(data[c("CB_TRANSECTID","TRANSECT_PLANAR_LENGTH_M")])
    transect_id <- unique(transect_summary_pairs$CB_TRANSECTID)
    transect_length <- unique(transect_summary_pairs$TRANSECT_PLANAR_LENGTH_M)
  }

  if (method_name == "SfM") {
    transect_summary <- summarize_transect(data,
                                          transect_id,
                                          transect_length,
                                          method_name = "SfM")
  }
  if (method_name == "IPRB") {
    transect_summary <- summarize_transect(data,
                                          transect_id,
                                          transect_length,
                                          method_name = "IPRB")
  }

  if (method_name == "Chords") {
    transect_summary <- summarize_transect(data,
                                           transect_id,
                                           transect_length,
                                           method_name = "Chords")
  }

  data$SUBSTRATE_CLASS <- NA
  data$SUBSTRATE_NAME <- NA
  data$TAXA_LEVEL <- NA

  for (i in sjmisc::seq_row(data)) {

    substrate_code_i <- data$SUBSTRATE_CODE[i]

    data$SUBSTRATE_NAME[i] <-
      unique(prod_dbase$SUBSTRATE_NAME[
        prod_dbase$SUBSTRATE_CODE == substrate_code_i])

    data$SUBSTRATE_CLASS[i] <-
      unique(prod_dbase$SUBSTRATE_CLASS[
        prod_dbase$SUBSTRATE_CODE == substrate_code_i])

    data$TAXA_LEVEL[i] <-
      unique(prod_dbase$TAXA_LEVEL[
        prod_dbase$SUBSTRATE_CODE == substrate_code_i])

    if (data$SUBSTRATE_CLASS[i] == "CORAL") {

      if (data$TAXA_LEVEL[i] == "SPECIES") {

        if (data$SUBSTRATE_CODE[i] == "MCAP") {

          data$MORPHOLOGYCODE[i] <- data$MORPHOLOGYCODE[i]

        } else {

          data$MORPHOLOGYCODE[i] <-
            as.character(paste0(unique(
              prod_dbase$MORPHOLOGYCODE[
                prod_dbase$SUBSTRATE_CODE == substrate_code_i])))
        }
      }
    }

    data$CORAL_GROUP[i] <-
      unique(prod_dbase$CORAL_GROUP[
        prod_dbase$SUBSTRATE_CODE == substrate_code_i &
          prod_dbase$MORPHOLOGYCODE == data$MORPHOLOGYCODE[i]])

    data$CORAL_GROUP_NAME[i] <-
     unique(prod_dbase$CORAL_GROUP_NAME[
        prod_dbase$SUBSTRATE_CODE == substrate_code_i &
          prod_dbase$MORPHOLOGYCODE == data$MORPHOLOGYCODE[i]])

    data$MORPHOLOGY[i] <-
     unique(prod_dbase$MORPHOLOGY[
        prod_dbase$SUBSTRATE_CODE == substrate_code_i &
          prod_dbase$MORPHOLOGYCODE == data$MORPHOLOGYCODE[i]])

    calc_i <- calc_prod(
      data$SUBSTRATE_CLASS[i],
      data$SUBSTRATE_CODE[i],
      data$MORPHOLOGYCODE[i],
      data$SUBSTRATE_COVER_CM[i],
      data$REGIONCODE[i],
      prod_dbase
    )

    data$COLONY_PROD_G_YR[i] <- calc_i$cp_i
    data$COLONY_PROD_G_YR_L95[i] <- calc_i$cp_l95_i
    data$COLONY_PROD_G_YR_U95[i] <- calc_i$cp_u95_i

  }

  data$OCC_SITEID_TRANSECT <-
    paste(data$OCC_SITEID, data$CB_TRANSECTID, sep = "-")

  data$TRANSECT_PLANAR_LENGTH_M <-
    transect_summary$TRANSECT_PLANAR_LENGTH_M[match(
      data$OCC_SITEID_TRANSECT,
      transect_summary$OCC_SITEID_TRANSECT)]

  data$TRANSECT_TOTAL_SUBSTRATE_COVER_M <-
    transect_summary$TRANSECT_TOTAL_SUBSTRATE_COVER_M[match(
      data$OCC_SITEID_TRANSECT,
      transect_summary$OCC_SITEID_TRANSECT)]

  return(list(
    data = data,
    transect_summary = transect_summary
    )
  )
}
