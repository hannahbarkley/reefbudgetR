#' Summarize benthic transect metadata
#'
#'@author Hannah Barkley
#'
#'@param data Production data set.
#'@param transect_id String of transect names (e.g., ("A1", "A2", "A3", "B1", "B2", "B3")).
#'@param transect_length String of transect lengths in meters (e.g., c(10, 10, 10, 10, 10, 10)).
#'@param data_type Type of data collection ("In water" or "SfM").
#'
#'@import dplyr
#'@importFrom rlang .data
#'
#'@export summarize_transect
#'
summarize_transect <-
  function(data,
           transect_id,
           transect_length,
           data_type = c("In water", "SfM")) {

    options(dplyr.summarise.inform = FALSE)

    if (data_type == "In water") {
      transects <- data.frame(transect_id, transect_length)

      data$TRANSECT_PLANAR_LENGTH_M <-
        transects$transect_length[match(data$CB_TRANSECTID,
                                        transects$transect_id)]
      transect_summary <-
        data %>%
        dplyr::group_by(
          .data$REGIONCODE,
          .data$LOCATIONCODE,
          .data$OCC_SITEID,
          .data$CB_METHOD,
          .data$CB_TRANSECTID
        ) %>%
        summarize(
          TRANSECT_PLANAR_LENGTH_M = mean(.data$TRANSECT_PLANAR_LENGTH_M),
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

      return(transect_summary)

    }
    if (data_type == "SfM") {

      transect_summary <-
        data %>%
        dplyr::group_by(.data$REGIONCODE,
                 .data$LOCATIONCODE,
                 .data$OCC_SITEID,
                 .data$CB_TRANSECTID) %>%
        summarize(
          TRANSECT_PLANAR_LENGTH_M = sum(.data$LINEAR_METER),
          TRANSECT_TOTAL_SUBSTRATE_COVER_M  =
            sum(.data$SUBSTRATE_COVER_CM / 100, na.rm = TRUE),
          TRANSECT_RUGOSITY =
            .data$TRANSECT_TOTAL_SUBSTRATE_COVER_M /
            .data$TRANSECT_PLANAR_LENGTH_M
        )

      transect_summary$OCC_SITEID_TRANSECT <-
        paste(transect_summary$OCC_SITEID,
              transect_summary$CB_TRANSECTID,
              sep = "-")

      return(transect_summary)
    }
  }
