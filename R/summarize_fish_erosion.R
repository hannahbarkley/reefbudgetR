#' Summarize site-level parrotfish bioerosion rates
#'
#' @author Rebecca Weible
#'
#' @param species_table Species table with average density, biomass, and bioerosion metrics,
#' outputs of `summarize_fish_metrics`.
#' @param full_summary Return full summary of results (default is TRUE)
#'
#' @import dplyr
#'
#' @export summarize_fish_erosion
#'
#' @examples
#'density_average <- summarize_fish_metrics(
#'  data = calc_eros_fish_output,
#'  metric = "density",
#'  level = "transect",
#'  summarize_by = "species"
#')

#'biomass_average <- suppressWarnings(
#'  summarize_fish_metrics(
#'    data = calc_eros_fish_output,
#'    metric = "biomass",
#'    level = "transect",
#'    summarize_by = "species"
#'  )

#'bioerosion_average <- suppressWarnings(
#'  summarize_fish_metrics(
#'    data = calc_eros_fish_output,
#'    metric = "bioerosion",
#'    level = "transect",
#'    summarize_by = "species"
#'  )

#'species_table <-
#'  rbind(density_average, biomass_average, bioerosion_average)

#' summary_fish_erosion <- summarize_fish_erosion(species_table, full_summary)
#'


summarize_fish_erosion <- function(species_table,
                                   full_summary = TRUE) {
  options(scipen = 999) #prevent scientific notation

  # Final bioerosion calculations per transect
  fish_erosion_transect_wide <- species_table %>%
    # Add in grazing type for final CB erosion rates
    left_join(.,
              fish_grazing_types %>%
                select(FISH_sciname,
                       SPECIES,
                       GRAZ_TYPE,
                       GENUS),
              by = "SPECIES") %>%
    # Re-order so Grazing Type comes first in column orders
    select(GRAZ_TYPE, everything()) %>%
    # Collapse dataframe to be workable with group_by()
    gather(
      .,
      "CB_TRANSECT_ID",
      "value",-c(GRAZ_TYPE:SPECIES),
      -GENUS,-FISH_sciname,-METRIC
    ) %>%
    dplyr::group_by(
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      OCC_SITENAME,
      LATITUDE,
      LONGITUDE,
      CB_METHOD,
      CB_TRANSECT_ID,
      GRAZ_TYPE,
      METRIC
    ) %>%
    # Summarize all species by Graze Type and Transect
    dplyr::summarize(.sum = sum(value)) %>%
    # Spread back out by Transect to look like CP Excel outputs
    spread(., CB_TRANSECT_ID, .sum, fill = NA) %>%
    # reorder Transect 10 to the end
        select(-c("TRANSECT_10"), everything()) %>%
    #exclude "other" if only want to look at excavators and scrapers
    # that contribute to bioerosion
    #filter(!GRAZ_TYPE %in% "Other") %>%
    bind_rows(
      species_table %>%
        group_by(
          REGION,
          REGIONCODE,
          YEAR,
          CRUISE_ID,
          LOCATION,
          LOCATIONCODE,
          OCC_SITEID,
          OCC_SITENAME,
          LATITUDE,
          LONGITUDE,
          CB_METHOD,
          METRIC
        ) %>%
        dplyr::summarise(across(
          TRANSECT_1:TRANSECT_10, ~ sum(.x, na.rm = T)
        )) %>%
        # add Total Row in dataframe for final CB calculation
        mutate(GRAZ_TYPE = "All")
    )



  # Calculate bioerosion per site and metric
  fish_erosion_transect <- fish_erosion_transect_wide %>%
    pivot_longer(.,
                 cols = 14:23,
                 names_to = "CB_TRANSECTID",
                 values_to = "Values") %>%
    select(
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      OCC_SITENAME,
      LATITUDE,
      LONGITUDE,
      CB_METHOD,
      METRIC,
      everything()
    ) %>%
    spread(., METRIC, Values, fill = 0)

  fish_erosion_transect$CB_TRANSECTID <- as.factor(str_split(fish_erosion_transect$CB_TRANSECTID, "\\_", simplify=T)[,2])

  fish_erosion_transect$CB_TRANSECTID <- factor(fish_erosion_transect$CB_TRANSECTID, levels = seq(1,10,1))

  fish_erosion_transect <- fish_erosion_transect[
    with(fish_erosion_transect, order(REGION, LOCATION, CB_METHOD, CB_TRANSECTID, GRAZ_TYPE)),
  ]

  # final bioerosion calculations per site and metric
  fish_erosion_site_long <- fish_erosion_transect_wide %>%
    select(
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      OCC_SITENAME,
      LATITUDE,
      LONGITUDE,
      CB_METHOD,
      METRIC,
      everything()
    ) %>%
    rowwise() %>%
    dplyr::mutate(MEAN = mean(c_across(TRANSECT_1:TRANSECT_10))) %>%
    dplyr::mutate(SD = sd(c_across(TRANSECT_1:TRANSECT_10))) %>%
    dplyr::mutate(SE = SD / sqrt(10)) %>% #10 is the number of total Transects
    dplyr::mutate(L95 = MEAN - (SE * 1.97)) %>%
    dplyr::mutate(U95 = MEAN + (SE * 1.97)) %>%
    select(-c(TRANSECT_1:TRANSECT_10)) #remove unnecessary columns

  fish_erosion_site_long$L95[fish_erosion_site_long$L95 < 0] <-
    0

  # format output dataframe for NCEI
  fish_erosion_site <- fish_erosion_site_long %>%
    select(
      REGION,
      REGIONCODE,
      YEAR,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      OCC_SITENAME,
      LATITUDE,
      LONGITUDE,
      CB_METHOD,
      METRIC,
      everything()
    ) %>%
    pivot_longer(.,
                 cols = 14:18,
                 names_to = "Variables",
                 values_to = "Values") %>%
    unite("METRIC", METRIC:GRAZ_TYPE) %>% #combine columns with underscore
    unite("METRIC", METRIC:Variables) %>% #combine columns with underscore
    spread(., METRIC, Values, fill = 0)

  names(fish_erosion_site) <- toupper(names(fish_erosion_site))

  if (full_summary == TRUE) {
    return(list(
      fish_erosion_transect = fish_erosion_transect,
      fish_erosion_site = fish_erosion_site)
    )
  }

  if (full_summary == FALSE) {
    return(fish_erosion_site)
  }


}
