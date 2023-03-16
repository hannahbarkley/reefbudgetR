#' Calculate metrics for density, biomass, and bioerosion at the transect
#' and site level
#'
#'@author Rebecca Weible
#'
#'@param data output from `calc_eros_fish_output`.
#'@param metric Metric to summarize ("count", "biomass", or "bioerosion").
#'@param level Level to summarize at ("transect" or "site" level).
#'@param summarize_by Grouping factor to summarize by ("size class",
#'"species", or "overall").
#'
#'@import Rmisc
#'@import dplyr
#'
#'@export summarize_fish_metrics
#'
#'@examples
#'calc_eros_fish_output <- calc_eros_fish(data, rates_dbase = rates_dbase)
#'
#'@examples
#'density_average <- summarize_fish_metrics(
#'  data = calc_eros_fish_output,
#'  metric = "density",
#'  level = "transect",
#'  summarize_by = "species"
#')

#'biomass_average <- summarize_fish_metrics(
#'  data = calc_eros_fish_output,
#'  metric = "biomass",
#'  level = "transect",
#'  summarize_by = "species"
#'  )

#'bioerosion_average <- summarize_fish_metrics(
#'  data = calc_eros_fish_output,
#'  metric = "bioerosion",
#'  level = "transect",
#'  summarize_by = "species"
#'  )


summarize_fish_metrics <- function(data,
                                   metric = c("density",
                                              "biomass",
                                              "bioerosion"),
                                   level = c("transect",
                                             "site"),
                                   summarize_by = c("size class",
                                                    "species",
                                                    "overall")) {
  if (metric == "density") {
    metric <- "COUNT"
  }
  if (metric == "biomass") {
    metric <- "BIOMASS_KG_HA"
  }
  if (metric == "bioerosion") {
    metric <- "FISH_EROSION_KG_M2_YR"
  }


  # Summarize SIZE CLASS at TRANSECT level ---------------------------------------

  summary_transect_sizeclass <- data %>%
    select(CRUISE_ID:SIZE_CLASS, all_of(metric)) %>%
    ungroup() %>%
    dplyr::group_by(
      REGION,
      REGIONCODE,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      OCC_SITENAME,
      LATITUDE,
      LONGITUDE,
      CB_METHOD,
      TRANSECT,
      SIZE_CLASS,
      PHASE) %>%
    #calculate density/biomass/bioerosion per size class in individuals per
    # hectare (converted from m^2 to hectare by /10000)
    dplyr::reframe('SIZE_CLASS_SUM' = (sum(!!sym(metric)) /
                                              (AREA_M2 / 10000))) %>%
    distinct() %>%
    #organize dataframe by Transect and size bin
    pivot_wider(names_from = TRANSECT, values_from = SIZE_CLASS_SUM, names_prefix = paste0("TRANSECT_")) %>%
    mutate_at(-c(1:12), ~ replace_na(., 0))


  # Summarize SIZE CLASS at SITE level --------------------------------------

  summary_site_sizeclass <- summary_transect_sizeclass %>%
    gather(., "TRANSECT", "value", -c(REGION:PHASE)) %>%
    #calculate density by species and sizeclass
    Rmisc::summarySE(
      .,
      measurevar = "value",
      groupvars = c("OCC_SITEID", "OCC_SITENAME", "SIZE_CLASS", "PHASE"),
      na.rm = T
    ) %>%
    # Order rows by Chlorurus, Scarus, than Other
    arrange(match(PHASE, c("J", "I", "T"))) %>%
    # Reorder factors by order displayed in dataframe
    mutate(PHASE = factor(PHASE, levels = c("J", "I", "T")))


  # Summarize SPECIES at TRANSECT LEVEL -------------------------------------

  summary_transect_species <- data %>%
    select(CRUISE_ID:SIZE_CLASS, all_of(metric)) %>%
    ungroup() %>%
    dplyr::group_by(
      REGION,
      REGIONCODE,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      OCC_SITENAME,
      LATITUDE,
      LONGITUDE,
      CB_METHOD,
      TRANSECT,
      SPECIES) %>%
    # calculate metric by species
    dplyr::reframe('SPECIES_SUM' = sum(!!sym(metric)) / (AREA_M2 / 10000)) %>%
    distinct() %>%
    # organize dataframe by Transect and species
    pivot_wider(names_from = TRANSECT, values_from = SPECIES_SUM, names_prefix = paste0("TRANSECT_")) %>%
    mutate_at(-c(1:11), ~ replace_na(., 0))

  # Summarize SPECIES at SITE LEVEL -------------------------------------

  summary_site_species <- data %>%
    select(CRUISE_ID:SIZE_CLASS, all_of(metric)) %>%
        dplyr::group_by(
          REGION,
          REGIONCODE,
          CRUISE_ID,
          LOCATION,
          LOCATIONCODE,
          OCC_SITEID,
          OCC_SITENAME,
          LATITUDE,
          LONGITUDE,
          CB_METHOD
    ) %>%
    reframe(AREA_AVG = mean(AREA_M2)) %>%
    right_join(., data, by = colnames(.)[colnames(.) %in% colnames(data)]) %>%
    # Create Standard Deviation table (put together all the transects)
    Rmisc::summarySE(
      .,
      measurevar = metric,
      groupvars = c(
        "REGION",
        "REGIONCODE",
        "CRUISE_ID",
        "LOCATION",
        "LOCATIONCODE",
        "OCC_SITEID",
        "OCC_SITENAME",
        "LATITUDE",
        "LONGITUDE",
        "CB_METHOD",
        "SPECIES",
        "AREA_AVG",
        "PHASE",
        "SIZE_CLASS"
      ),
      na.rm = T
    ) %>%
    replace(is.na(.), 0) %>%
    dplyr::group_by(
      REGION,
      REGIONCODE,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      OCC_SITENAME,
      LATITUDE,
      LONGITUDE,
      CB_METHOD,
      SPECIES) %>%
    dplyr::reframe(
      'value' = (sum(!!sym(metric)) / (AREA_AVG / 10000)),
      #calculate average density per sp
      'sd' = (sqrt(sum(sd ^
                         2)) / AREA_AVG),
      # calculate standard deviation per sp
      'se' = sd ^ 2
    ) %>%
    distinct() %>% # eliminate duplicate rows
    # Attach Grazing Types Classifications for Excavator, Scraper, etc.
    left_join(., fish_grazing_types, by = "SPECIES") %>%
    # Order rows by Chlorurus, Scarus, then Other
    arrange(match(GENUS, c(
      "Chlorurus", "Scarus", "Calotomus", "Parrotfish"
    )))

  # Summarize OVERALL at TRANSECT LEVEL -------------------------------------

  summary_transect <- summary_transect_sizeclass %>%
    gather(., "TRANSECT", "value",-c(REGION:PHASE)) %>%
    dplyr::group_by(
      REGION,
      REGIONCODE,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      OCC_SITENAME,
      LATITUDE,
      LONGITUDE,
      CB_METHOD,
      TRANSECT) %>%
    dplyr::reframe('TRANSECT_SUM' = sum(value)) %>%
    spread(., TRANSECT, "TRANSECT_SUM")

  # Summarize OVERALL at SITE LEVEL -------------------------------------

  summary_site <- summary_transect_sizeclass %>%
    gather(., "TRANSECT", "value",-c(REGION:PHASE)) %>%
    dplyr::group_by(
      REGION,
      REGIONCODE,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      OCC_SITENAME,
      LATITUDE,
      LONGITUDE,
      CB_METHOD,
      TRANSECT
    ) %>%
    reframe('TRANSECT_SUM' = sum(value)) %>%
    dplyr::group_by(
      REGION,
      REGIONCODE,
      CRUISE_ID,
      LOCATION,
      LOCATIONCODE,
      OCC_SITEID,
      OCC_SITENAME,
      LATITUDE,
      LONGITUDE,
      CB_METHOD
    ) %>%
    reframe(
      TOTAL = sum(TRANSECT_SUM),
      MEAN = mean(TRANSECT_SUM, na.rm = TRUE),
      SD = sd(TRANSECT_SUM),
      SE = sd(TRANSECT_SUM) / sqrt(length(TRANSECT_SUM)),
      N = length(TRANSECT_SUM)
    )


  if (summarize_by == "size class" && level == "transect") {
    return(summary_transect_sizeclass)
  }
  if (summarize_by == "species" && level == "transect") {
    return(summary_transect_species)
  }
  if (summarize_by == "overall" && level == "transect") {
    return(summary_transect)
  }
  if (summarize_by == "size class" && level == "site") {
    return(summary_site_sizeclass)
  }
  if (summarize_by == "species" && level == "site") {
    return(summary_site_species)
  }
  if (summarize_by == "overall" && level == "site") {
    return(summary_site)
  }

}

