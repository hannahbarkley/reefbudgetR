#' Calculate parrotfish erosion rates from belt data
#'
#'@author Rebecca Weible
#'
#'@param data fish survey data.
#'@param rates_dbase Erosion rates database to use.
#'
#'@import Rmisc
#'@import tidyverse
#'@import tibble
#'@import dplyr
#'
#'@export calc_fish_belt
#'
#'@examples
#'fish_data <- read.csv("CB_FishBelt_alldata.csv", na = "", check.names = FALSE)
#'
#'fish_belt <- calc_fish_belt(data = fish_data, rates_dbase = rates_dbase,
#'full_summary = TRUE)

calc_fish_belt <- function(data, 
                           rates_dbase, 
                           fixed_metadata = sites_metadata,
                           full_summary = TRUE) {
  
  # 0. Clean duplicate columns ----------------------------------------------
  # Safely removes duplicate columns (like "TRANSECT") that cause pivot_longer to crash
  data <- data[, !duplicated(colnames(data)), drop = FALSE]
  
  # 1. Pre-process data -----------------------------------------------------
  data_prepped <- prep_spc(data, fixed_metadata)
  calc_eros_fish_output <- calc_eros_fish(data_prepped, rates_dbase)
  
  # 2. Calculate metrics ----------------------------------------------------
  metrics_map <- list(
    density    = "FISH_DENSITY_ABUNDANCE_HA",
    biomass    = "FISH_BIOMASS_KG_HA",
    bioerosion = "FISH_EROSION_KG_M2_YR"
  )
  
  # Iterate over the map and bind results into one dataframe
  species_table <- purrr::imap_dfr(metrics_map, function(label, metric_key) {
    suppressWarnings(
      summarize_fish_metrics(
        data = calc_eros_fish_output,
        metric = metric_key,
        level = "transect",
        summarize_by = "species"
      )
    ) %>%
      dplyr::mutate(METRIC = label)
  })
  
  # 3. Add back zero-count transects safely ---------------------------------
  # Identify the core grouping columns that exist in the calculated species table
  id_cols <- intersect(c("OCC_SITEID", "SITE", "CB_TRANSECTID", "TRANSECT", "TRANSECT_ID", "REPLICATE"), names(species_table))
  
  # Extract the true master roster of all transects from the prepped data
  base_transects <- calc_eros_fish_output %>%
    dplyr::select(dplyr::all_of(id_cols)) %>%
    dplyr::distinct()
  
  # Find exactly which transects were dropped (the empty ones)
  missing_transects <- dplyr::anti_join(base_transects, species_table, by = id_cols)
  
  if (nrow(missing_transects) > 0) {
    # Create blank rows for the missing transects crossed with our 3 metrics
    zero_rows <- tidyr::expand_grid(
      missing_transects,
      METRIC = unique(species_table$METRIC)
    )
    
    # Safely assign "NONE" to whichever species columns your package happens to output
    sp_cols <- intersect(c("TAXON_CODE", "TAXON_NAME", "SPECIES", "SPECIES_NAME"), names(species_table))
    for (sc in sp_cols) {
      zero_rows[[sc]] <- "NONE"
    }
    
    # Bind the empty transects to the bottom of the table and fill metrics with 0
    species_table <- dplyr::bind_rows(species_table, zero_rows) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ tidyr::replace_na(., 0)))
  }
  
  # 4. Final Summary --------------------------------------------------------
  summary_belt_erosion <- summarize_fish_erosion(species_table, full_summary)
  
  # 5. Return Logic ---------------------------------------------------------
  if (full_summary) {
    return(list(
      fish_erosion_transect = summary_belt_erosion$fish_erosion_transect,
      fish_erosion_site     = summary_belt_erosion$fish_erosion_site
    ))
  } else {
    return(summary_belt_erosion$fish_erosion_site)
  }
}