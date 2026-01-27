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


summarize_fish_metrics <- function(data,
                                   metric = c("density", "biomass", "bioerosion"),
                                   level = c("transect", "site"),
                                   summarize_by = c("size class", "species", "overall")) {
  
  # Setup Parameters --------------------------------------------------------
  metric_arg <- match.arg(metric)
  level_arg  <- match.arg(level)
  sum_by_arg <- match.arg(summarize_by)
  
  # Map generic metric names to column names
  col_name <- switch(metric_arg,
                     "density"    = "COUNT",
                     "biomass"    = "BIOMASS_PER_FISH_KG",
                     "bioerosion" = "FISH_EROSION_KG_M2_YR")
  
  # Define Area Normalization Factor
  # Density/Biomass = per Hectare (Area / 10,000)
  # Bioerosion = per m2 (Area / 1)
  area_denom_factor <- if (metric_arg == "bioerosion") 1 else 10000
  
  # Define Grouping Columns based on user input
  grp_cols <- switch(sum_by_arg,
                     "size class" = c("SIZE_CLASS", "PHASE"),
                     "species"    = c("FXN_GRP", "SPECIES"),
                     "overall"    = character(0)) 
  
  # Define Metadata Columns (to keep in output)
  meta_cols <- c("REGION", "REGIONCODE", "CRUISE_ID", "LOCATION", 
                 "LOCATIONCODE", "OCC_SITEID", "LATITUDE", "LONGITUDE", "METHOD")
  
  
  # Base Calculation (Transect Level) ---------------------------------------
  # Calculate the metric per transect first. 
  
  # Aggregate
  base_transect <- data %>%
    dplyr::group_by(dplyr::across(all_of(c(meta_cols, "TRANSECT", grp_cols)))) %>%
    dplyr::summarise(
      # Sum the metric for the group, then divide by Area (converted if needed)
      # We take the mean Area per transect group (should be unique per transect)
      VAL = sum(!!sym(col_name), na.rm = TRUE) / (mean(AREA_M2) / area_denom_factor),
      .groups = "drop"
    )
  
  # Zero Filling 
  
  # Define standard transect list
  transect_levels <- paste0("TRANSECT_", 1:10)
  
  base_filled <- base_transect %>%
    dplyr::mutate(TRANSECT = factor(paste0("TRANSECT_", TRANSECT), levels = transect_levels)) %>%
    tidyr::complete(
      tidyr::nesting(!!!syms(c(meta_cols, grp_cols))), # Keep existing metadata/groups together
      TRANSECT, 
      fill = list(VAL = 0)
    )
  
  
  # Output Generation -----------------------------------------------
  
  # --- Transect Level Output  ---
  if (level_arg == "transect") {
    
    output <- base_filled %>%
      tidyr::pivot_wider(
        names_from = TRANSECT,
        values_from = VAL,
        values_fill = 0
      ) %>%
      dplyr::select(all_of(meta_cols), any_of(grp_cols), everything())
    
    return(output)
  }
  
  # --- Site Level Output  ---
  if (level_arg == "site") {
    
    output <- base_filled %>%
      dplyr::group_by(dplyr::across(all_of(c(meta_cols, grp_cols)))) %>%
      dplyr::summarise(
        MEAN = mean(VAL, na.rm = TRUE),
        SD   = sd(VAL, na.rm = TRUE),
        SE   = SD / sqrt(n()), # n() will be 10 due to complete() above
        N    = n(),
        .groups = "drop"
      )
    
    # If summarizing by SIZE CLASS: Sort by Phase
    if (sum_by_arg == "size class") {
      output <- output %>%
        dplyr::mutate(PHASE = factor(PHASE, levels = c("J", "I", "T"))) %>%
        dplyr::arrange(PHASE)
    }
    
    # If summarizing by SPECIES: Attach Genus/Grazing info
    if (sum_by_arg == "species") {
      # Assumes 'fish_grazing_types' exists globally, as per original code
      if (exists("fish_grazing_types")) {
        output <- output %>%
          dplyr::left_join(
            fish_grazing_types %>% dplyr::select(SPECIES, GENUS), 
            by = "SPECIES"
          ) %>%
          dplyr::arrange(match(GENUS, c("Chlorurus", "Scarus", "Calotomus", "Parrotfish"))) %>%
          dplyr::select(all_of(meta_cols), any_of(grp_cols), GENUS, everything())
      }
    }
    
    return(output)
  }
}