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
#'@import tidyr
#'
#'@export summarize_fish_metrics

summarize_fish_metrics <- function(data,
                                   metric = c("density", "biomass", "bioerosion"),
                                   level = c("transect", "site"),
                                   summarize_by = c("size class", "species", "overall")) {
  
  # Setup Parameters --------------------------------------------------------
  metric_arg <- match.arg(metric)
  level_arg  <- match.arg(level)
  sum_by_arg = match.arg(summarize_by)
  
  # Map generic metric names to column names
  col_name <- switch(metric_arg,
                     "density"    = "COUNT",
                     "biomass"    = "BIOMASS_PER_FISH_KG",
                     "bioerosion" = "FISH_EROSION_KG_M2_YR")
  
  # Define Area Normalization Factor
  area_denom_factor <- if (metric_arg == "bioerosion") 1 else 10000
  
  # Define Grouping Columns based on user input
  grp_cols <- switch(sum_by_arg,
                     "size class" = c("SIZE_CLASS", "PHASE"),
                     "species"    = c("FXN_GRP", "SPECIES"),
                     "overall"    = character(0)) 
  
  # Define Metadata Columns
  meta_cols <- c("REGION", "REGIONCODE", "CRUISE_ID", "LOCATION", 
                 "LOCATIONCODE", "OCC_SITEID", "LATITUDE", "LONGITUDE", "METHOD")
  
  # Base Calculation (Transect Level) ---------------------------------------
  base_transect <- data %>%
    dplyr::group_by(dplyr::across(all_of(c(meta_cols, "TRANSECT", grp_cols)))) %>%
    dplyr::summarise(
      VAL = sum(!!sym(col_name), na.rm = TRUE) / (mean(AREA_M2) / area_denom_factor),
      .groups = "drop"
    )
  
  # Dynamic Zero Filling Per Site -------------------------------------------
  base_filled <- base_transect %>%
    # FIX: Select only the structural columns needed for expansion.
    # This strips metadata out temporarily so left_join doesn't create .x/.y suffixes.
    dplyr::select(OCC_SITEID, TRANSECT, any_of(grp_cols), VAL) %>%
    dplyr::group_by(OCC_SITEID) %>%
    dplyr::group_modify(~ {
      site_max <- max(as.numeric(as.character(.x$TRANSECT)), na.rm = TRUE)
      site_transect_levels <- paste0("TRANSECT_", 1:site_max)
      
      .x %>%
        dplyr::mutate(TRANSECT = factor(paste0("TRANSECT_", TRANSECT), levels = site_transect_levels)) %>%
        tidyr::complete(
          tidyr::nesting(!!!syms(intersect(names(.x), grp_cols))), 
          TRANSECT,
          fill = list(VAL = 0)
        )
    }) %>%
    dplyr::ungroup() %>%
    # FIX: Cleanly map the full metadata columns back onto every row by OCC_SITEID
    dplyr::left_join(
      base_transect %>% dplyr::select(all_of(meta_cols)) %>% dplyr::distinct(),
      by = "OCC_SITEID"
    )
  
  # Output Generation -----------------------------------------------
  
  # --- Transect Level Output  ---
  if (level_arg == "transect") {
    output <- base_filled %>%
      tidyr::pivot_wider(
        names_from = TRANSECT,
        values_from = VAL
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
        SE   = SD / sqrt(n()), 
        N    = n(),
        .groups = "drop"
      )
    
    if (sum_by_arg == "size class") {
      output <- output %>%
        dplyr::mutate(PHASE = factor(PHASE, levels = c("J", "I", "T"))) %>%
        dplyr::arrange(PHASE)
    }
    
    if (sum_by_arg == "species") {
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