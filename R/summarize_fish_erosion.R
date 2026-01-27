#' Summarize site-level parrotfish bioerosion rates
#'
#' @author Rebecca Weible
#'
#' @param species_table Species table with average density, biomass, and bioerosion metrics,
#' outputs of `summarize_fish_metrics`.
#' @param full_summary Return full summary of results (default is TRUE)
#'
#' @import stringr
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


summarize_fish_erosion <- function(species_table, full_summary = TRUE) {
  
  options(scipen = 999) # prevent scientific notation
  
  # Base Summary (Long Format) ----------------------------------------------
  
  base_long <- species_table %>%
    dplyr::select(-SPECIES) %>%
    tidyr::pivot_longer(
      cols = starts_with("TRANSECT_"),
      names_to = "TRANSECT",
      values_to = "VALUE"
    )
  
  # Calculate Sums for each Functional Group
  sum_by_fxn <- base_long %>%
    dplyr::group_by(
      REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, 
      OCC_SITEID, LATITUDE, LONGITUDE, METHOD, METRIC, TRANSECT, FXN_GRP
    ) %>%
    dplyr::summarise(VALUE = sum(VALUE, na.rm = TRUE), .groups = "drop")
  
  # Calculate Sums for "All" Groups combined
  sum_all <- base_long %>%
    dplyr::group_by(
      REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, 
      OCC_SITEID, LATITUDE, LONGITUDE, METHOD, METRIC, TRANSECT
    ) %>%
    dplyr::summarise(VALUE = sum(VALUE, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(FXN_GRP = "All")
  
  # Combine Specific Groups + "All"
  combined_long <- dplyr::bind_rows(sum_by_fxn, sum_all)
  
  
  # Transect Level Summary ----------------------------------------
  
  fish_erosion_transect <- combined_long %>%
    tidyr::pivot_wider(
      names_from = METRIC,
      values_from = VALUE,
      values_fill = 0
    ) %>%
    dplyr::mutate(
      TRANSECT_NUM = as.integer(sub("TRANSECT_", "", TRANSECT))
    ) %>%
    dplyr::arrange(REGION, LOCATION, METHOD, TRANSECT_NUM, FXN_GRP) %>%
    dplyr::select(-TRANSECT_NUM) %>%
    dplyr::select(
      REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, 
      OCC_SITEID, LATITUDE, LONGITUDE, METHOD, TRANSECT, FXN_GRP,
      everything()
    )
  
  
  # Site Level Summary --------------------------------------------
  
  fish_erosion_site_long <- combined_long %>%
    dplyr::group_by(
      REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, 
      OCC_SITEID, LATITUDE, LONGITUDE, METHOD, METRIC, FXN_GRP
    ) %>%
    dplyr::summarise(
      MEAN = mean(VALUE, na.rm = TRUE),
      SD   = sd(VALUE, na.rm = TRUE),
      n    = n(), # Should be 10, but safe to calculate
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      SD  = dplyr::coalesce(SD, 0),
      SE  = SD / sqrt(10), # Assuming n=10 fixed as per original code
      L95 = pmax(0, MEAN - (SE * 1.97)), # Vectorized 'max' logic
      U95 = MEAN + (SE * 1.97)
    )
  
  # Format for NCEI 
  fish_erosion_site <- fish_erosion_site_long %>%
    tidyr::pivot_longer(
      cols = c(MEAN, SD, SE, L95, U95),
      names_to = "STAT",
      values_to = "VAL"
    ) %>%
    tidyr::unite("COL_NAME", c(METRIC, FXN_GRP, STAT), sep = "_") %>%
    tidyr::pivot_wider(
      names_from = COL_NAME,
      values_from = VAL,
      values_fill = 0
    ) 
  
  # Ensure uppercase names
  names(fish_erosion_site) <- toupper(names(fish_erosion_site))
  
  # Reorder columns to ensure metadata is first
  fish_erosion_site <- fish_erosion_site %>%
    dplyr::select(
      REGION, REGIONCODE, CRUISE_ID, LOCATION, LOCATIONCODE, 
      OCC_SITEID, LATITUDE, LONGITUDE, METHOD,
      everything()
    )
  
  # Return ------------------------------------------------------------------
  if (full_summary) {
    return(list(
      fish_erosion_transect = fish_erosion_transect,
      fish_erosion_site = fish_erosion_site
    ))
  } else {
    return(fish_erosion_site)
  }
}
