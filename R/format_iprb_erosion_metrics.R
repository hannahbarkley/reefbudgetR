#' Format Indo-Pacific ReefBudget parrotfish erosion rates
#'
#'@author Rebecca Weible
#'
#'@param iprb_rates Indo-Pacific ReefBudget parrotfish erosion rates
#'("Equations" tab in IP_Parrotfish_erosion_rates_database_v1.3.xlsx,
#'downloaded from https://geography.exeter.ac.uk/reefbudget/indopacific/).
#'
#'@import dplyr
#'@importFrom sjmisc seq_row
#'
#'@export format_iprb_erosion_metrics
#'

format_iprb_erosion_metrics <- function(iprb_rates, output = c("all", "finalerosion")) {
  
  output <- match.arg(output)
  
  rates_raw <- iprb_rates %>%
    dplyr::filter(X.2 != "") %>%
    dplyr::rename(
      GrazType = "X.1",
      TAXONNAME = "X.2",
      `I.11-20cm` = "Initial.Phase",
      `I.21-30cm` = "X.3",
      `I.31-40cm` = "X.4",
      `I.41-50cm` = "X.5",
      `T.11-20cm` = "Terminal.Phase",
      `T.21-30cm` = "X.6",
      `T.31-40cm` = "X.7",
      `T.41-50cm` = "X.8",
      `T.51-60cm` = "X.9",
      I.x20 = "X.11", I.x30 = "X.12", I.x40 = "X.13", I.x50 = "X.14",
      T.x20 = "X.15", T.x30 = "X.16", T.x40 = "X.17", T.x50 = "X.18", T.x60 = "X.19"
    )
  

  split_indices <- cumsum(1:nrow(rates_raw) %in% c(26, 45, 65, 84))
  rates_list <- split(rates_raw, split_indices)
  names(rates_list) <- c("SizeClassErosionRates", "PropBiteScars", "BiteRates", "VolumeRemoved", "MassRemoved")
  
  rates_processed <- lapply(rates_list, function(df) {
    df %>%
      dplyr::filter(TAXONNAME != "SPECIES") %>%
      dplyr::select(-any_of(c("X", "X.10"))) %>% 
      dplyr::mutate(across(-c(GrazType, TAXONNAME), ~if_else(. == "", "0", .))) %>%
      # Pivot Size Bins Long
      tidyr::pivot_longer(
        cols = -c(GrazType, TAXONNAME),
        names_to = "size_raw",
        values_to = "count"
      ) %>%
      tidyr::separate(size_raw, c("PHASE", "sizebin"), extra = "merge", fill = "left") %>%
      dplyr::filter(PHASE != "X") %>%
      # Pivot wide (Sizebin -> Count)
      tidyr::pivot_wider(names_from = sizebin, values_from = count, values_fill = "0") %>%
      dplyr::select(-matches("^x[0-9]+"))
  })
  
  if (output == "finalerosion") {
    return(
      rates_processed[["SizeClassErosionRates"]] %>%
        tidyr::pivot_longer(cols = matches("cm"), names_to = "sizeclass", values_to = ".value") %>%
        dplyr::select(PHASE, sizeclass, TAXONNAME, .value)
    )
  }
  
  species_map <- c(
    "ocellatus" = "Cetoscarus ocellatus",
    "microrhinos" = "Chlorurus microrhinos",
    "gibbus" = "Chlorurus gibbus",
    "spilurus" = "Chlorurus spilurus",
    "bleekeri" = "Chlorurus bleekeri",
    "capistratoides" = "Chlorurus capistratoides",
    "japanensis" = "Chlorurus japanensis",
    "longiceps" = "Hipposcarus longiceps",
    "falcipinnis" = "Scarus falcipinnis",
    "altipinnis" = "Scarus altipinnis",
    "persicus" = "Scarus persicus",
    "forsteni" = "Scarus forsteni",
    "flavipectoralis" = "Scarus flavipectoralis",
    "fuscopurpureus" = "Scarus fuscopurpureus",
    "russelii" = "Scarus russelii",
    "rivulatus" = "Scarus rivulatus",
    "chameleon" = "Scarus chameleon",
    "festivus" = "Scarus festivus",
    "oviceps" = "Scarus oviceps",
    "dimidiatus" = "Scarus dimidiatus",
    "viridifucatus" = "Scarus viridifucatus",
    "xanthopleura" = "Scarus xanthopleura"
  )
  
  output_metrics <- dplyr::bind_rows(rates_processed, .id = "metric") %>%
    dplyr::filter(!metric %in% c("SizeClassErosionRates", "MassRemoved")) %>%
    dplyr::mutate(
      metric = dplyr::case_match(metric,
                                 "PropBiteScars" ~ "Proportion of Scars",
                                 "BiteRates" ~ "Bite Rate",
                                 "VolumeRemoved" ~ "Bite Volume"
      ),
      FISH_sciname = strsplit(TAXONNAME, "/")
    ) %>%
    tidyr::unnest(FISH_sciname) %>%
    dplyr::mutate(
      FISH_sciname = ifelse(FISH_sciname %in% names(species_map), species_map[FISH_sciname], FISH_sciname),
      type = "CPM"
    ) %>%
    dplyr::mutate(across(matches("[0-9]+-[0-9]+cm"), ~round(as.numeric(.), 5))) %>%
    dplyr::select(-any_of("10")) %>%
    # Gather size bins
    tidyr::pivot_longer(
      cols = matches("cm$"), 
      names_to = "sizeclass", 
      values_to = "value"
    ) %>%
    tidyr::unite("sizeclass_combined", PHASE, sizeclass, sep = ".") %>%
    dplyr::filter(FISH_sciname != "Substrate density (g cm-3)") %>%
    tidyr::pivot_wider(names_from = metric, values_from = value, values_fill = 0) %>%
    dplyr::right_join(
      GrazTypes %>% dplyr::select(GRAZ_TYPE, SPECIES, sister_sp, FISH_sciname),
      by = "FISH_sciname"
    ) %>%
    tidyr::complete(sizeclass_combined, type, nesting(FISH_sciname, GRAZ_TYPE, SPECIES, sister_sp)) %>%
    dplyr::filter(!is.na(sizeclass_combined)) %>%
    dplyr::mutate(across(c(`Bite Rate`, `Bite Volume`, `Proportion of Scars`), ~tidyr::replace_na(., 0))) %>%
    dplyr::arrange(FISH_sciname) %>%
    dplyr::filter(!is.na(type)) %>%
    dplyr::select(SPECIES, FISH_sciname, sizeclass = sizeclass_combined, everything(), -sister_sp, -type, -GRAZ_TYPE, -TAXONNAME)
  
  return(output_metrics)
}