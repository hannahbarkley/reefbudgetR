#' Create fish associated SPC sites in relation to fixed sites
#'
#'@author Thomas Oliver 
#'@author Kisei Tanaka
#'@author Rebecca Weible
#'
#'@param data fish survey data.
#'@param subset_distance_m Assigned associated site distances, in meters, from fixed SPC/OCC site to all other fish SPC sites, based on parrotfish foraging boundaries. 

#'
#'@import Rmisc
#'@import terra
#'@import sf
#'@import tidyverse
#'@import dplyr
#'
#'@export create_fish_assoc_sites
#'
#'@examples
#'fish_data <- read.csv("ESD_CarbBudget_SPC_OAHU_2021.csv", na = "", check.names = FALSE)
#'
#'sites_associated_dbase <- create_fish_assoc_sites(data = fish_data, subset_distance_m = 2000)


create_fish_assoc_sites <- function(data, fixed_metadata, subset_distance_m) {
  
  #  Download/Load Shapefile -----------------------
  if ("American Samoa" %in% unique(data$REGION)) {
    if (!exists("shape_file"))
      stop("'shape_file' not found in environment for American Samoa.")
    shape_file <- shape_file
    
  } else {
    dl_dir <- tempdir()
    gshhg_path <- file.path(dl_dir, "gshhg-bin-2.3.7/GSHHS_shp/f/GSHHS_f_L1.shp")
    
    if (!file.exists(gshhg_path)) {
      message("Downloading GSHHG shapefiles...") 
      zip_path <- file.path(dl_dir, "gshhg.zip")
      
      op <- getOption("timeout")
      options(timeout = max(600, getOption("timeout")))
      
      tryCatch({
        download.file(
          "https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.7.zip",
          destfile = zip_path,
          mode = "wb"
        )
        unzip(zip_path, exdir = dl_dir)
      }, finally = {
        options(timeout = op)
        if (file.exists(zip_path))
          unlink(zip_path)
      })
    }
    shape_file <- st_read(gshhg_path, quiet = TRUE)
  }
  
  # Format Data -------------------------------------------------------------
  
  sf_use_s2(FALSE)
  
  # Don't clean the whole world here (too slow). We clean just the crop later.
  ALLisl_poly <- shape_file
  
  if (!"OCC_SITEID" %in% names(data))
    data$OCC_SITEID <- NA
  
  # Filter Fish Data
  ptsFISH <- data %>%
    filter(METHOD == "nSPC") %>%
    dplyr::select(
      REGION:REA_SITEID,
      LATITUDE:LONGITUDE,
      METHOD,
      HABITAT_CODE:HABITAT_TYPE,
      OCC_SITEID
    ) %>%
    distinct() %>%
    st_as_sf(coords = c("LONGITUDE", "LATITUDE"),
             crs = st_crs(ALLisl_poly))
  
  # Filter Fixed Data
  ptsOCC <- fixed_metadata %>%
    filter(LOCATIONCODE %in% ptsFISH$LOCATIONCODE) %>%
    dplyr::select(-any_of(c("LOCALDATE", "SITE_DEPTH_M"))) %>%
    mutate(REA_SITEID = "NA", METHOD = "NA") %>%
    st_as_sf(coords = c("LONGITUDE", "LATITUDE"),
             crs = st_crs(ALLisl_poly))
  
  # Calculate Distances ----------------------------------------
  uLC <- unique(ptsOCC$LOCATIONCODE)
  
  dist_list <- lapply(uLC, function(thisisl) {
    ptsFISH_ISL <- ptsFISH %>% filter(LOCATIONCODE == thisisl)
    ptsOCC_ISL <- ptsOCC %>% filter(LOCATIONCODE == thisisl)
    
    if (nrow(ptsFISH_ISL) == 0 ||
        nrow(ptsOCC_ISL) == 0)
      return(NULL)
    
    pts_ext <- terra::ext(bind_rows(ptsFISH_ISL, ptsOCC_ISL)) * 1.5
    
    suppressWarnings({
      poly_ISL <- st_crop(x = ALLisl_poly, y = st_bbox(pts_ext))
      poly_ISL <- st_make_valid(poly_ISL)
    })
    # ------------------------------------------------------------------
    
    D <- calc_fish_distance_matrix(
      pointsOCC = ptsOCC_ISL,
      pointsFISH = ptsFISH_ISL,
      island_poly = poly_ISL,
      resolution_m = 50
    )
    
    as.data.frame(as.table(D)) %>%
      dplyr::rename(OCCSITE = Var1,
                    REASITE = Var2,
                    DISTANCE_m = Freq)
  })
  
  out <- bind_rows(dist_list)
  
  # Create Associations -----------------------------------------------------
  subset_by_dist <- out %>%
    filter(DISTANCE_m <= subset_distance_m) %>% # Optimization: Filter early
    mutate(OCCSITE = as.character(OCCSITE), REASITE = as.character(REASITE)) %>%
    dplyr::select(OCCSITE, REASITE)
  
  # Prepare tables for Join
  fixed_dat <- fixed_metadata %>%
    filter(LOCATIONCODE %in% ptsFISH$LOCATIONCODE) %>%
    dplyr::select(-any_of(c("LOCALDATE", "SITE_DEPTH_M"))) %>%
    mutate(
      REA_SITEID = "NA",
      METHOD = "NA",
      value = 1,
      ASSOC_OCCSITEID = OCC_SITEID
    )
  
  fish_dat <- data %>%
    dplyr::select(REGION:LONGITUDE, METHOD, HABITAT_CODE:HABITAT_TYPE) %>%
    distinct()
  
  fish_expanded <- fish_dat %>%
    left_join(subset_by_dist, by = c("REA_SITEID" = "REASITE")) %>%
    dplyr::rename(ASSOC_OCCSITEID = OCCSITE) %>%
    mutate(value = ifelse(!is.na(ASSOC_OCCSITEID), -1, 0))
  
  finaldf <- bind_rows(fixed_dat, fish_expanded)
  
  # Apply Filters ------------------------------------------------
  non_forage <- c("AGR", "APR", "APS", "PPR", "ROB", "SAG", "WAL")
  
  finaldf <- finaldf %>%
    mutate(ASSOC_OCCSITEID = gsub(" ", "", ASSOC_OCCSITEID)) %>%
    mutate(value = case_when(value == -1 &
                               HABITAT_CODE %in% non_forage ~ 0, TRUE ~ value))
  
  fixedhab <- finaldf %>%
    filter(value == 1) %>%
    dplyr::select(OCC_SITEID, FIXED_HAB = HABITAT_CODE) %>%
    distinct()
  
  output <- finaldf %>%
    left_join(fixedhab, by = c("ASSOC_OCCSITEID" = "OCC_SITEID")) %>%
    mutate(value = case_when(value == -1 &
                               HABITAT_CODE != FIXED_HAB ~ 0, TRUE ~ value)) %>%
    dplyr::select(-FIXED_HAB)
  
  # Sample Size Summary -----------------------------------------------------
  strs_samplesize <- output %>%
    dplyr::select(ASSOC_OCCSITEID, HABITAT_CODE, value) %>%
    group_by(ASSOC_OCCSITEID, value) %>%
    count() %>%
    ungroup() %>%
    mutate(
      value_label = case_when(
        value == -1 ~ "Associated",
        value == 0  ~ "Not_Associated",
        value == 1  ~ "Fixed"
      )
    ) %>%
    dplyr::select(-value) %>%
    tidyr::spread(value_label, n, fill = 0)
  
  return(list(output = output, surveysamplesize = strs_samplesize))
}
  