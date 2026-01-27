#' Calculate the shortest distance between each fixed SPC/OCC site and each fish REA site 'as a fish swims'
#'
#'@author Thomas Oliver
#'@author Rebecca Weible
#'
#'@param pointsOCC Lat/Long of fixed OCC site, or fixed SPC data, aka target sites.
#'@param pointsFISH Lat/Long of fish REA sites, or stratified random sampled SPC surveys, aka all other survey sites.
#'@param island_poly Simple Feature polygon Shapefile where surveys were conducted.
#'@param resolution_m Resolution of the output raster plot displaying distance from target point to survey point of interest.

#'
#'@import terra
#'@import sf
#'@import tidyverse
#'@import dplyr
#'
#'@export calc_fish_distance_matrix
#'
#'@examples
#'distance_matrix <- calc_fish_distance_matrix(pointsOCC = ptsOCC_ISL, pointsFISH = ptsFISH_ISL, island_poly = poly_ISL, resolution_m = 50)


calc_fish_distance_matrix <- function(pointsOCC,
                                      pointsFISH,
                                      island_poly,
                                      resolution_m = 50) {
  # 1. Setup & Pre-processing --------------------------------------------------
  # Calculate resolution in degrees
  res_d <- resolution_m / 111111
  
  # Convert inputs to terra SpatVectors once (much faster than sf for raster ops)
  v_occ <- terra::vect(pointsOCC)
  v_fish <- terra::vect(pointsFISH)
  v_island <- terra::vect(island_poly)
  
  # Create raster covering combined extent
  # c() combines the vectors to get total extent, * 1.2 adds buffer
  e <- terra::ext(c(v_occ, v_fish)) * 1.2
  r <- terra::rast(e, res = res_d, crs = terra::crs(v_occ))
  
  # 2. Rasterize Island (Land = NA, Water = 99) --------------------------------
  # Initialize raster with 99 (Water)
  island_rast <- terra::init(r, 99)
  
  # Burn island polygon into raster as NA (Land/Barrier)
  # 'touches=TRUE' ensures small islands are caught
  island_rast <- terra::rasterize(
    v_island,
    island_rast,
    field = NA,
    update = TRUE,
    touches = TRUE
  )
  
  # 3. Force Points to be Water ------------------------------------------------
  # Get cell indices for all points once (Fastest access method)
  cells_occ <- terra::cellFromXY(island_rast, terra::crds(v_occ))
  cells_fish <- terra::cellFromXY(island_rast, terra::crds(v_fish))
  
  # Force these specific cells to be traversable (99)
  # This handles cases where a point might slightly overlap a land pixel due to resolution
  island_rast[cells_occ] <- 99
  island_rast[cells_fish] <- 99
  
  # Initialize Output Matrix
  dist_matrix <- matrix(0, nrow = length(v_occ), ncol = length(v_fish))
  rownames(dist_matrix) <- pointsOCC$OCC_SITEID
  colnames(dist_matrix) <- pointsFISH$REA_SITEID
  
  # 4. Calculate Distances (Optimized Loop) ------------------------------------
  
  # Loop through Source (OCC) points
  for (i in 1:length(cells_occ)) {
    target_idx <- cells_occ[i]
    
    # Skip if point falls outside raster extent (rare safety check)
    if (is.na(target_idx))
      next
    
    # Modify raster in-place: Set Target cell to 0
    # (gridDist calculates distance TO the nearest cell with value 0)
    prev_val <- island_rast[target_idx] # Store previous value (99)
    island_rast[target_idx] <- 0
    
    # Calculate accumulated distance grid
    # gridDist respects NAs (land) as barriers
    dist_rast <- terra::gridDist(island_rast, target = 0)
    
    # --- PLOTTING (Bottleneck - commented out for efficiency) ---
    # plot(dist_rast)
    # plot(island_poly$geometry, add = TRUE, col = "white")
    # plot(pointsFISH$geometry[1], add=TRUE, col="blue", pch=19)
    # plot(pointsOCC$geometry[i], add=TRUE, col="red", pch=19)
    # title(main= paste(pointsOCC$OCC_SITEID[i]))
    # ------------------------------------------------------------
    
    # Extract distances for ALL fish points at once (Vectorized)
    # Direct vector subsetting [cells_fish] is instantaneous compared to the original 'extract' loop
    vals <- dist_rast[cells_fish]
    
    # Assign row to matrix
    # vals is a data.frame, extract first column
    dist_matrix[i, ] <- as.numeric(vals[, 1])
    
    # Revert the target cell back to original value for the next iteration
    # This avoids copying the whole raster object every loop
    island_rast[target_idx] <- prev_val
  }
  
  return(dist_matrix)
}