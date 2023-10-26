#' Calculate the shortest distance between each fixed SPC/OCC site and each fish REA site 'as a fish swims'
#'
#'@author Thomas Oliver 
#'@author Rebecca Weible
#'
#'@param pointsOCC Lat/Long of fixed OCC site, or fixed SPC data, aka target sites.
#'@param pointsFISH Lat/Long of fish REA sites, or stratified random sampled SPC surveys, aka all other survey sites.
#'@param island_poly Simple Feature polygon Shapefile of Pacific Islands where surveys were conducted.
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




calc_fish_distance_matrix <- function(pointsOCC, pointsFISH, island_poly, resolution_m=50) {
  
  res_d=resolution_m/111111
  # Create a raster covering the combined extent of the island and points
  all_features <- rbind(pointsOCC, pointsFISH)
  r <- rast(ext(all_features)*1.2, res=res_d)
  
  # Rasterize the island
  island_rast <- terra::rasterize(x=vect(island_poly), y=r, field=1)
  #set water to 99
  island_rast[is.na(island_rast)]=99
  #set land to NA
  island_rast[island_rast==1]=NA
  
  #also force any points of interest to not be NA
  island_rast[terra::extract(island_rast, vect(pointsOCC), cells=TRUE)$cell]=99
  island_rast[terra::extract(island_rast, vect(pointsFISH), cells=TRUE)$cell]=99
  
  # Placeholder matrix for distances
  dist_matrix <- matrix(0, nrow=nrow(pointsOCC), ncol=nrow(pointsFISH))
  rownames(dist_matrix) <- pointsOCC$OCC_SITEID
  colnames(dist_matrix) <- pointsFISH$REA_SITEID
  
  # For each OCC point (fewer of them), calculate the grid distance avoiding the island
  for (i in 1:nrow(pointsOCC)) {
    target_point <- pointsOCC[i,]
    # Identify the cell in the raster where the point lies
    target_idx <- terra::extract(island_rast, vect(target_point), cells=TRUE)$cell

    # Set the value of that cell to a desired value, e.g., 999
    target_rast <- island_rast
    target_rast[target_idx] <- 0
    dist_rast <- terra::gridDist(x = target_rast)

   
     #plot(target_rast)
     plot(dist_rast)
     plot(island_poly$geometry, add = TRUE, col = "white")
     plot(pointsFISH[1], add=TRUE, col="blue",pch=19)
     plot(target_point[1], add=TRUE, col="red",pch=19)

    # For each point in set B, extract the distance from the current distance grid
    for (j in 1:nrow(pointsFISH)) {
      #j=36
      dist_matrix[i, j] <- as.numeric(terra::extract(dist_rast, pointsFISH[j,])$layer)
    }
  }
  
  return(dist_matrix)
}
