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


create_fish_assoc_sites <- function(data, fixed_metadata, subset_distance_m){
  
  # Download shapefile from GSHHG to users Downloads folder---------------------------------------------------
  if(unique(data$REGION) == "American Samoa") {
    #load("islands_shp.Rdata") #loaded as shape_file
    shape_file <- shape_file
    
    } else{
      
    data_path <- paste0("C:/Users/", Sys.info()[7], "/Downloads") # can use stringr::str_extract(getwd(), "^.{3}") before "/Users/" if you need to paste working directory letter.
    
    if (!file.exists(paste(data_path, "gshhg-bin-2.3.7", sep = "/"))) {
      
      op <- getOption("timeout")
      options(timeout = max(600, getOption("timeout")))
      
      tmp <- tempfile(fileext = ".zip")
      
      download.file(
        "https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.7.zip",
        tmp
      )
      
      unzip(tmp, 
            exdir = paste(data_path, 
                          "gshhg-bin-2.3.7", 
                          sep = "/"))
      unlink(tmp)
      
      options(timeout = op)
      
    }
    
    shape_file <- paste(data_path, "gshhg-bin-2.3.7/GSHHS_shp/f/GSHHS_f_L1.shp", sep = "/")
    shape_file <- st_read(shape_file)
  }
  

  # Format data and shapefiles--------------------------------------------------------------------------
  # Format Shapefiles 
  sf_use_s2(FALSE)
  ALLisl_poly=shape_file %>% st_make_valid() # could take a couple minutes; st_make_valid() is necessary to catch any errors in points
  
  if("OCC_SITEID" %in% colnames(data)== FALSE) {
    data$OCC_SITEID <- NA
  }
  
  # Filter spc field data (without fixed site SPC)
  ptsFISH = data %>%
    dplyr::select(REGION:REA_SITEID,LATITUDE:LONGITUDE, METHOD, HABITAT_CODE:HABITAT_TYPE, OCC_SITEID) %>%
    filter(METHOD %in% "nSPC") %>% #nSPC indicates SPC field data (excluding fixed sites)
    distinct(.)
  
  # Filter spc fixed site, aka occ fixed site, data only
  ptsOCC = fixed_metadata %>%
            filter(LOCATIONCODE %in% ptsFISH$LOCATIONCODE) %>%
            select(-LOCALDATE, -SITE_DEPTH_M) %>%
            mutate(REA_SITEID = "NA",
                   METHOD = "NA")
            
  
  
  # Reformat spc and field and fixed site dataframes into simple features point format
  ptsOCC=st_as_sf(ptsOCC,coords=c("LONGITUDE","LATITUDE"),crs=st_crs(ALLisl_poly))
  ptsFISH=st_as_sf(ptsFISH,coords=c("LONGITUDE","LATITUDE"),crs=st_crs(ALLisl_poly))
  

  
  
  # Create distance matrix from each fixed target OCC point to every fish SPC survey -------------------
  
  uLC=unique(ptsOCC$LOCATIONCODE) # create vector of unique islands where data was collected from
  Fish_OCC_pts=NULL #create empty vector/df
  
  #Loop Through each island's survey sites and create shortest distances "as a fish swims"
  for(i_isl in 1:length(uLC)){
    #i_isl=1
    thisisl=uLC[i_isl] # for each island...
    ptsFISH_ISL=ptsFISH %>% filter(LOCATIONCODE==uLC[i_isl]) # ...filter data by island
    ptsOCC_ISL=ptsOCC %>% filter(LOCATIONCODE==uLC[i_isl]) # ...filter data by island
  
    #clip shape file to rasterize
    pts_ext=terra::ext(bind_rows(ptsFISH_ISL,ptsOCC_ISL))*1.5 # 1.5 is to zoom out of the spatial extent...can adjust to zoom out more or less
    options(warn = 1)
    poly_ISL=st_crop(x = ALLisl_poly,y = st_bbox(pts_ext))
    options(warn = 0) # more conservative
    
    # if need to change projection of coordinates...AMSM is UTM2
    # ptsOCC_ISL_utm <- st_transform(ptsOCC_ISL, crs = 32702) # 32702 UTMZone2 code OR sprintf("+proj=utm +zone=2 +south +datum=WGS84 units=km")
    # ptsFISH_ISL_utm <- st_transform(ptsFISH_ISL, crs = 32702)
    # poly_ISL_utm <- st_transform(poly_ISL, crs = 32702)
    
    par(mfrow = c(1,1)) #par(mfrow = c(2, length(ptsOCC_ISL$OCC_SITEID)/2))
    
    D=calc_fish_distance_matrix(pointsOCC = ptsOCC_ISL,
                                pointsFISH = ptsFISH_ISL,
                                island_poly = poly_ISL,
                                resolution_m = 50)
  
    OUT <- as.data.frame(as.table(D))
    
    Fish_OCC_pts <- rbind(Fish_OCC_pts,OUT) 
    
    out <- Fish_OCC_pts %>% 
                      dplyr::rename(OCCSITE = Var1,
                                    REASITE = Var2,
                                    DISTANCE_m = Freq)
  }
  
  
  
  # Create Associated fish SPC site labels based on distance from each fixed SPC/OCC site ----------------
  
  # subset distance matrix/dataframe for fish SPC sites within 6000m of each fixed SPC/OCC site
  subset_by_dist <- out %>%
                      mutate(value = case_when(DISTANCE_m <= subset_distance_m ~ "-1",
                                               DISTANCE_m > subset_distance_m ~ "0",))
  
  # number of associated sites with each fixed SPC/OCC site
  #subset_by_dist %>% group_by(OCCSITE, value) %>% count() %>% spread(value, n)
  
  
  
  
  # create final dataframe with associated sites-----------------------------------------------------------
  fixed_dat <- fixed_metadata %>%
                filter(LOCATIONCODE %in% ptsFISH$LOCATIONCODE) %>%
                select(-LOCALDATE, -SITE_DEPTH_M) %>%
                mutate(REA_SITEID = "NA",
                       METHOD = "NA")
  
  finaldf <- data %>%
              dplyr::select(REGION:LONGITUDE, METHOD, HABITAT_CODE:HABITAT_TYPE) %>%
              distinct()
  
  finaldf <- bind_rows(fixed_dat, finaldf)
    
  # filter all fixed SPC/OCC sites within distance chosen (i.e. 1500m) of each fish SPC site
  for (i in 1:nrow(finaldf)) {
    sub.i <- left_join(finaldf[i,], subset_by_dist %>% filter(value != 0), by = c("REA_SITEID" = "REASITE"), relationship = "many-to-many")
    
    finaldf$ASSOC_OCCSITEID[i] =  ifelse(nrow(sub.i) == 0, NA, paste0(sub.i$OCCSITE, collapse=", ")) # paste multiple OCC_SITEID buffers for each point into one cell value separated by ","
  
    finaldf$value[i] = ifelse(finaldf$ASSOC_OCCSITEID[i] %in% "NA", 0, -1) 
    
    }
  
  
  finaldf <- separate_rows(finaldf, ASSOC_OCCSITEID, sep = ", ") # duplicate rows with multiple ASSOC_OCCSITEIDs 
  finaldf$ASSOC_OCCSITEID <- gsub(" ", "", finaldf$ASSOC_OCCSITEID)
  
  # Assign "value" column to 1
  finaldf <- finaldf %>%
              mutate(value = case_when(OCC_SITEID != "" ~ 1,
                                       TRUE ~ value)) %>% #copy OCC_SITEID into ASSOC_OCCSITE for fixed sites only
              mutate(ASSOC_OCCSITEID = case_when(value == "1" ~ OCC_SITEID,
                                                 TRUE ~ ASSOC_OCCSITEID)) %>% #copy OCC_SITEID into ASSOC_OCCSITE for fixed sites only
              mutate(value = case_when(value == -1 & !HABITAT_CODE %in% c("AGR", "APR", "APS", "PPR", "ROB", "SAG", "WAL") ~ 0,
                                       TRUE ~ value)) #change all associated sites (-1) that have pavement or rubble structure to 0 b/c parrotfish do nor forage in these habitat types

  #finaldf %>% select(REA_SITEID, HABITAT_CODE, value) %>% group_by(HABITAT_CODE, value) %>% count() %>% spread(value, n)
  
  # subset associated OCCSITEID to include only fish REA SPC sites that match the habitat types of their respective fixed SPC/OCC site
  fixedhab <- finaldf %>% dplyr::select(OCC_SITEID, HABITAT_CODE) %>% distinct() %>% filter(OCC_SITEID != "")

  output <- left_join(finaldf, fixedhab, by = c("ASSOC_OCCSITEID" = "OCC_SITEID"), relationship = "many-to-many") %>%
              mutate(value = case_when(value == 0 & HABITAT_CODE.x == HABITAT_CODE.y ~ -1, #assign all associated fish SPC sites that do not match the habitat type of their respective fixed SPC/OCC site to 0
                                       TRUE ~ value)) %>%
              dplyr::select(-HABITAT_CODE.y) %>%
              dplyr::rename(HABITAT_CODE = HABITAT_CODE.x)

  # number of associated sites with each fixed SPC/OCC site
  strs_samplesize <- finaldf %>% 
                      dplyr::select(ASSOC_OCCSITEID, HABITAT_CODE, value) %>%
                      dplyr::group_by(ASSOC_OCCSITEID, value) %>% 
                      dplyr::count() %>% 
                      tidyr::spread(value, n) %>% 
                      dplyr::rename(Associated = `-1`) %>% 
                      dplyr::rename(Not_Associated = `0`) %>% 
                      dplyr::rename(Fixed = `1`)

  
  return(list(output = finaldf, surveysamplesize = strs_samplesize))

  }
  
  
  
  
  
  
  
  
#   # Convert to sf (simple features) object and assign coordinates. Define spatial reference as 4326 (WGS 84 geographic coordinate system, units in degrees).
#   sites.lookup.sf = st_as_sf(sites.lookup, coords = c("LONGITUDE", "LATITUDE")) %>% st_set_crs(4326)
#   
#   # Transform to 3857 (projected coordinate system, units in meters)
#   sites.lookup.proj = sites.lookup.sf %>% st_transform(3857)
#   
#   # Add a 12.5km buffer (25km diameter) around each site
#   proj.buffer = st_buffer(sites.lookup.proj, 12500)
#   buffer = proj.buffer %>% st_transform(4326)
#   
#   
#   
#   # Convert to sf, set the coordinates and projection
#   assoc.sites.sf = st_as_sf(assoc.sites, coords = c("LONGITUDE", "LATITUDE")) %>% st_set_crs(4326) %>% st_transform(3857)
#   
#   
#   for (i in 1:nrow(assoc.sites)) {
#     join.i = st_join(assoc.sites.sf[i,], proj.buffer, join = st_within)
#     
#     assoc.sites$ASSOC_OCCSITEID[i] =  ifelse(nrow(join.i) == 0, NA, paste(join.i$OCC_SITEID.y, collapse=", ")) # paste multiple OCC_SITEID buffers for each point into one cell value separated by ","
#     assoc.sites$value[i] = ifelse(assoc.sites$ASSOC_OCCSITEID[i] %in% "NA", 0, -1) 
#     
#   }
#   
#   assoc.sites <- separate_rows(assoc.sites,ASSOC_OCCSITEID,sep = ",") # duplicate rows with multiple ASSOC_OCCSITEIDs 
#   
#   # merge back the fixed site metadata and assign "ASSOCIATED" column to 1
#   final <- sites.lookup %>%
#     mutate(ASSOC_OCCSITEID = NA) %>%
#     mutate(value = 1) %>% # assign value of 1 to fixed site SPC data
#     rbind(assoc.sites) %>%
#     mutate(ASSOC_OCCSITEID = case_when(value == "1" ~ OCC_SITEID,
#                                        TRUE ~ ASSOC_OCCSITEID)) #copy OCC_SITEID into ASSOC_OCCSITE for fixed sites only
#   final$ASSOC_OCCSITEID <- gsub('\\s+', '', final$ASSOC_OCCSITEID)
#   
#   
#   return(final)
# }
