#' Create fish associated SPC sites in relation to fixed sites
#'
#'@author Rebecca Weible
#'
#'@param data fish survey data.
#'
#'@import Rmisc
#'@import sf
#'@import rgeos
#'@import tidyverse
#'@import dplyr
#'
#'@export create_fish_assoc_sites
#'
#'@examples
#'fish_data <- read.csv("CB_FishBelt_alldata.csv", na = "", check.names = FALSE)
#'
#'sites_associated_dbase <- create_fish_assoc_sites(data = fish_data)



create_fish_assoc_sites <- function(data){
  
  # Filter spc field data (without fixed site SPC) -----------------------------------------------------
  assoc.sites = data %>%
    select(REGION:LONGITUDE, CB_METHOD, HABITAT_CODE:HABITAT_TYPE) %>%
    filter(CB_METHOD %in% "nSPC") %>% #nSPC indicates SPC field data (excluding fixed sites)
    distinct(.)
  

  # Filter spc fixed site data only -----------------------------------------------------
  sites.lookup = data %>%
    select(REGION:LONGITUDE, CB_METHOD, HABITAT_CODE:HABITAT_TYPE) %>%
    filter(CB_METHOD %in% "IPRB") %>% #IPRB indicates fixed OCC sites
    distinct(.)
  
  
  # Create 12.5km buffer around fixed OCC sites ---------------------------------
  
  # Convert to sf (simple features) object and assign coordinates. Define spatial reference as 4326 (WGS 84 geographic coordinate system, units in degrees).
  sites.lookup.sf = st_as_sf(sites.lookup, coords = c("LONGITUDE", "LATITUDE")) %>% st_set_crs(4326)
  
  # Transform to 3857 (projected coordinate system, units in meters)
  sites.lookup.proj = sites.lookup.sf %>% st_transform(3857)
  
  # Add a 12.5km buffer (25km diameter) around each site
  proj.buffer = st_buffer(sites.lookup.proj, 12500)
  buffer = proj.buffer %>% st_transform(4326)
  
  
  
  # Convert to sf, set the coordinates and projection
  assoc.sites.sf = st_as_sf(assoc.sites, coords = c("LONGITUDE", "LATITUDE")) %>% st_set_crs(4326) %>% st_transform(3857)
  
  
  for (i in 1:nrow(assoc.sites)) {
    join.i = st_join(assoc.sites.sf[i,], proj.buffer, join = st_within)
    
    assoc.sites$ASSOC_OCCSITEID[i] =  ifelse(nrow(join.i) == 0, NA, paste(join.i$OCC_SITEID.y, collapse=", ")) # paste multiple OCC_SITEID buffers for each point into one cell value separated by ","
    assoc.sites$value[i] = ifelse(assoc.sites$ASSOC_OCCSITEID[i] %in% "NA", 0, -1) 
    
  }
  
  assoc.sites <- separate_rows(assoc.sites,ASSOC_OCCSITEID,sep = ",") # duplicate rows with multiple ASSOC_OCCSITEIDs 
  
  # merge back the fixed site metadata and assign "ASSOCIATED" column to 1
  final <- sites.lookup %>%
    mutate(ASSOC_OCCSITEID = NA) %>%
    mutate(value = 1) %>% # assign value of 1 to fixed site SPC data
    rbind(assoc.sites) %>%
    mutate(ASSOC_OCCSITEID = case_when(value == "1" ~ OCC_SITEID,
                                       TRUE ~ ASSOC_OCCSITEID)) #copy OCC_SITEID into ASSOC_OCCSITE for fixed sites only
  final$ASSOC_OCCSITEID <- gsub('\\s+', '', final$ASSOC_OCCSITEID)
  
  
  return(final)
}
