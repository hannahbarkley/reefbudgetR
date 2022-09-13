#'Convert all individually saved Chris Perry's Excel Entry sheets for each site into one large dataframe
#'@author Rebecca Weible
#'@param cpm_csv .csv file in Chris Perry's Excel sheet format without any changes
#'@param METADATA .csv file with all of the metadata from Chris Perry's Excel sheet saved in a dataframe format
#'@param ADDLOCATION create column for unique, general site name where the data comes from (i.e. "Kaneohe)
#'@param OCC.SITEID create column for unique OCC site ID (i.e. OCC-OAH-005)
#'@import dplyr
#'@export convert_cpm

convert_cpm <- function(cpm_csv, cpm_metadata, add_location, OCC_SITEID) {

  data <- read.csv(file = cpm_csv) %>% #call in CPM fish data entry spreadsheet
    dplyr::filter(!TRANSECT.1=="") %>% #remove blank rows
    dplyr::rename(SPECIES = TRANSECT.1, #rename column headings to size bins
                  `J.0-10cm` = X,
                  `I.11-20cm` = 'X.1',
                  `I.21-30cm` = 'X.2',
                  `I.31-40cm` = 'X.3',
                  `I.41-50cm` = 'X.4',
                  `T.11-20cm` = 'X.5',
                  `T.21-30cm` = 'X.6',
                  `T.31-40cm` = 'X.7',
                  `T.41-50cm` = 'X.8',
                  `T.51-60cm` = 'X.9') %>%
    dplyr::filter(!SPECIES == "Species") %>% #remove unnecessary rows
    split(., cumsum(1:nrow(.)%in% c(20,40,60,80,100,120,140,160,180))) %>% #split dataframe into list of dataframes by transect
    lapply(., function(x) x %>% replace(., .=="", 0)) %>% #replace blank values with zeros
    lapply(., function(x) x %>% filter(!SPECIES %in% c("TRANSECT 2", "TRANSECT 3", "TRANSECT 4", "TRANSECT 5", "TRANSECT 6", "TRANSECT 7", "TRANSECT 8", "TRANSECT 9", "TRANSECT 10"))) %>% #remove more unnecessary rows
    setNames(., c("1","2","3","4","5","6","7","8","9","10")) %>%#rename each dataframe in the list by transect number
    bind_rows(., .id="Transect") %>%
    filter(!SPECIES %in% "Total") %>%
    #complete(SPECIES) #no need b/c CPM Excel sheet is already complete
    gather(., sizebin, count, -Transect, -SPECIES) %>% #set up dataframe to create separate phase column
    separate(sizebin, c("PHASE", "sizebin"), extra = "merge", fill = "left") %>% #pull phase letters out of sizebin column (J, I, or T)
    spread(., sizebin, count, fill = 0) %>% #spread back out like data entry format
    select(Transect:SPECIES, `0-10cm`:`51-60cm`, PHASE) %>% #reorder columns
    mutate(CRUISEID = "MP2108") %>% #here and below, add columns of metadata that are missing
    mutate(LOCATIONCODE = "OAH") %>%
    mutate(OCC_SITEID = OCC_SITEID) %>%
    mutate(CB_METHOD = "CPM_BELT") %>%
    mutate(LOCATION = add_location) %>%
    dplyr::rename(TRANSECT = Transect) %>%
    mutate_at(vars(TRANSECT), as.factor) %>%
    mutate(QC1 = "") %>%
    mutate(NOTES = "")


  finaldata <- read.csv(file = cpm_metadata) %>% #import rest of metadata from site description tab in CPM excel sheets
    mutate_at(vars(TRANSECT), as.factor) %>%
    left_join(., data, by = c("OCC_SITEID", "TRANSECT")) %>% #merge with dataframe created above
    ungroup() %>%
    select(CRUISEID, LOCALDATE, LOCATIONCODE, OCC_SITEID, LOCATION, CB_METHOD, DIVER:HABITAT_TYPE, everything())

  return(finaldata)

}







