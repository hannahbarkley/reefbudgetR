#'Convert all individual Indo-Pacific ReefBudget (IPRB) Excel
#'Entry sheets for each site into one large data frame
#'
#'@author Rebecca Weible
#'
#'@param iprb_csv .csv file in Indo-Pacific ReefBudget Excel sheet format,
#'downloaded as is.
#'@param iprb_metadata .csv file with all of the metadata from
#'Indo-Pacific ReefBudget Excel sheet, saved in a dataframe format.
#'@param add_location create column for unique, general site name
#'where the data comes from ("Kaneohe Bay", "Tumon Bay", ...).
#'@param occ_siteid OCC site ID ("OCC-OAH-005", "OCC-GUA_015", ...).
#'
#'@import dplyr
#'@importFrom sjmisc seq_row
#'
#'@export convert_iprb

convert_iprb <- function(iprb_csv, iprb_metadata, add_location, occ_siteid) {

  data <-

    # Call in Indo-Pacific ReefBudget fish data entry spreadsheet
    read.csv(file = iprb_csv) %>%

    # Remove blank rows
    dplyr::filter(!TRANSECT.1 == "") %>%

    # Rename column headings to size bins
    dplyr::rename(SPECIES = TRANSECT.1,
                  `J.0-10cm` = X,
                  `I.11-20cm` = "X.1",
                  `I.21-30cm` = "X.2",
                  `I.31-40cm` = "X.3",
                  `I.41-50cm` = "X.4",
                  `T.11-20cm` = "X.5",
                  `T.21-30cm` = "X.6",
                  `T.31-40cm` = "X.7",
                  `T.41-50cm` = "X.8",
                  `T.51-60cm` = "X.9") %>%
    # Remove unnecessary rows
    dplyr::filter(!SPECIES == "Species") %>%

    # Split dataframe into list of dataframes by transect
    split(., cumsum(1:seq_row(.) %in%
                      c(20, 40, 60, 80, 100, 120, 140, 160, 180))) %>%

    # Replace blank values with zeros
    lapply(., function(x) x %>%
             replace(., . == "", 0)) %>%

    # Remove more unnecessary rows
    lapply(., function(x) x %>%
             filter(!SPECIES %in% c("TRANSECT 2",
                                    "TRANSECT 3",
                                    "TRANSECT 4",
                                    "TRANSECT 5",
                                    "TRANSECT 6",
                                    "TRANSECT 7",
                                    "TRANSECT 8",
                                    "TRANSECT 9",
                                    "TRANSECT 10"))) %>%

    # Rename each dataframe in the list by transect number
    setNames(., c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) %>%

    bind_rows(., .id = "Transect") %>%

    filter(!SPECIES %in% "Total") %>%

    # Set up dataframe to create separate phase column
    gather(., sizebin, count, -Transect, -SPECIES) %>%

    # Pull phase letters out of sizebin column (J, I, or T)
    separate(sizebin, c("PHASE", "sizebin"), extra = "merge", fill = "left") %>%

    # Spread back out like data entry format
    spread(., sizebin, count, fill = 0) %>%

    # Reorder columns
    select(Transect:SPECIES, `0-10cm`:`51-60cm`, PHASE) %>%

    #Add columns of metadata that are missing
    mutate(CRUISEID = "MP2108") %>%
    mutate(LOCATIONCODE = "OAH") %>%
    mutate(OCC_SITEID = occ_siteid) %>%
    mutate(CB_METHOD = "IPRB_BELT") %>%
    mutate(LOCATION = add_location) %>%
    dplyr::rename(TRANSECT = Transect) %>%
    mutate_at(vars(TRANSECT), as.factor) %>%
    mutate(QC1 = "") %>%
    mutate(NOTES = "")

# Import rest of metadata from site description tab in IPRB excel sheets
  finaldata <- read.csv(file = iprb_metadata) %>%
    mutate_at(vars(TRANSECT), as.factor) %>%

    # Merge with dataframe created above
    left_join(., data, by = c("OCC_SITEID", "TRANSECT")) %>%
    ungroup() %>%
    select(CRUISEID,
           LOCALDATE,
           LOCATIONCODE,
           OCC_SITEID,
           LOCATION,
           CB_METHOD,
           DIVER:HABITAT_TYPE,
           everything())

  return(finaldata)

}
