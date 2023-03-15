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
#'@export format_iprb_erosion_rates
#'

format_iprb_erosion_rates <- function(iprb_rates, output = c("all", "finalerosion")) {

  rates <- iprb_rates %>%

    # Remove blank rows
    dplyr::filter(!X.2 == "") %>%

    # Rename column headings to size bins
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
      I.x20 = "X.11",
      I.x30 = "X.12",
      I.x40 = "X.13",
      I.x50 = "X.14",
      T.x20 = "X.15",
      T.x30 = "X.16",
      T.x40 = "X.17",
      T.x50 = "X.18",
      T.x60 = "X.19"
    ) %>%

    # Split data frame into list of data frames by type of equations
    split(., cumsum(seq_row(.) %in% c(26, 45, 65, 84))) %>%

    lapply(., function(x) x %>%

        # Remove unnecessary rows
        dplyr::filter(!TAXONNAME == "SPECIES")) %>%

    # Rename each data frame in the list by type of equations
    setNames(
      .,
      c(
        "SizeClassErosionRates",
        "PropBiteScars",
        "BiteRates",
        "VolumeRemoved",
        "MassRemoved"
      )
    ) %>%

    # Remove columns that are not necessary
    lapply(., function(x) select(x, -X, -13, -c(X.20:X.25))) %>%

    # Replace blank values in equation list with zeros
    lapply(., function(x) x %>%
        mutate_each(funs(replace(., . == "", 0)), -GrazType, -TAXONNAME, -X.10)) %>%

    # Set up dataframe to create separate phase column
    lapply(., function(x) gather(x, sizebin, count, -GrazType, -TAXONNAME)) %>%

    # Pull phase letters out of sizebin column (J, I, or T)
    lapply(., function(x)
      x %>% separate(
        sizebin,
        c("PHASE", "sizebin"),
        extra = "merge",
        fill = "left"
      )) %>%

    lapply(., function(x) spread(x, sizebin, count, fill = 0)) %>%

    lapply(., function(x) x %>%
             filter(PHASE != "X")) %>%

    # Remove excess columns (these will be calculated later)
    lapply(., function(x) x %>%
             select(-c(x20:x60)))

  output <- rates[1][[1]] %>%
    gather(., "sizeclass", ".value", -c(GrazType:PHASE)) %>%
    select(PHASE, sizeclass, TAXONNAME, .value) %>%
    filter(!(sizeclass %in% "10"))

  # Tell function what to print
  if (output == "all") {
    return(rates)
  }

  if (output == "finalerosion") {
    return(output)
  }



}
