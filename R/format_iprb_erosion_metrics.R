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

format_iprb_erosion_metrics <- function(iprb_rates, 
                                      output = c("all", "finalerosion")) {

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
    split(., cumsum(1:nrow(.)%in% c(26,45,65,84))) %>% #split dataframe into list of dataframes by type of equations
    lapply(., function(x) x %>% dplyr::filter(!TAXONNAME == "SPECIES")) %>% #remove unnecessary rows
    setNames(., c("SizeClassErosionRates","PropBiteScars","BiteRates","VolumeRemoved", "MassRemoved")) %>%

    # Remove columns that are not necessary
    lapply(., function(x) select(x, -X, -13)) %>%

    # Replace blank values in equation list with zeros
    lapply(., function(x) x %>% mutate_at(vars(-GrazType, -TAXONNAME, -X.10), ~(replace(., .=="", 0)))) %>%

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

  
  
  # format output dataframes
  output.metrics <- rates %>%
                      map_df(., ~as.data.frame(.x), .id = "metric") %>% # unlist
                      filter(!metric %in% c("SizeClassErosionRates", "MassRemoved")) %>%
                      mutate(metric = case_when(metric == 'PropBiteScars' ~ "Proportion of Scars",
                                                metric == 'BiteRates' ~ "Bite Rate",
                                                metric == 'VolumeRemoved' ~ "Bite Volume")) %>%
                      dplyr::rename(FISH_sciname = "TAXONNAME") %>%
                      mutate(FISH_sciname = strsplit(FISH_sciname, "/")) %>% #split grouped species in cpm metrics 
                      unnest(FISH_sciname) %>% # separate the split of names into separate rows
                      mutate(FISH_sciname=replace(FISH_sciname, FISH_sciname=="ocellatus", "Cetoscarus ocellatus"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="microrhinos", "Chlorurus microrhinos"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="gibbus", "Chlorurus gibbus"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="spilurus", "Chlorurus spilurus"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="bleekeri", "Chlorurus bleekeri"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="capistratoides", "Chlorurus capistratoides"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="japanensis", "Chlorurus japanensis"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="longiceps", "Hipposcarus longiceps"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="falcipinnis", "Scarus falcipinnis"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="altipinnis", "Scarus altipinnis"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="persicus", "Scarus persicus"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="forsteni", "Scarus forsteni"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="flavipectoralis", "Scarus flavipectoralis"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="fuscopurpureus", "Scarus fuscopurpureus"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="russelii", "Scarus russelii"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="rivulatus", "Scarus rivulatus"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="chameleon", "Scarus chameleon"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="festivus", "Scarus festivus"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="oviceps", "Scarus oviceps"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="dimidiatus", "Scarus dimidiatus"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="viridifucatus", "Scarus viridifucatus"),
                             FISH_sciname=replace(FISH_sciname, FISH_sciname=="xanthopleura", "Scarus xanthopleura")) %>% # complete fish scientific names
                      mutate(type = "CPM") %>% # add column for whos data this belongs to            
                      mutate_at(vars(c("11-20cm":"51-60cm")), as.numeric) %>% # all sizeclasses to numbers
                      mutate(across(where(is.numeric), round, 5)) %>% #round numeric values to 5 decimal places
                      select(-`10`) %>% # remove 0-10cm size bin
                      gather(., "sizeclass", "value", -c(metric:PHASE), -type) %>% # gather size bins into single column
                      unite("sizeclass", c(PHASE, sizeclass)) %>% # combine phase and sizeclass into one column
                      filter(!FISH_sciname %in% "Substrate density (g cm-3)") %>% #remove substrate density rows
                      spread(., metric, value, fill = "0") %>% # spread by three metrics for calculations down the road
                      mutate_at(vars(c("Bite Rate":"Proportion of Scars")), as.numeric) %>% # convert characters to numeric for next steps
                      right_join(., GrazTypes %>% select(GRAZ_TYPE, SPECIES, sister_sp, FISH_sciname), by = "FISH_sciname") %>%
                      select(-GrazType) %>%
                      complete(., sizeclass, type, nesting(FISH_sciname, GRAZ_TYPE, SPECIES, sister_sp)) %>%
                      filter(!sizeclass == "NA") %>%
                      mutate_at(vars(`Bite Rate`, `Bite Volume`, `Proportion of Scars`), ~replace_na(., 0)) %>%
                      arrange(FISH_sciname) %>%
                      filter(!type %in% NA) %>%
                      select(SPECIES, FISH_sciname, sizeclass, everything(.), -sister_sp, -type, -GRAZ_TYPE)
                      
    
  output.rates <- rates[1][[1]] %>%
                    gather(., "sizeclass", ".value", -c(GrazType:PHASE)) %>%
                    select(PHASE, sizeclass, TAXONNAME, .value) %>%
                    filter(!(sizeclass %in% "10"))

  # Tell function what to print
  if (output == "all") {
    return(output.metrics)
  }

  if (output == "finalerosion") {
    return(output.rates)
  }



}
