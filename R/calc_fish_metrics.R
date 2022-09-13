#'Calculate bite rate, volume, and proportion of scars from bioerosion equations
#'@author Rebecca Weible
#'@param bio_equations dataframe with Tye's equations of bite rate,
#'bite volume, and prop of scars
#'@param tl2fl_conversion dataframe with Tyes TL to FL conversion constants
#'@param avg_sst Input average Sea Surface Temperature (Celsius) of
#'location where data was collected (HI, MARIAN, AMSM, PRIAS, etc.)
#'@param avg_time Input average time period (HHMM) of location where data were
#'collected (e.g. if collected between 8am and 1pm, average will be 11am)
#'@param avg_depth Input average depth (Meters) at location where data were
#'collected
#'@import tidyverse
#'@import dplyr
#'@export calc_fish_metrics

calc_fish_metrics <- function(bio_equations,
                              tl2fl_conversion,
                              avg_sst,
                              avg_time,
                              avg_depth) {

  # Format conversion df and convert TL to FL below
  constants <- tl2fl_conversion %>%
    gather(., range, metric, -c(SPECIES:TLtoFL)) %>%
    # Pull phase letters out of sizebin column (J, I, or T)
    separate(range, c("phase", "range"), extra = "merge", fill = "left") %>%
    # Create column where the average size will be used for each bin for the TL
    # to FL conversion in a later step
    mutate(avgSIZE = case_when(
      range == "11-20cm" ~ 15,
      range == "21-30cm" ~ 25,
      range == "31-40cm" ~ 35,
      range == "41-50cm" ~ 45,
      range == "51-60cm" ~ 55
    )) %>%
    # Convert TL to FL
    mutate(SIZE_FL = ifelse(conv == "TLtoFL", TLtoFL_a + (TLtoFL_b * avgSIZE),
                            (avgSIZE - FLtoTL_a) / FLtoTL_b)) %>%
    # Remove unnecessary calculation metrics for TL to FL
    select(SPECIES:FISH_sciname, phase:range, avgSIZE:SIZE_FL) %>%
    # Input SST of Hawaii (~79 deg F or 26 deg C)
    mutate(SST_INPUT = avg_sst) %>%
    # Input time of day the surveys were conducted (on average around 11am)
    mutate(TIME_INPUT = avg_time) %>%
    # Input average depth of surveys (~45 ft or 13.7m)
    mutate(DEPTH_INPUT = avg_depth) %>%
    # Convert all columns that should be numeric
    mutate_at(vars(c("SST_INPUT":"DEPTH_INPUT")), as.numeric)

  # Bring in bioerosion equations and calculate metrics for PROP SCARS,
  # BITE RATES, BITE VOLUMES

  output <- bio_equations %>%
    merge(constants, ., by = "FISH_sciname") %>%
    # Replace NA values in select columns to 0 so that we can do the
    # calculations for all species in Tye's list
    mutate_at(vars(SIZE_cm_FL, SST_degC, TIME_day_num, DEPTH_m),
              ~ replace_na(., 0)) %>%
    mutate(metric = case_when(RESP_var == "prop scars" ~
                                (exp((Intercept + SIZE_cm_FL * SIZE_FL) /
                                       (1 + exp(Intercept +
                                                  SIZE_cm_FL * SIZE_FL))
                                )),
                              TRUE ~ (
                                Intercept + (SIZE_cm_FL * SIZE_FL) +
                                  (SST_degC * SST_INPUT) +
                                  (TIME_day_num * TIME_INPUT) +
                                  (DEPTH_m * DEPTH_INPUT)
                              ))) %>%
    # Remove combination of bite rate rows with NA bite analysis
    filter(RESP_var != "bite rate" |
             biterate_analysis != "NA") %>%
    select(FISH_sciname:range, RESP_var, metric) %>%
    spread(., range, metric, fill = "T") %>%
    filter(!`11-20cm` == "T") %>%
    # Change specific columns to factors
    mutate_at(vars(c("FISH_sciname":"RESP_var")), as.factor) %>%
    # Make numeric values b/c need to do math in a sec
    mutate_if(is.character, as.numeric) %>%
    # Round numeric values to 3 decimal places
    mutate(across(where(is.numeric), round, 5)) %>%
    # Select only species that are distinct in the dataframe
    distinct(.) %>%
    right_join(., fish_grazing_types %>% select(-c(available_number:notes)),
               by = "FISH_sciname") %>%
    distinct() %>%
    # Add column for whos data this belongs to
    mutate(type = "Tye") %>%
    mutate_at(vars(c("11-20cm":"51-60cm")), as.numeric)

  return(output)
}
