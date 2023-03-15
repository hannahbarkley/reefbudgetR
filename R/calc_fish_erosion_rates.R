#'Calculate erosion rate metrics from bite rate, volume, and proportion of scars
#'
#'@author Rebecca Weible
#'
#'@param data Data frame with bite rate, bite volume, and prop scars.
#'@param substrate_density Substrate density in g cm^-3; default to 1.47.
#'@param perc_day_feeding Percent of day that Chlorurus gibbus and
#'large parrotfish spend eating (from Bellwood et al. 1995); default to 83.3.
#'
#'@import dplyr
#'@importFrom rlang .data
#'
#'@export calc_fish_erosion_rates

calc_fish_erosion_rates <- function(data,
                                      substrate_density = 1.47,
                                      perc_day_feeding = 83.3) {

  df %>%
    # Calculate Bites Leaving Scars per Minute
    mutate(BITES_SCARS_MIN = .data$`Bite Rate` * .data$`Proportion of Scars`) %>%
    # Calculate Volume Removed per Day
    mutate(VOLUME_ERODED_DAY = BITES_SCARS_MIN *
             `Bite Volume` * 60 * (12 * (perc_day_feeding / 100))) %>%
    # Calculate Mass removed per day (converted g to kg with /1000)
    mutate(MASS_ERODED_KG_DAY = (VOLUME_ERODED_DAY *
                                   substrate_density) / 1000) %>%
    # Calculate Mass (kg) removed per year
    mutate(MASS_ERODED_KG_YEAR = MASS_ERODED_KG_DAY * 365) %>%
    # Final erosion rates = mass (kg) removed per year
    mutate(EROSION_RATE_KG_YEAR = MASS_ERODED_KG_YEAR)

}
