% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_fish_erosion_rates.R
\name{calc_fish_erosion_rates}
\alias{calc_fish_erosion_rates}
\title{Calculate erosion rate metrics from bite rate, volume, and proportion of scars}
\usage{
calc_fish_erosion_rates(
  data,
  substrate_density = 1.47,
  perc_day_feeding = 83.3
)
}
\arguments{
\item{data}{Data frame with bite rate, bite volume, and prop scars.}

\item{substrate_density}{Substrate density in g cm^-3; default to 1.47.}

\item{perc_day_feeding}{Percent of day that Chlorurus gibbus and
large parrotfish spend eating (from Bellwood et al. 1995); default to 83.3.}
}
\description{
Calculate erosion rate metrics from bite rate, volume, and proportion of scars
}
\author{
Rebecca Weible
}
