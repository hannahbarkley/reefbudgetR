% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_urchins.R
\name{process_urchins}
\alias{process_urchins}
\title{Calculate urchin erosion rates from urchin census data}
\usage{
process_urchins(
  data,
  transect_id = NULL,
  transect_length = NULL,
  method_name = c("IPRB", "Chords"),
  full_summary = TRUE
)
}
\arguments{
\item{data}{Urchin observation data set.}

\item{transect_id}{String of transect names (e.g., ("A1", "A2", "A3", "B1", "B2", "B3")). Defaults to NULL.}

\item{transect_length}{String of transect lengths in meters (e.g., c(10, 10, 10, 10, 10, 10)). Defaults to NULL.}

\item{method_name}{Transect design by which data were collected ("IPRB" or "Chords").}
}
\description{
Calculate urchin erosion rates from urchin census data
}
\examples{
urchin_data <- read.csv("ESD_CarbBudget_Urchins_OAHU_2021.csv",
  na = "", check.names = FALSE)

urch_iprb <- process_urchins(
  data = urchin_data[urchin_data$CB_METHOD == "IPRB", ],
  method_name = "IPRB"
)

urch_chords <- process_urchins(
  data = urchin_data[urchin_data$CB_METHOD == "Chords", ],
  method_name = "Chords"
)
}
\author{
Hannah Barkley
}
