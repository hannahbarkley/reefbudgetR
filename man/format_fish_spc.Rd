% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_fish.R
\name{format_fish_spc}
\alias{format_fish_spc}
\title{Format Stationary Point Count data}
\usage{
format_fish_spc(data, rates_dbase_ = c("IPRB", "Kindinger"))
}
\arguments{
\item{data}{all stationary point count data.}

\item{rates_dbase_}{Erosion rates database to use. Choose either Indo-Pacific
#'ReefBudget ("rates_dbase_ = "IPRB") or U.S. Pacific Islands rates developed
#'by Tye Kindinger, NOAA PIFSC ("rates_dbase_ = "Kindinger").}
}

\description{
Format Stationary Point Count data
}
\examples{
fish_data <- read.csv("CB_FishSPC_alldata.csv", na = "", check.names = FALSE)

fish_spc <- format_fish_spc(data = fish_data, rates_dbase_ = "Kindinger")
}
\author{
Rebecca Weible
}