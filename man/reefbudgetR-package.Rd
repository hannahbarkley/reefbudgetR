\name{reefbudgetR-package}
\alias{reefbudgetR-package}
\alias{reefbudgetR}
\docType{package}

\title{reefbudgetR: Process ReefBudget carbonate budget data}

\description{ The \pkg{reefbudgetR} package tools for working with ReefBudget carbonate budget data, with functions to process field-based and SfM-derived benthic, urchin, and parrotfish census data and calculate carbonate production and erosion.  }

\details{The functions in the \pkg{reefbudgetR} package contain tools for carbonate budget analysis.
}

\examples{
## Process benthic production data
prod = process_prod(
  data = benthic_cpm,
  transect_id = c("A1", "A2", "A3", "B1", "B2", "B3"),
  transect_length = c(10, 10, 10, 10, 10, 10),
  dbase_type = "NCRMP",
  method_name = "CPM",
  data_type = "In water",
  full_summary = TRUE
)

## Process urchin erosion data
urch = process_urchins(
  data = urchins_cpm,
  method_name = "CPM",
  data_type = "In water",
  transect_id = c("A1", "A2", "A3", "B1", "B2", "B3"),
  transect_length = c(10, 10, 10, 10, 10, 10),
  full_summary = TRUE
)
}

