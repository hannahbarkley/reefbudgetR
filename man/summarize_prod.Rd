% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_prod.R
\name{summarize_prod}
\alias{summarize_prod}
\title{Summarize carbonate production rates at transect and site level}
\usage{
summarize_prod(
  data,
  transect_summary,
  dbase_type = c("IPRB", "NCRMP"),
  summarize_by = c("substrate code", "substrate class", "coral group", "overall"),
  level = c("transect", "site"),
  macro_rate = 0.209,
  macro_rate_ci = 0.129,
  micro_rate = 0.262,
  micro_rate_ci = 0.18,
  ...
)
}
\arguments{
\item{data}{Data set to summarize; product of \code{run_calc_prod}.}

\item{transect_summary}{Transect sumamry; product of \code{run_calc_prod}.}

\item{dbase_type}{Production database to use, either Indo-Pacific ReefBudget ("IPRB")
or U.S. Pacific Islands NCRMP-specific database ("NCRMP"). The Indo-Pacific ReefBudget
database is derived from "IP Calcification and bioerosion rates database v.1.3",
downloaded from https://geography.exeter.ac.uk/reefbudget/indopacific/.}

\item{summarize_by}{Grouping factor to summarize by ("substrate code",
"substrate class", "coral group",or "overall")}

\item{level}{Summarize at "transect" or "site" level.}

\item{macro_rate}{Rate of macrobioerosion in kg/cm2/yr. Default is 0.209.}

\item{macro_rate_ci}{Confidence interval for rate of macrobioerosion in kg/cm2/yr, Default is 0.129.}

\item{micro_rate}{Rate of microbioerosion in kg/cm2/yr. Default is 0.262.}

\item{micro_rate_ci}{Confidence interval for rate of microbioerosion in kg/cm2/yr. Default is 0.180.}
}
\description{
Summarize carbonate production rates at transect and site level
}
\examples{
calc_prod_output <- run_calc_prod(
    data = data,
    method_name = "IPRB",
    dbase_type = "NCRMP")

data <- calc_prod_output$data
transect_summary <- calc_prod_output$transect_summary

prod_transect_substratecode <-
  summarize_prod(data,
                 transect_summary,
                 dbase_type,
                 summarize_by = "substrate code",
                 level = "transect")

prod_transect_substrateclass <-
  summarize_prod(data,
                 transect_summary,
                 dbase_type,
                 summarize_by = "substrate class",
                 level = "transect")

prod_transect_coral <-
  summarize_prod(data,
                 transect_summary,
                 dbase_type,
                 summarize_by = "coral group",
                 level = "transect")

prod_transect <-
  summarize_prod(data,
                 transect_summary,
                 dbase_type,
                 summarize_by = "overall",
                 level = "transect")
}
\author{
Hannah Barkley
}
