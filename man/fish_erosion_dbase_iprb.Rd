\name{fish_erosion_dbase_iprb}
\alias{fish_erosion_dbase_iprb}
\docType{data}
\title{Fish erosion rates database from Indo-Pacific ReefBudget Methodology
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("fish_erosion_dbase_iprb")}
\format{
  \describe{
    \item{\code{PHASE}}{Life phase, initial ("I") or terminal ("T").}
    \item{\code{SIZE_CLASS}}{Fish size class, in 10cm bins.}
    \item{\code{TAXON_NAME}}{Fish name.}
    \item{\code{EROSION_RATE}}{Fish erosion rates.}
  }
}
\details{ Indo-Pacific ReefBudget parrotfish erosion rates are derived from the
"Equations" tab in "IP_Parrotfish_erosion_rates_database_v1.3.xlsx",
downloaded from the Indo-Pacific ReefBudget website
(https://geography.exeter.ac.uk/reefbudget/indopacific/).
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(fish_erosion_dbase_iprb)
}
\keyword{datasets}
