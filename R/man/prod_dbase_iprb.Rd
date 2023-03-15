\name{prod_dbase_iprb}
\alias{prod_dbase_iprb}
\docType{data}
\title{Indo-Pacific ReefBudget carbonate production database
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("prod_dbase_iprb")}
\format{
  \describe{
    \item{\code{REGIONCODE}}{Region, Main Hawaiian Islands "MHI" or Marianas Archipelago "MARIAN".}
    \item{\code{SUBSTRATE_CODE}}{Substrate code for taxa.}
    \item{\code{CORAL_GROUP_NAME}}{Coral genera group.}
    \item{\code{CORAL_GROUP}}{Coral genera group code.}
    \item{\code{SUBSTRATE_CLASS}}{Substrate classification.}
    \item{\code{SUBSTRATE_NAME}}{Substrate classification name.}
    \item{\code{MORPHOLOGY}}{Morphology (corals only).}
    \item{\code{MORPHOLOGYCODE}}{Morphology code (corals only).}
    \item{\code{TAXA_LEVEL}}{Taxa levels for corals, "GENUS" or "SPECIES".}
    \item{\code{EXTENSION_CM_YR}}{Mean extension rate in cm/yr.}
    \item{\code{EXTENSION_CM_YR_CI}}{Confidence interval for extension rate in cm/yr}
    \item{\code{DENSITY_G_CM3}}{Mean skeletal density in g/cm3.}
    \item{\code{DENSITY_G_CM3_CI}}{Confidence interval for skeletal density in g/cm3}
    \item{\code{CONVERSION_FACTOR}}{Conversion factor for branching morphologies.}
  }
}
\details{
 Indo-Pacific ReefBudget parrotfish erosion rates are derived from the
"Coral taxa average rates" tab in "IP Calcification and bioerosion rates database v.1.3.xlsx",
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
data(prod_dbase_iprb)
}
\keyword{datasets}
