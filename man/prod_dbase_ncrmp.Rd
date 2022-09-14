\name{NCRMP Carbonate production rates database}
\alias{prod_dbase_ncrmp}
\docType{data}
\title{NCRMP U.S. Pacific Islands carbonate production rates database
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("prod_dbase_ncrmp")}
\format{
  A data frame with 209 observations on the following 17 variables.
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
    \item{\code{EXTENSION_SUBSTITUTE}}{Substitute used for extension.}
    \item{\code{DENSITY_SUBSTITUTE}}{Substitute used for density.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}

\keyword{datasets}
