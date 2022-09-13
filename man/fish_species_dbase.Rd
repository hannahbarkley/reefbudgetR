\name{fish_species_dbase}
\alias{fish_species_dbase}
\docType{data}
\title{Fish species database
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("fish_species_dbase")}
\format{
  A data frame with 740 observations on the following 8 variables.
  \describe{
  \item{\code{SPECIES}}{Species code}
  \item{\code{TAXON_NAME}}{Taxon name}
  \item{\code{FAMILY}}{Family}
  \item{\code{COMMONFAMILYALL}}{Common name}
    \item{\code{TROPHIC_MONREP}}{a factor with levels \code{PISCIVORE}, \code{PLANKTIVORE}, \code{PRIMARY}, \code{SECONDARY}, \code{UNKNOWN}}
    \item{\code{LW_A}}{a numeric vector}
    \item{\code{LW_B}}{a numeric vector}
    \item{\code{LENGTH_CONVERSION_FACTOR}}{a numeric vector}
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
\examples{
data(fish_species_dbase)
## maybe str(fish_species_dbase) ; plot(fish_species_dbase) ...
}
\keyword{datasets}
