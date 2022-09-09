\name{IUPAC_Isotopes}
\alias{IUPAC_Isotopes}
\docType{data}
\title{
IUPAC Isotopes
}
\description{
This data consists of element, mass, abundance and valence of 289 stable isotopes for 84 elements.
}
\usage{data("IUPAC_Isotopes")}
\format{
  A data frame with 289 observations on the following 4 variables.
  \describe{
    \item{\code{element}}{a character vector}
    \item{\code{mass}}{a numerical vector}
    \item{\code{abundance}}{a numerical vector}
    \item{\code{valence}}{a numerical vector}
  }
}
\references{
Kim S, Chen J, Cheng T, Gindulyte A, He J, He S, Li Q, Shoemaker BA, Thiessen PA, Yu B, Zaslavsky L, Zhang J, Bolton EE. PubChem in 2021: new data content and improved web interfaces. Nucleic Acids Res. 2021 Jan 8; 49(D1):D1388–D1395. doi:10.1093/nar/gkaa971.
}
\note{
The PubChem source for isotopes abundance and mass data is IUPAC.
}
\examples{
data(IUPAC_Isotopes)
}
\keyword{datasets}
