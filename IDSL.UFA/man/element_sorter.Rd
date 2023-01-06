\name{element_sorter}
\alias{element_sorter}
\title{
Element Sorter
}
\description{
This module sorts 84 non-labeled and 14 labeled elements in the periodic table for molecular formula deconvolution and isotopic profile calculation.
}
\usage{
element_sorter(ElementList = "all", alphabeticalOrder = TRUE)
}
\arguments{
  \item{ElementList}{
A string vector of elements needed for isotopic profile calculation. The default value for this parameter is a vector string of entire elements.
}
  \item{alphabeticalOrder}{
`TRUE` should be used to sort the elements for elemental deconvolution (default value), `FALSE` should be used to keep the input order.
}
}
\value{
  \item{Elements}{
A string vector of elements (alphabetically sorted or unsorted)
}
  \item{massAbundanceList}{
A list of isotopic mass and abundance of elements.
}
  \item{Valence}{
A vector of electron valences.
}
}
\examples{
EL_mass_abundance_val <- element_sorter()
}
