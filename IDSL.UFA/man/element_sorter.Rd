\name{element_sorter}
\alias{element_sorter}
\title{
Element Sorter
}
\description{
This function sorts 84 elements in the periodic table for molecular formula deconvolution and isotopic profile calculation.
}
\usage{
element_sorter(ElementList = "all", ElementOrder = "alphabetical")
}
\arguments{
  \item{ElementList}{
A string vector of elements needed for isotopic profile calculation. The default value for this parameter is a vector string of entire elements.
}
  \item{ElementOrder}{
ElementOrder = c("alphabetical", "same") where "alphabetical" should be used to sort the elements for elemental deconvolution (default value), "same" should be used to keep the input order.
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
