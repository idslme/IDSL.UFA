\name{monoisotopic_mass_calculator}
\alias{monoisotopic_mass_calculator}
\title{
Monoisotopic Mass Calculator
}
\description{
This function calculates monoisotopic mass of a molecular formula
}
\usage{
monoisotopic_mass_calculator(MoleFormVec, Elements_mass_abundance)
}
\arguments{
  \item{MoleFormVec}{
A numerical vector molecular formula
}
  \item{Elements_mass_abundance}{
A list of isotopic mass and abundance of elements obtained from the element_sorter function
}
}
\value{
The monoisotopic mass
}
\examples{
Elements <- c("C", "H", "O")
MoleFormVec <- c(2, 6, 1) # C2H6O
EL_mass_abundance <- element_sorter(ElementList = Elements, ElementOrder = "alphabetical")
Elements_mass_abundance <- EL_mass_abundance[[2]]
MImass <- monoisotopic_mass_calculator(MoleFormVec, Elements_mass_abundance)
}
