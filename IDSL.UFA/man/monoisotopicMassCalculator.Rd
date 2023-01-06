\name{monoisotopicMassCalculator}
\alias{monoisotopicMassCalculator}
\title{
Monoisotopic Mass Calculator
}
\description{
This function calculates monoisotopic mass of a molecular formula
}
\usage{
monoisotopicMassCalculator(MoleFormVec, massAbundanceList,
LElements = length(massAbundanceList))
}
\arguments{
  \item{MoleFormVec}{
A numerical vector molecular formula
}
  \item{massAbundanceList}{
A list of isotopic mass and abundance of elements obtained from the `element_sorter` function
}
  \item{LElements}{
length of elements
}
}
\value{
The monoisotopic mass
}
\examples{
Elements <- c("C", "H", "O")
MoleFormVec <- c(2, 6, 1) # C2H6O
EL_mass_abundance <- element_sorter(ElementList = Elements, alphabeticalOrder = FALSE)
massAbundanceList <- EL_mass_abundance[["massAbundanceList"]]
MImass <- monoisotopicMassCalculator(MoleFormVec, massAbundanceList)
}
