\name{formula_vector_generator}
\alias{formula_vector_generator}
\title{
Molecular Formula Vector Generator
}
\description{
This function convert a molecular formulas into a numerical vector
}
\usage{
formula_vector_generator(molecular_formula, Elements, LElements = length(Elements),
allowedRedundantElements = FALSE)
}
\arguments{
  \item{molecular_formula}{
molecular formula
}
  \item{Elements}{
a string vector of elements. This value must be driven from the `element_sorter` function.
}
  \item{LElements}{
number of elements. To speed up loop calculations, consider calculating the number of elements outside of the loop.
}
  \item{allowedRedundantElements}{
`TRUE` should be used to deconvolute molecular formulas with redundant elements (e.g. CO2CH3O), and `FALSE` should be used to skip such complex molecular formulas.(default value)
}
}
\value{
a numerical vector for the molecular formula. This function returns a vector of -Inf values when the molecular formula has elements not listed in the `Elements` string vector.
}
\examples{
molecular_formula <- "[13]C2C12H2Br5Cl3O"
Elements_molecular_formula <- c("[13]C", "C", "H", "O", "Br", "Cl")
EL <- element_sorter(ElementList = Elements_molecular_formula, alphabeticalOrder = TRUE)
Elements <- EL[["Elements"]]
LElements <- length(Elements)
##
mol_vec <- formula_vector_generator(molecular_formula, Elements, LElements,
allowedRedundantElements = TRUE)
##
regenerated_molecular_formula <- hill_molecular_formula_printer(Elements, mol_vec)
}
