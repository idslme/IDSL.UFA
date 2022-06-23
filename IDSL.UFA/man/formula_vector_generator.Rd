\name{formula_vector_generator}
\alias{formula_vector_generator}
\title{
Molecular Formula Vector Generator
}
\description{
This function convert a molecular formulas into a numerical vector
}
\usage{
formula_vector_generator(molecular_formula, Elements, L_Elements = length(Elements))
}
\arguments{
  \item{molecular_formula}{
molecular formula
}
  \item{Elements}{
a string vector of elements. This value must be driven from the `element_sorter` function.
}
  \item{L_Elements}{
number of elements. To speed up loop calculations, consider to calculate number of elements outside of the loop.
}
}
\value{
a numerical vector for the molecular formula. This function returns a vector of -Inf values when the molecular formula has elements not listed in the `Elements` string vector.
}
\examples{
molecular_formula <- "C12H2Br5Cl3O"
Elements_molecular_formula <- c("C", "H", "O", "Br", "Cl")
EL <- element_sorter(ElementList = Elements_molecular_formula)
Elements <- EL[[1]]
L_Elements <- length(Elements)
mol_vec <- formula_vector_generator(molecular_formula, Elements, L_Elements)
}
