\name{hill_molecular_formula_printer}
\alias{hill_molecular_formula_printer}
\title{
Print Hill Molecular Formula
}
\description{
This function produces molecular formulas from a list numerical vectors in the Hill notation system
}
\usage{
hill_molecular_formula_printer(Elements, MolVecMat, number_processing_threads = 1)
}
\arguments{
  \item{Elements}{
A vector string of the used elements.
}
  \item{MolVecMat}{
A matrix of numerical vectors of molecular formulas in each row.
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A vector of molecular formulas
}
\examples{
Elements <- c("C", "H", "O", "N", "Br", "Cl")
MoleFormVec1 <- c(2, 6, 1, 0, 0, 0) # C2H6O
MoleFormVec2 <- c(8, 10, 2, 4, 0 ,0) # C8H10N4O2
MoleFormVec3 <- c(12, 2, 1, 0, 5, 3) # C12H2Br5Cl3O
MolVecMat <- rbind(MoleFormVec1, MoleFormVec2, MoleFormVec3)
H_MolF <- hill_molecular_formula_printer(Elements, MolVecMat)
}
