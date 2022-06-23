\name{molecular_formula_library_generator}
\alias{molecular_formula_library_generator}
\title{
Molecular Formula Database Producer
}
\description{
This function produce an efficient database for molecular formula matching against a database.
}
\usage{
molecular_formula_library_generator(entire_molecular_formulas)
}
\arguments{
  \item{entire_molecular_formulas}{
A string vector of molecular formulas (redundancy is allowed)
}
}
\value{
A vector of frequency of molecular formulas in the databse.
}
\examples{
entire_molecular_formulas <- c("C2H6O", "C2H6O", "C2H6O", "C2H6O", "CH4O", "CH4O",
"CH4O", "NH4", "C6H12O6")
db <- molecular_formula_library_generator(entire_molecular_formulas)
freq <- db[c("C6H12O6", "CH4O")]
}
