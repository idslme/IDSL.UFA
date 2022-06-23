\name{UFAx_molecular_formula_library_search}
\alias{UFAx_molecular_formula_library_search}
\title{
Searching Molecular Formula in a Library of Known Compounds
}
\description{
This function searches for detected molecular formula ions in a library of known molecular formulas.
}
\usage{
UFAx_molecular_formula_library_search(molecular_formula_ions, IonPathways, Elements,
MF_library, number_processing_threads = 1)
}
\arguments{
  \item{molecular_formula_ions}{
A string vector of molecular formula ions.
}
  \item{IonPathways}{
A vector of ionization pathways. Pathways should be like [Coeff*M+ADD1-DED1+...] where "Coeff" should be an integer between 1-9 and ADD1 and DED1 may be ionization pathways. ex: 'IonPathways <- c("[M]+", "[M+H]+", "[2M-Cl]-", "[3M+CO2-H2O+Na-KO2+HCl-NH4]-")'
}
  \item{Elements}{
A vector string of the used elements.
}
  \item{MF_library}{
A library of molecular formulas generated using the 'molecular_formula_library_generator' module.
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A matrix of matched compounds in the reference library.
}