\name{molecular_formula_library_search}
\alias{molecular_formula_library_search}
\title{
Molecular Formula Library Search
}
\description{
This function attempts to match candidate molecular formulas against a library of molecular formulas using a set of ionization pathways.
}
\usage{
molecular_formula_library_search(MolecularFormulaAnnotationTable,
IPDB, MF_library, IonPathways, number_processing_threads = 1)
}
\arguments{
  \item{MolecularFormulaAnnotationTable}{
A molecular formula annotation table from the 'molecular_formula_annotator' module.
}
  \item{IPDB}{
A list of isotopic profiles database for the targeted compounds.
}
  \item{MF_library}{
A library of molecular formulas generated using the 'molecular_formula_library_generator' module.
}
  \item{IonPathways}{
A vector of ionization pathways. Pathways should be like [Coeff*M+ADD1-DED1+...] where "Coeff" should be an integer between 1-9 and ADD1 and DED1 may be ionization pathways. ex: 'IonPathways <- c("[M]+", "[M+H]+", "[2M-Cl]-", "[3M+CO2-H2O+Na-KO2+HCl-NH4]-")'
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
