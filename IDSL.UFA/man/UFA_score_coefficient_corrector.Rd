\name{UFA_score_coefficient_corrector}
\alias{UFA_score_coefficient_corrector}
\title{
Score Coefficient Corrector for MolecularFormulaAnnotationTable
}
\description{
This function updates ranking orders of the individual MolecularFormulaAnnotationTable when score coefficients changed.
}
\usage{
UFA_score_coefficient_corrector(input_annotated_molf_address,
output_annotated_molf_address, IPDB_address, maxNEME, Score_coeff,
number_processing_threads = 1)
}
\arguments{
  \item{input_annotated_molf_address}{
Address to load the individual MolecularFormulaAnnotationTables.
}
  \item{output_annotated_molf_address}{
Address to save the individual MolecularFormulaAnnotationTables.
}
  \item{IPDB_address}{
Address of the IPDB (.Rdata).
}
  \item{maxNEME}{
Maximum value for Normalized Euclidean Mass Error (NEME) in mDa
}
  \item{Score_coeff}{
A vector of five numbers representing coefficients of the identification score function.
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded computations.
}
}
