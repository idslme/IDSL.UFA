\name{UFA_score_coefficients_corrector}
\alias{UFA_score_coefficients_corrector}
\title{
Score Coefficients Corrector for MolecularFormulaAnnotationTable
}
\description{
This function updates ranking orders of the individual MolecularFormulaAnnotationTable when score coefficients changed.
}
\usage{
UFA_score_coefficients_corrector(input_annotated_molf_address,
output_annotated_molf_address, scoreCoefficients, number_processing_threads = 1)
}
\arguments{
  \item{input_annotated_molf_address}{
Address to load the individual MolecularFormulaAnnotationTables.
}
  \item{output_annotated_molf_address}{
Address to save the individual MolecularFormulaAnnotationTables.
}
  \item{scoreCoefficients}{
A vector of five numbers representing coefficients of the identification score function.
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded computations.
}
}
