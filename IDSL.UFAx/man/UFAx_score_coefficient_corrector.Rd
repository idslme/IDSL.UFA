\name{UFAx_score_coefficient_corrector}
\alias{UFAx_score_coefficient_corrector}
\title{
Score Coefficient MolecularFormulaAnnotationTable Corrector
}
\description{
This function updates ranking orders of the individual MolecularFormulaAnnotationTable when score coefficients changed. This function creates files with similar names to the original files with "_updated" estensions.
}
\usage{
UFAx_score_coefficient_corrector(annotated_molf_address, maxNEME,
Score_coeff, number_processing_threads = 1)
}
\arguments{
  \item{annotated_molf_address}{
Address of the individual MolecularFormulaAnnotationTables (.Rdata).
}
  \item{maxNEME}{
Maximum value for Normalized Euclidean Mass Error (NEME) in mDa
}
  \item{Score_coeff}{
A vector of five numbers representing coefficients of the identification score
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded computations.
}
}
\value{
The MolecularFormulaAnnotationTable is saved with an "_updated" estension.
}