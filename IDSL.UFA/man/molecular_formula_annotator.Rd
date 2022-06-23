\name{molecular_formula_annotator}
\alias{molecular_formula_annotator}
\title{
Molecular Formula Annotator
}
\description{
This module annotate candidate molecular formulas in the
peaklists from the IDSL.IPA pipeline using isotopic profiles.
}
\usage{
molecular_formula_annotator(IPDB, spectraList, peaklist,
mass_accuracy, maxNEME, minPCS, minNDCS, minRCS, Score_coeff,
number_processing_threads)
}
\arguments{
  \item{IPDB}{
An isotopic profile database produced by the IDSL.UFA functions.
}
  \item{spectraList}{
spectraList from the 'MS_deconvoluter' function of the IDSL.IPA package
}
  \item{peaklist}{
Peaklist from the IDSL.IPA pipeline
}
  \item{mass_accuracy}{
Mass accuracy in Da
}
  \item{maxNEME}{
Maximum value for Normalized Euclidean Mass Error (NEME) in mDa
}
  \item{minPCS}{
Minimum value for Profile Cosine Similarity (PCS)
}
  \item{minNDCS}{
Minimum value for Number of Detected Chromatogram Scans (NDCS)
}
  \item{minRCS}{
Minimum value for Ratio of Chromatogram Scans (RCS) in percentage
}
  \item{Score_coeff}{
A vector of five numbers representing coefficients of the identification score
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A dataframe of candidate molecular formulas
}
