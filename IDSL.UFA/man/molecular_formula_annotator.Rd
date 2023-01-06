\name{molecular_formula_annotator}
\alias{molecular_formula_annotator}
\title{
Molecular Formula Annotator
}
\description{
This module annotates candidate molecular formulas in the peaklists from the IDSL.IPA pipeline using isotopic profiles.
}
\usage{
molecular_formula_annotator(IPDB, spectraList, peaklist, selectedIPApeaks,
massAccuracy, maxNEME, minPCS, minNDCS, minRCS, scoreCoefficients, RTtolerance = NA,
correctedRTpeaklist = NULL, exportSpectraParameters = NULL, number_processing_threads = 1)
}
\arguments{
  \item{IPDB}{
An isotopic profile database produced by the IDSL.UFA functions.
}
  \item{spectraList}{
a list of mass spectra in each chromatogram scan.
}
  \item{peaklist}{
Peaklist from the IDSL.IPA pipeline
}
  \item{selectedIPApeaks}{
selected IPA peaklist
}
  \item{massAccuracy}{
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
  \item{scoreCoefficients}{
A vector of five numbers representing coefficients of the identification score
}
  \item{RTtolerance}{
Retention time tolerance (min)
}
  \item{correctedRTpeaklist}{
corrected retention time peaklist
}
  \item{exportSpectraParameters}{
Parameters for export MS/MS match figures
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A dataframe of candidate molecular formulas
}
