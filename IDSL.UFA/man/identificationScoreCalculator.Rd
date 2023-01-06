\name{identificationScoreCalculator}
\alias{identificationScoreCalculator}
\title{
Multiplicative Identification Score for the IDSL.UFA pipeline
}
\description{
This function calculates the score values to rank candidate molecular formulas for a mass spectrometry-chromatography peak.
}
\usage{
identificationScoreCalculator(scoreCoefficients, nIisotopologues, PCS, RCS, NEME,
R13C_PL, R13C_IP)
}
\arguments{
  \item{scoreCoefficients}{
A vector of seven numbers equal or greater than 0
}
  \item{nIisotopologues}{
Number of isotopologues in the theoretical isotopic profiles.
}
  \item{PCS}{
PCS (per mille)
}
  \item{RCS}{
RCS (percentage)
}
  \item{NEME}{
NEME (mDa)
}
  \item{R13C_PL}{
R13C of the peak from IDSL.IPA peaklists
}
  \item{R13C_IP}{
R13C from theoretical isotopic profiles
}
}
