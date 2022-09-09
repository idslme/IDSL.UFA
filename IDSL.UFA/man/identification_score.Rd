\name{identification_score}
\alias{identification_score}
\title{
Multiplicative Identification Score for the IDSL.UFA pipeline
}
\description{
This function calculates the score values to rank candidate molecular formulas for a mass spectrometry-chromatography peak.
}
\usage{
identification_score(Score_coefficients, N_isotopologues, PCS, RCS, NEME,
maxNEME, R13C_PL, R13C_IP)
}
\arguments{
  \item{Score_coefficients}{
A vector of seven numbers equal or greater than 0
}
  \item{N_isotopologues}{
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
  \item{maxNEME}{
maximum NEME (mDa)
}
  \item{R13C_PL}{
R13C of the peak from IDSL.IPA peaklists
}
  \item{R13C_IP}{
R13C from theoretical isotopic profiles
}
}
