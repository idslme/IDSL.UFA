\name{extendedSENIORrule}
\alias{extendedSENIORrule}
\title{
extended SENIOR rule check
}
\description{
This function checks whether a molecular formula follows the extended SENIOR rule.
}
\usage{
extendedSENIORrule(mol_vec, valence_vec, ionization_correction = 0)
}
\arguments{
  \item{mol_vec}{
A vector of the deconvoluted molecular formula
}
  \item{valence_vec}{
A vector of the valences from the molecular formula. Valences may be acquired from the 'IUPAC_Isotopes' data.
}
  \item{ionization_correction}{
A number to compensate for the ionization losses/gains. For example, `-1` for [M+H/K/Na] ionization pathways and `+1' for [M-H] ionization pathway.
}
}
\value{
\item{rule2}{
TURE for when the molecular formula passes the rule and FALSE for when the molecular formula fails to pass the rule.
}
}
