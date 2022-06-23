\name{formula_adduct_calculator}
\alias{formula_adduct_calculator}
\title{
Formula Adduct Calculator
}
\description{
a function that takes a formula and an vecotr of ionization pathways and returns the adduct formulas.
}
\usage{
formula_adduct_calculator(molecular_formula, IonPathways)
}
\arguments{
  \item{molecular_formula}{
molecular formula
}
  \item{IonPathways}{
A vector of ionization pathways. Pathways should be like [Coeff*M+ADD1-DED1+...] where "Coeff" should be an integer between 1-9 and ADD1 and DED1 may be ionization pathways. ex: 'IonPathways <- c("[M]+", "[M+H]+", "[2M-Cl]-", "[3M+CO2-H2O+Na-KO2+HCl-NH4]-")'
}
}
\value{
A vector of adduct formulas
}
\examples{
molecular_formula = "C15H10O7"
IonPathways = c("[M+]","[M+H]","[M+H20+H]","[M+Na]")
Formula_adducts <- formula_adduct_calculator(molecular_formula, IonPathways)
}