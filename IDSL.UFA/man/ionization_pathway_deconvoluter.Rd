\name{ionization_pathway_deconvoluter}
\alias{ionization_pathway_deconvoluter}
\title{
Ionization Pathway Deconvoluter
}
\description{
This function deconvolutes ionization pathways into a coefficient and a numerical vector to simplify prediction ionization pathways.
}
\usage{
ionization_pathway_deconvoluter(IonPathways, Elements)
}
\arguments{
  \item{IonPathways}{
A vector of ionization pathways. Pathways should be like [Coeff*M+ADD1-DED1+...] where "Coeff" should be an integer between 1-9 and ADD1 and DED1 may be ionization pathways. ex: 'IonPathways <- c("[M]+", "[M+H]+", "[2M-Cl]-", "[3M+CO2-H2O+Na-KO2+HCl-NH4]-")'
}
  \item{Elements}{
A vector string of the used elements
}
}
\value{
A list of adduct calculation values for each ionization pathway.
}
\examples{
Elements <- element_sorter()[["Elements"]]
IonPathways <- c("[M]+", "[M+H]+", "[2M-Cl]-", "[3M+CO2-H2O+Na-KO2+HCl-NH4]-")
Ion_DC <- ionization_pathway_deconvoluter(IonPathways, Elements)
}
