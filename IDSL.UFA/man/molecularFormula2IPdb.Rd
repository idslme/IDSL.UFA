\name{molecularFormula2IPdb}
\alias{molecularFormula2IPdb}
\title{
Molecular Formula to IPDB
}
\description{
A function to calculate IPDBs from a vector of molecular formulas
}
\usage{
molecularFormula2IPdb(molecularFormulaDatabase, retentionTime = NULL, peak_spacing = 0,
intensity_cutoff_str = 1, IonPathways = "[M]+", number_processing_threads = 1,
UFA_IP_memeory_variables = c(1e30, 1e-12, 100), allowedMustRunCalculation = FALSE,
allowedVerbose = TRUE)
}
\arguments{
  \item{molecularFormulaDatabase}{
A vector string of molecular formulas OR a list of elements and molecular formula matrix
}
  \item{retentionTime}{
retention time
}
  \item{peak_spacing}{
A maximum space between isotopologues in Da to merge neighboring isotopologues.
}
  \item{intensity_cutoff_str}{
A minimum intensity threshold for isotopic profiles in percentage. This parameter may be a string piece of R commands using c, b, br, cl, k, s, se, and si variables corresponding to the same elements.
}
  \item{IonPathways}{
A vector of ionization pathways. Pathways should be like [Coeff*M+ADD1-DED1+...] where "Coeff" should be an integer between 1-9 and ADD1 and DED1 may be ionization pathways. ex: 'IonPathways <- c("[M]+", "[M+H]+", "[2M-Cl]-", "[3M+CO2-H2O+Na-KO2+HCl-NH4]-")'
}
  \item{number_processing_threads}{
number of processing cores for multi-threaded computations.
}
  \item{UFA_IP_memeory_variables}{
A vector of three variables. Default values are c(1e30, 1e-12, 100) to manage memory usage. UFA_IP_memeory_variables[1] is used to control the overall size of isotopic combinations. UFA_IP_memeory_variables[2] indicates the minimum relative abundance (RA calculated by eq(1) in the reference [1]) of an isotopologue to include in the isotopic profile calculations. UFA_IP_memeory_variables[3] is the maximum elapsed time to calculate the isotopic profile on the `setTimeLimit` function of the base R.
}
  \item{allowedMustRunCalculation}{
c(TRUE, FALSE). A `TRUE` allowedMustRunCalculation applies a brute-force method to calculate complex isotopic profiles. When `TRUE`, this option may significantly reduce the speed for multithreaded processing.
}
  \item{allowedVerbose}{
c(TRUE, FALSE). A `TRUE` allowedVerbose provides messages about the flow of the function.
}
}
\value{
An IPDB list of isotopic profiles
}
\references{
[1] Fakouri Baygi, S., Crimmins, B.S., Hopke, P.K. Holsen, T.M. (2016). Comprehensive emerging chemical discovery: novel polyfluorinated compounds in Lake Michigan trout. \emph{Environmental Science and Technology}, 50(17), 9460-9468, \doi{10.1021/acs.est.6b01349}.
}
\seealso{
\url{https://ipc.idsl.me/}
}
\examples{
\donttest{
library(IDSL.UFA, attach.required = TRUE)
molecular_formula <- c("C13F8N8O2", "C20H22", "C8HF16ClSO3", "C12Cl10", "C123H193N35O37")
peak_spacing <- 0.005 # in Da for QToF instruments
# Use this piece of code for intensity cutoff to preserve significant isotopologues
intensity_cutoff_str <- "if (s>0 & si>0) {min(c(c, 10, si*3, s*4))}
else if (s>0 & si==0) {min(c(c, 10, s*4))}
else if (s==0 & si>0) {min(c(c, 10, si*3))}
else if (s==0 & si==0) {min(c(c, 10))}"
UFA_IP_memeory_variables <- c(1e30, 1e-12, 100)
IonPathways <- c("[M+H]+", "[M+Na]+", "[M-H2O+H]+")
number_processing_threads <- 2
listIsoProDataBase <- molecularFormula2IPdb(molecular_formula, retentionTime = NULL,
peak_spacing, intensity_cutoff_str, IonPathways, number_processing_threads,
UFA_IP_memeory_variables, allowedMustRunCalculation = FALSE, allowedVerbose = TRUE)
save(listIsoProDataBase, file = "listIsoProDataBase.Rdata")
}
}
