\name{isotopic_profile_calculator}
\alias{isotopic_profile_calculator}
\title{
Isotopic Profile Calculator
}
\description{
This function was designed to calculate isotopic profile distributions for small molecules with masses <= 1200 Da. Nonetheless, this function may suit more complicated tasks with complex biological compounds. Details of the equations used in this function are available in the reference[1]. In this function, neighboring isotopologues are merged using the satellite clustering merging (SCM) method described in the reference[2].
}
\usage{
isotopic_profile_calculator(MoleFormVec, massAbundanceList, peak_spacing,
intensity_cutoff, UFA_IP_memeory_variables = c(1e30, 1e-12, 100))
}
\arguments{
  \item{MoleFormVec}{
A numerical vector of the molecular formula
}
  \item{massAbundanceList}{
A list of isotopic mass and abundance of elements obtained from the `element_sorter` function
}
  \item{peak_spacing}{
A maximum space between two isotopologues in Da
}
  \item{intensity_cutoff}{
A minimum intensity threshold for isotopic profiles in percentage
}
  \item{UFA_IP_memeory_variables}{
A vector of three variables. Default values are c(1e30, 1e-12, 100) to manage memory usage. UFA_IP_memeory_variables[1] is used to control the overall size of isotopic combinations. UFA_IP_memeory_variables[2] indicates the minimum relative abundance (RA calculated by eq(1) in the reference [1]) of an isotopologue to include in the isotopic profile calculations. UFA_IP_memeory_variables[3] is the maximum elapsed time to calculate the isotopic profile on the `setTimeLimit` function of base R.
}
}
\value{
A matrix of isotopic profile. The first and second column represents the mass and intensity profiles, respectively.
}
\references{
[1] Fakouri Baygi, S., Crimmins, B.S., Hopke, P.K. Holsen, T.M. (2016). Comprehensive emerging chemical discovery: novel polyfluorinated compounds in Lake Michigan trout. \emph{Environmental Science and Technology}, 50(17), 9460-9468, \doi{10.1021/acs.est.6b01349}.

[2] Fakouri Baygi, S., Fernando, S., Hopke, P.K., Holsen, T.M. and Crimmins, B.S. (2019). Automated Isotopic Profile Deconvolution for High Resolution Mass Spectrometric Data (APGC-QToF) from Biological Matrices. \emph{Analytical chemistry}, 91(24), 15509-15517, \doi{10.1021/acs.analchem.9b03335}.
}
\seealso{
\url{https://ipc.idsl.me/}
}
\examples{
EL <- element_sorter(alphabeticalOrder = TRUE)
Elements <- EL[["Elements"]]
massAbundanceList <- EL[["massAbundanceList"]]
peak_spacing <- 0.005 # mDa
intensity_cutoff <- 1 # (in percentage)
MoleFormVec <- formula_vector_generator("C8H10N4O2", Elements)
IP <- isotopic_profile_calculator(MoleFormVec, massAbundanceList, peak_spacing,
intensity_cutoff)
}
