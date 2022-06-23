\name{isotopic_profile_calculator}
\alias{isotopic_profile_calculator}
\title{
Isotopic Profile Calculator
}
\description{
This function was designed to calculate isotopic profile distributions for small molecules with masses <= 1200 Da. Details of the equations used in this function are available in the reference[1]. In this function, neighboring isotopologues are merged using the satellite clustering merging (SCM) method described in the reference[2].
}
\usage{
isotopic_profile_calculator(MoleFormVec, Elements_mass_abundance, peak_spacing,
intensity_cutoff, UFA_IP_memeory_variables = c(1e30, 1e-12))
}
\arguments{
  \item{MoleFormVec}{
A numerical vector of the molecular formula
}
  \item{Elements_mass_abundance}{
A list of isotopic mass and abundance of elements obtained from the `element_sorter` function
}
  \item{peak_spacing}{
A maximum space between two isotopologues in Da
}
  \item{intensity_cutoff}{
A minimum intesnity threshold for isotopic profiles in percentage
}
  \item{UFA_IP_memeory_variables}{
A vector of two variables. Default values are c(1e30, 1e-12). Memory may be an issue when the entire isotopologues are calculated; therefore, memeory_variables[1] is used to adjust memory usage. memeory_variables[2] indicates the minimum relative abundance (RA calcaulated by eq(1) in the reference [1]) of an isotopologue to include in the isotopic profile calculations.
}
}
\value{
A matrix of isotopic profile. The first and second column represents the mass
and instensity profiles, respectively.
}
\references{
[1] Fakouri Baygi, S., Crimmins, B.S., Hopke, P.K. Holsen, T.M. (2016). Comprehensive emerging chemical discovery: novel polyfluorinated compounds in Lake Michigan trout. \emph{Environmental Science and Technology}, 50(17), 9460-9468, \doi{10.1021/acs.est.6b01349}.

[2] Fakouri Baygi, S., Fernando, S., Hopke, P.K., Holsen, T.M. and Crimmins, B.S. (2019). Automated Isotopic Profile Deconvolution for High Resolution Mass Spectrometric Data (APGC-QToF) from Biological Matrices. \emph{Analytical chemistry}, 91(24), 15509-15517, \doi{10.1021/acs.analchem.9b03335}.
}
\seealso{
\url{https://ipc.idsl.me/}
}
\examples{
EL <- element_sorter()
Elements <- EL[[1]]
Elements_mass_abundance <- EL[[2]]
peak_spacing <- 0.005 # mDa
intensity_cutoff <- 1 # (in percentage)
MoleFormVec <- formula_vector_generator("C8H10N4O2", Elements)
IP <- isotopic_profile_calculator(MoleFormVec, Elements_mass_abundance,
                                  peak_spacing, intensity_cutoff)
}
