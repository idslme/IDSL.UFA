\name{detect_formula_sets }
\alias{detect_formula_sets }
\title{
Organic Class Detection by Repeated Unit Patterns
}
\description{
This function sorts a vector of molecular formulas to detect organic compound class with reapeted/non-repeated substructure units. This function only works for molecular formulas with following elements: c("As", "Br", "Cl", "Na", "Se", "Si", "B", "C", "F", "H", "I", "K", "N", "O", "P", "S")
}
\usage{
detect_formula_sets(molecular_formulas, ratio_delta_HBrClFI_C,
mixed.HBrClFI.allowed, min_molecular_formula_class, max_number_formula_class,
number_processing_threads = 1)
}
\arguments{
  \item{molecular_formulas}{
a vector of molecular formulas
}
  \item{ratio_delta_HBrClFI_C}{
c(2, 1/2, 0).
2 to detect structures with linear carbon chains such as PFAS, lipids, chlorinated paraffins, etc.
1/2 to detect structures with cyclic chains such as PAHs.
0 to detect molecular formulas with a fixed structures but changing H/Br/Cl/F/I atoms similar to PCBs, PBDEs, etc.
}
  \item{mixed.HBrClFI.allowed}{
mixed.HBrClFI.allowed = c(TRUE, FALSE). Select `FALSE` to detect halogenated-saturated compounds such as PFOS or select `TRUE` to detect mixed halogenated compounds with hydrogen.
}
  \item{min_molecular_formula_class}{
minimum number of molecular formulas in each class. This number should be greater than or equal to 2.
}
  \item{max_number_formula_class}{
maximum number of molecular formulas in each class
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded computations.
}
}
\value{
A matrix of clustered classes of organic molecular formulas.
}
\examples{
molecular_formulas <- c("C3F7O3S", "C4F9O3S", "C5F11O3S", "C6F9O3S", "C8F17O3S",
"C9F19O3S", "C10F21O3S", "C7ClF14O4", "C10ClF20O4", "C11ClF22O4", "C11Cl2F21O4",
"C12ClF24O4")
##
ratio_delta_HBrClFI_C <- 2 # to class polymeric classes
mixed.HBrClFI.allowed <- FALSE # To detect only halogen saturated classes
min_molecular_formula_class <- 2
max_number_formula_class <- 20
##
classes <- detect_formula_sets(molecular_formulas, ratio_delta_HBrClFI_C,
mixed.HBrClFI.allowed, min_molecular_formula_class, max_number_formula_class,
number_processing_threads = 1)
}
