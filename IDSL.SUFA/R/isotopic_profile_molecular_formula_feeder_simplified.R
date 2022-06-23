isotopic_profile_molecular_formula_feeder_simplified <- function(molecular_formula, IonPathways = "[M]+", peak_spacing = 0, intensity_cutoff = 1, UFA_IP_memeory_variables = c(1e30, 1e-12)) {
  IPP <- c()
  ##
  EL <- element_sorter()
  Elements <- EL[[1]]
  Elements_mass_abundance <- EL[[2]]
  L_Elements <- length(Elements)
  ## IonPathway = [Coeff*M+CO2-H2O+Na-KO2+HCl-...] # Coeff should be an integer between 1-9
  IonPW_DC <- ionization_pathway_deconvoluter(IonPathways, Elements)
  ##
  FormulaVector <- formula_vector_generator(molecular_formula, Elements, L_Elements)
  #
  Ion_coeff <- IonPW_DC[[1]][[1]]
  Ion_adduct <- IonPW_DC[[1]][[2]]
  MoleFormVec <- Ion_coeff*FormulaVector + Ion_adduct
  x_neg <- which(MoleFormVec < 0)
  if (length(x_neg) > 0) {
    MoleFormVec <- c()
  }
  #
  if (!is.null(MoleFormVec)) {
    IPP <- isotopic_profile_calculator(MoleFormVec, Elements_mass_abundance, peak_spacing, intensity_cutoff, UFA_IP_memeory_variables)
    IPP[, 1] <- round(IPP[, 1], 7)
    IPP[, 2] <- round(IPP[, 2], 5)
  } else {
    warning("Molecular formula is not consistent with the ionization pathway!")
  }
  return(IPP)
}
