formula_adduct_calculator <- function(molecular_formula, IonPathways) {
  ##
  molecular_adducts <- NA
  ##
  EL <- element_sorter()
  Elements <- EL[[1]]
  Elements_mass_abundance <- EL[[2]]
  L_Elements <- length(Elements)
  ## IonPathway = [Coeff*M+CO2-H2O+Na-KO2+HCl-...] # Coeff should be an integer between 1-9
  L_PW <- length(IonPathways)
  IonPW_DC <- ionization_pathway_deconvoluter(IonPathways, Elements)
  ##
  FormulaVector <- formula_vector_generator(molecular_formula, Elements, L_Elements)
  ##
  MolVecMat <- do.call(rbind, lapply(1:L_PW, function(pathway) {
    IonPW <- IonPW_DC[[pathway]]
    Ion_coeff <- IonPW[[1]]
    Ion_adduct <- IonPW[[2]]
    MoleFormVec <- Ion_coeff*FormulaVector + Ion_adduct
    x_neg <- which(MoleFormVec < 0)
    if (length(x_neg) > 0) {
      MoleFormVec <- c()
    }
    MoleFormVec
  }))
  ##
  if (length(MolVecMat) > 0) {
    molecular_adducts <- hill_molecular_formula_printer(Elements, MolVecMat, number_processing_threads = 1)
  }
  return(molecular_adducts)
}