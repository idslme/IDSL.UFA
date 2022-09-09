formula_adduct_calculator <- function(molecular_formula, IonPathways) {
  ##
  molecular_adducts <- NA
  ##
  EL <- element_sorter()
  Elements <- EL[["Elements"]]
  L_Elements <- length(Elements)
  ## IonPathway = [Coeff*M+CO2-H2O+Na-KO2+HCl-...] # Coeff should be an integer between 1-9
  IonPW_DC <- ionization_pathway_deconvoluter(IonPathways, Elements)
  ##
  FormulaVector <- formula_vector_generator(molecular_formula, Elements, L_Elements)
  ##
  MolVecMat <- do.call(rbind, lapply(IonPW_DC, function(pathway) {
    Ion_coeff <- pathway[[1]]
    Ion_adduct <- pathway[[2]]
    MoleFormVec <- Ion_coeff*FormulaVector + Ion_adduct
    x_neg <- which(MoleFormVec < 0)
    if (length(x_neg) > 0) {
      MoleFormVec <- NULL
    }
    MoleFormVec
  }))
  ##
  if (length(MolVecMat) > 0) {
    molecular_adducts <- hill_molecular_formula_printer(Elements, MolVecMat, number_processing_threads = 1)
  }
  return(molecular_adducts)
}