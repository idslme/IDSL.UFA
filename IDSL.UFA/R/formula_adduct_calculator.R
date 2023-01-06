formula_adduct_calculator <- function(molecular_formula, IonPathways) {
  ##
  Elements <- element_sorter()[["Elements"]]
  LElements <- length(Elements)
  ##
  IonPW_DC <- ionization_pathway_deconvoluter(IonPathways, Elements)
  ##
  FormulaVector <- formula_vector_generator(molecular_formula, Elements, LElements, allowedRedundantElements = TRUE)
  ##
  MolVecMat <- do.call(rbind, lapply(IonPW_DC, function(pathway) {
    Ion_coeff <- pathway[[1]]
    Ion_adduct <- pathway[[2]]
    MoleFormVec <- Ion_coeff*FormulaVector + Ion_adduct
    xNeg <- which(MoleFormVec < 0)
    if (length(xNeg) > 0) {
      MoleFormVec <- NULL
    }
    MoleFormVec
  }))
  ##
  if (!is.null(MolVecMat)) {
    MolVecMat <- matrix(MolVecMat, ncol = LElements)
    molecular_adducts <- hill_molecular_formula_printer(Elements, MolVecMat)
  } else {
    molecular_adducts <- NA
  }
  return(molecular_adducts)
}
