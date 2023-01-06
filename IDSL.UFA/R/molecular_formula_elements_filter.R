molecular_formula_elements_filter <- function(molecularFormulaMatrix, Elements) {
  ##
  xElementNon0 <- do.call('c', lapply(1:length(Elements), function(i) {
    xNon0 <- which(molecularFormulaMatrix[, i] > 0)
    if (length(xNon0) > 0) {
      i
    }
  }))
  ##
  LxElementNon0 <- length(xElementNon0)
  if (LxElementNon0 > 0) {
    Elements <- Elements[xElementNon0]
    molecularFormulaMatrix <- matrix(molecularFormulaMatrix[, xElementNon0], ncol = LxElementNon0)
  } else {
    Elements <- NULL
    molecularFormulaMatrix <- NULL
  }
  ##
  elementSorterList <- element_sorter(ElementList = Elements, alphabeticalOrder = FALSE)
  ##
  molecularFormulaMatrixElementSorterList <- list(molecularFormulaMatrix, elementSorterList)
  names(molecularFormulaMatrixElementSorterList) <- c("molecularFormulaMatrix", "elementSorterList")
  ##
  return(molecularFormulaMatrixElementSorterList)
}
