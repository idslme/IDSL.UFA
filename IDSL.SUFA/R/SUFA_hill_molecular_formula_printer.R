SUFA_hill_molecular_formula_printer <- function(Elements, MolVecMat) {
  HillElements <- sort(Elements, decreasing = FALSE)
  x_c <- which(HillElements == "C")
  x_h <- which(HillElements == "H")
  if ((length(x_c) == 1) & (length(x_h) == 1)) {
    x_c_h <- c(x_c, x_h)
    HillElements <- HillElements[-x_c_h]
    HillElements <- c("C", "H", HillElements)
  } else if ((length(x_c) == 1) & (length(x_h) == 0)) {
    HillElements <- HillElements[-x_c]
    HillElements <- c("C", HillElements)
  }
  ##
  L_Hill_Elements <- length(HillElements)
  x_alpha_hill <- sapply(1:L_Hill_Elements, function(i) {
    which(Elements == HillElements[i])
  })
  ##
  MolVecMat <- matrix(MolVecMat, ncol = L_Hill_Elements)
  ##
  MolFormList <- do.call(rbind, lapply(1:dim(MolVecMat)[1], function(k) {
    MolFor <- c()
    MolVec <- MolVecMat[k, ]
    for (i in 1:L_Hill_Elements) {
      x_vec <- x_alpha_hill[i]
      if (MolVec[x_vec] > 0) {
        if (MolVec[x_vec] == 1) {
          MolFor <- paste0(MolFor, HillElements[i])
        } else {
          MolFor <- paste0(MolFor, HillElements[i], MolVec[x_vec])
        }
      }
    }
    if (is.null(MolFor)) {
      MolFor <- NA
    }
    MolFor
  }))
  return(MolFormList)
}
