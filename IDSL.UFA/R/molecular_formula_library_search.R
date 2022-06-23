molecular_formula_library_search <- function(MolecularFormulaAnnotationTable, IPDB, MF_library, IonPathways, number_processing_threads = 1) {
  MolVecMat0 <- IPDB[[2]]
  Elements <- MolVecMat0[[1]]
  MolVecMat_DB <- MolVecMat0[[2]]
  L_Elements <- length(Elements)
  FormulaID <- as.numeric(MolecularFormulaAnnotationTable[, 2])
  L_Table <- length(FormulaID)
  MolVecMat <- MolVecMat_DB[FormulaID, ]
  L_MolVecMat <- dim(MolVecMat)[1]
  ##
  IonPW_DC <- ionization_pathway_deconvoluter(IonPathways, Elements)
  ##
  MF_cbind <- do.call(cbind, lapply(1:length(IonPW_DC), function(p) {
    IonPW <- IonPW_DC[[p]]
    Ion_coeff <- IonPW[[1]]
    Ion_adduct <- IonPW[[2]]
    if (!is.infinite(Ion_adduct[1])) {
      MolVecMat_deIonized <- do.call(rbind, lapply(1:L_MolVecMat, function(i) {
        MolVec_deIonized <- MolVecMat[i, ]/Ion_coeff - Ion_adduct
        x_neg <- which(MolVec_deIonized < 0)
        if (length(x_neg) > 0) {
          MolVec_deIonized <- rep(-Inf, L_Elements)
        }
        MolVec_deIonized
      }))
      MF_intact <- hill_molecular_formula_printer(Elements, MolVecMat_deIonized, number_processing_threads)
    } else {
      MF_intact <- rep(NA, L_MolVecMat)
    }
    freq_p <- MF_library[MF_intact]
    x_NA <- which(is.na(freq_p))
    freq_p[x_NA] <- 0
    names(freq_p) <- c()
    lib_MF_freq <- cbind(MF_intact, freq_p)
    lib_MF_freq <- data.frame(lib_MF_freq)
    colnames(lib_MF_freq) <- c(paste0("MolecluarFormula ", IonPathways[p]), paste0("Frequency_library_", p))
    lib_MF_freq
  }))
  MolecularFormulaAnnotationTable <- cbind(MolecularFormulaAnnotationTable, MF_cbind)
  return(MolecularFormulaAnnotationTable)
}
