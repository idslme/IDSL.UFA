molecular_formula_library_search <- function(MolecularFormulaAnnotationTable, MFlibrary, IonPathways, number_processing_threads = 1) {
  ##
  molecularFormula <- MolecularFormulaAnnotationTable$IonFormula
  LmolecularFormula <- length(molecularFormula)
  ##
  Elements <- element_sorter()[["Elements"]]
  LElements <- length(Elements)
  ##
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    MolVecMat <- do.call(rbind, lapply(molecularFormula, function(i) {
      formula_vector_generator(i, Elements, LElements, allowedRedundantElements = FALSE)
    }))
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    ############################################################################
    ##
    if (osType == "Linux") {
      MolVecMat <- do.call(rbind, mclapply(molecularFormula, function(i) {
        formula_vector_generator(i, Elements, LElements, allowedRedundantElements = FALSE)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
      ##########################################################################
      ##
    } else if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      MolVecMat <- foreach(i = molecularFormula, .combine = 'rbind', .verbose = FALSE) %dopar% {
        formula_vector_generator(i, Elements, LElements, allowedRedundantElements = FALSE)
      }
      ##
      stopCluster(clust)
      ##
    }
  }
  ##
  ##############################################################################
  ##
  IonPW_DC <- ionization_pathway_deconvoluter(IonPathways, Elements)
  ##
  MF_cbind <- do.call(cbind, lapply(1:length(IonPW_DC), function(p) {
    ##
    if (!is.null(IonPW_DC[[p]])) {
      Ion_coeff <- IonPW_DC[[p]][[1]]
      Ion_adduct <- IonPW_DC[[p]][[2]]
      MolVecMatdeIonized <- do.call(rbind, lapply(1:LmolecularFormula, function(i) {
        MolVecdeIonized <- MolVecMat[i, ]/Ion_coeff - Ion_adduct
        xNeg <- which(MolVecdeIonized < 0)
        if (length(xNeg) > 0) {
          MolVecdeIonized <- rep(0, LElements)
        }
        MolVecdeIonized
      }))
      ##
      molecularFormulaMatrixElementSorterList <- molecular_formula_elements_filter(MolVecMatdeIonized, Elements)
      MolVecMatdeIonized <- molecularFormulaMatrixElementSorterList[["molecularFormulaMatrix"]]
      Elements <- molecularFormulaMatrixElementSorterList[["elementSorterList"]][["Elements"]]
      ##
      MF_intact <- hill_molecular_formula_printer(Elements, MolVecMatdeIonized, number_processing_threads)
    } else {
      MF_intact <- rep(NA, LmolecularFormula)
    }
    freq_p <- MFlibrary[MF_intact]
    xNA <- which(is.na(freq_p))
    freq_p[xNA] <- 0
    names(freq_p) <- NULL
    lib_MF_freq <- cbind(MF_intact, freq_p)
    lib_MF_freq <- data.frame(lib_MF_freq)
    colnames(lib_MF_freq) <- c(paste0("MolecluarFormula ", IonPathways[p]), paste0("Frequency_library_", p))
    ##
    return(lib_MF_freq)
  }))
  ##
  MolecularFormulaAnnotationTable <- cbind(MolecularFormulaAnnotationTable, MF_cbind)
  ##
  return(MolecularFormulaAnnotationTable)
}
