UFAx_molecular_formula_library_search <- function(molecular_formula_ions, IonPathways, Elements, MF_library, number_processing_threads = 1) {
  L_Elements <- length(Elements)
  ##
  IonPW_DC <- ionization_pathway_deconvoluter(IonPathways, Elements)
  L_PW <- length(IonPathways)
  ##
  L_MolF <- length(molecular_formula_ions)
  ##
  osType <- Sys.info()[['sysname']]
  if (osType == "Windows") {
    clust <- makeCluster(number_processing_threads)
    registerDoParallel(clust)
    ##
    print("Initiated deconvoluting molecular formulas!!!")
    MolVecMat <- foreach(i = 1:L_MolF, .combine = 'rbind', .verbose = FALSE) %dopar% {
      formula_vector_generator(molecular_formula_ions[i], Elements, L_Elements)
    }
    print("Completed deconvoluting molecular formulas!!!")
    ##
    print("Initiated calculating molecular formulas from the ionization pathways!!!")
    MF_cbind <- do.call(cbind, lapply(1:L_PW, function(p) {
      IonPW <- IonPW_DC[[p]]
      Ion_coeff <- IonPW[[1]]
      Ion_adduct <- IonPW[[2]]
      MolVecMat <- MolVecMat/Ion_coeff - Ion_adduct
      if (!is.infinite(Ion_adduct)) {
        MolVecMat_deIonized <- foreach(i = 1:L_MolF, .combine = 'rbind', .verbose = FALSE) %dopar% {
          MolVec_deIonized <- MolVecMat[i, ]/Ion_coeff - Ion_adduct
          x_neg <- which(MolVec_deIonized < 0)
          if (length(x_neg) > 0) {
            MolVec_deIonized <- rep(-Inf, L_Elements)
          }
          MolVec_deIonized
        }
        MF_intact <- hill_molecular_formula_printer(Elements, MolVecMat_deIonized, number_processing_threads)
      } else {
        MF_intact <- rep(NA, L_MolF)
      }
      freq_p <- MF_library[MF_intact]
      x_NA <- which(is.na(freq_p))
      freq_p[x_NA] <- 0
      names(freq_p) <- c()
      lib_MF_freq <- cbind(MF_intact, freq_p)
      lib_MF_freq <- data.frame(lib_MF_freq)
      colnames(lib_MF_freq) <- c(paste0("MolFormula ", IonPathways[p]), "Frequency")
      lib_MF_freq
    }))
    print("Completed calculating molecular formulas from the ionization pathways!!!")
    ##
    stopCluster(clust)
    ##
  } else if (osType == "Linux") {
    ##
    print("Initiated deconvoluting molecular formulas!!!")
    MolVecMat <- do.call(rbind, mclapply(1:L_MolF, function(i) {
      formula_vector_generator(molecular_formula_ions[i], Elements, L_Elements)
    }, mc.cores = number_processing_threads))
    print("Completed deconvoluting molecular formulas!!!")
    ##
    print("Initiated calculating molecular formulas from the ionization pathways!!!")
    MF_cbind <- do.call(cbind, lapply(1:L_PW, function(p) {
      IonPW <- IonPW_DC[[p]]
      Ion_coeff <- IonPW[[1]]
      Ion_adduct <- IonPW[[2]]
      if (!is.infinite(Ion_adduct[1])) {
        MolVecMat_deIonized <- do.call(rbind, mclapply(1:L_MolF, function(i) {
          MolVec_deIonized <- MolVecMat[i, ]/Ion_coeff - Ion_adduct
          x_neg <- which(MolVec_deIonized < 0)
          if (length(x_neg) > 0) {
            MolVec_deIonized <- rep(-Inf, L_Elements)
          }
          MolVec_deIonized
        }, mc.cores = number_processing_threads))
        MF_intact <- hill_molecular_formula_printer(Elements, MolVecMat_deIonized, number_processing_threads)
      } else {
        MF_intact <- rep(NA, L_MolF)
      }
      freq_p <- MF_library[MF_intact]
      x_NA <- which(is.na(freq_p))
      freq_p[x_NA] <- 0
      names(freq_p) <- c()
      lib_MF_freq <- cbind(MF_intact, freq_p)
      lib_MF_freq <- data.frame(lib_MF_freq)
      colnames(lib_MF_freq) <- c(paste0("MolFormula ", IonPathways[p]), "Frequency")
      lib_MF_freq
    }))
    print("Completed calculating molecular formulas from the ionization pathways!!!")
    ##
    closeAllConnections()
  }
  ##
  MF_cbind <- data.frame(MF_cbind)
  colnames(MF_cbind) <- NULL
  rownames(MF_cbind) <- NULL
  ##
  return(MF_cbind)
}