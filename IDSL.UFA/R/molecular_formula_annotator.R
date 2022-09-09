molecular_formula_annotator <- function(IPDB, spectraList, peaklist, mass_accuracy, maxNEME, minPCS, minNDCS, minRCS, Score_coeff, number_processing_threads = 1) {
  MolecularFormulaAnnotationTable <- NULL
  ##
  molecular_formula_annotator_call <- function(k) {
    ##
    mzList.m <- NULL
    ##
    massPLround <- round(peaklist[k, 8], digits = 2)
    IDkIPDB <- do.call(c, lapply(c(-0.02, -0.01, 0, 0.01, 0.02), function(j) {
      IPDB[["AggregatedList"]][[as.character(massPLround + j)]]
    }))
    ##
    if (!is.null(IDkIPDB)) {
      x_mzF <- which(abs(peaklist[k, 8] - IPDB[["MassMAIso"]][IDkIPDB]) <= mass_accuracy)
      if (length(x_mzF) > 0) {
        PL_R13C <- peaklist[k, 11]
        RangeScan <- peaklist[k, 1]:peaklist[k, 2]
        NumberScans <- length(RangeScan)
        IDjIPDB <- IDkIPDB[x_mzF]
        ##
        mzList.m <- do.call(rbind, lapply(IDjIPDB, function(j) {
          ##
          IsotopicProfile <- IPDB[["IsotopicProfile"]][[j]]
          IPsize <- IPDB[["IPsize"]][j]
          MW_exp <- matrix(rep(0, IPsize*NumberScans), ncol = NumberScans)
          INT_exp <- MW_exp
          for (sc in 1:NumberScans) {
            PEAKS <- spectraList[[RangeScan[sc]]]
            for (Iso in 1:IPsize) {
              x_Iso <- which(abs(PEAKS[, 1] - IsotopicProfile[Iso, 1]) <= mass_accuracy)
              if (length(x_Iso) > 0) {
                if (length(x_Iso) > 1) {
                  x_Iso_min <- which.min(abs(PEAKS[x_Iso, 1] - IsotopicProfile[Iso, 1]))
                  x_Iso <- x_Iso[x_Iso_min[1]]
                }
                MW_exp[Iso, sc] <- PEAKS[x_Iso, 1]
                INT_exp[Iso, sc] <- PEAKS[x_Iso, 2]
              }
            }
          }
          sum_INT_exp <- rowSums(INT_exp)
          x_INT_0 <- which(sum_INT_exp == 0)
          if (length(x_INT_0) == 0) {
            PCS <- sum(sum_INT_exp*IsotopicProfile[, 2])/sqrt(sum(sum_INT_exp^2)*sum(IsotopicProfile[, 2]^2))*1000 # in per-mille
            if (PCS >= minPCS) {
              Ave_MW_exp <- rowSums(MW_exp*INT_exp)/sum_INT_exp
              NEME <- sqrt(sum((Ave_MW_exp - IsotopicProfile[, 1])^2)/IPsize)*1000 # in mDa
              if (NEME <= maxNEME) {
                MW_exp[which(MW_exp > 0)] <- 1
                nd <- colSums(MW_exp)
                x_100 <- IPDB[["IndexMAIso"]][j]
                Int_100 <- INT_exp[x_100, ]
                max_Int <- max(Int_100)
                x_80 <- which(Int_100/max_Int > 0.2)
                NDCS <- length(which(nd[x_80] == IPsize))
                if (NDCS >= minNDCS) {
                  L_80 <- x_80[length(x_80)] - x_80[1] + 1
                  RCS <- NDCS/L_80*100
                  if (RCS >= minRCS) {
                    IP_R13C <- IPDB[["R13C"]][j]
                    IdentificationScore <- identification_score(Score_coeff, IPsize, PCS, RCS, NEME, maxNEME, PL_R13C, IP_R13C)
                    c(k, j,
                      IsotopicProfile[x_100, 1],
                      peaklist[k, 8],
                      peaklist[k, 3],
                      max_Int,
                      NEME,
                      PCS,
                      PL_R13C,
                      IP_R13C,
                      NDCS,
                      RCS,
                      IdentificationScore)
                  }
                }
              }
            }
          }
        }))
        if (!is.null(mzList.m)) {
          mzList.m <- matrix(mzList.m, ncol = 13)
          mzList.m <- matrix(mzList.m[order(mzList.m[, 13], decreasing = TRUE), ], ncol = 13)
          mzList.m <- cbind(mzList.m, 1:nrow(mzList.m))
        }
      }
    }
    return(mzList.m)
  }
  ##
  n_peaks <- dim(peaklist)[1]
  if (n_peaks > 0) {
    ##
    if (number_processing_threads == 1) {
      mzList <- do.call(rbind, lapply(1:n_peaks, function(k) {
        molecular_formula_annotator_call(k)
      }))
    } else {
      osType <- Sys.info()[['sysname']]
      if (osType == "Windows") {
        clust <- makeCluster(number_processing_threads)
        registerDoParallel(clust)
        mzList <- foreach(k = 1:n_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
          molecular_formula_annotator_call(k)
        }
        stopCluster(clust)
        ##
      } else if (osType == "Linux") {
        mzList <- do.call(rbind, mclapply(1:n_peaks, function(k) {
          molecular_formula_annotator_call(k)
        }, mc.cores = number_processing_threads))
        closeAllConnections()
      }
    }
    ##
    if (!is.null(mzList)) {
      mzList <- matrix(mzList, ncol = 14)
      mzList <- matrix(mzList[, -13], ncol = 13)
      ##
      mzList[, 3] <- round(mzList[, 3], 5)
      mzList[, 4] <- round(mzList[, 4], 5)
      mzList[, 5] <- round(mzList[, 5], 3)
      mzList[, 6] <- round(mzList[, 6], 0)
      mzList[, 7] <- round(mzList[, 7], 2)
      mzList[, 8] <- round(mzList[, 8], 0)
      mzList[, 9] <- round(mzList[, 9], 2)
      mzList[, 10] <- round(mzList[, 10], 2)
      mzList[, 12] <- round(mzList[, 12], 2)
      ##
      Elements <- IPDB[["MolecularFormulaDB"]][["Elements"]]
      ##
      MolVecMat <- IPDB[["MolecularFormulaDB"]][["MolecularFormulaMatrix"]][mzList[, 2], ]
      MolFlist <- hill_molecular_formula_printer(Elements, MolVecMat, number_processing_threads)
      ##
      MolecularFormulaAnnotationTable <- data.frame(cbind(matrix(mzList[, 1:2], ncol = 2), MolFlist, matrix(mzList[, 3:13], ncol = 11)))
      colnames(MolecularFormulaAnnotationTable) <- c("PeakID", "ID_IonFormula", "IonFormula", "m/z IsotopicProfile", "m/z peaklist", "RT(min)", "PeakHeight", "NEME(mDa)", "PCS", "R13C peakList", "R13C IsotopicProfile", "NDCS", "RCS(%)", "Rank")
      rownames(MolecularFormulaAnnotationTable) <- NULL
    }
  }
  return(MolecularFormulaAnnotationTable)
}
