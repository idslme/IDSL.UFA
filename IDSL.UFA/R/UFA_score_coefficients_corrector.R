UFA_score_coefficients_corrector <- function(input_annotated_molf_address, output_annotated_molf_address, scoreCoefficients, number_processing_threads = 1) {
  ## To re-rank individual annotated molecular formula tables
  ann_peaklist <- dir(input_annotated_molf_address, pattern = ".Rdata$")
  file_names <- gsub(".Rdata$", "", ann_peaklist)
  ##
  if (!dir.exists(output_annotated_molf_address)) {
    dir.create(output_annotated_molf_address, recursive = TRUE)
  }
  ##
  call_UFA_score_coefficients_corrector <- function(i) {
    ##
    annotated_molf <- IDSL.IPA::loadRdata(paste0(input_annotated_molf_address, "/", ann_peaklist[i]))
    colnames_annotated_molf <- colnames(annotated_molf)
    ##
    size_IP <- as.numeric(annotated_molf$'isotopologueCount')
    NEME <- as.numeric(annotated_molf$'normEucMassError')
    PCS <- as.numeric(annotated_molf$'PCS')
    R13C_PL <- as.numeric(annotated_molf$'r13cPeaklist')
    R13C_IP <- as.numeric(annotated_molf$'r13cTheoretical')
    RCS <- as.numeric(annotated_molf$'ratioChromScan')
    ##
    IdentificationScore <- identificationScoreCalculator(scoreCoefficients, size_IP, PCS, RCS, NEME, R13C_PL, R13C_IP)
    ##
    IPApeaks <- as.numeric(annotated_molf$'IDSL.IPA_PeakID')
    xDiff <- c(0, which(abs(diff(IPApeaks)) > 0), length(size_IP))
    ##
    annotated_molf_updated <- do.call(rbind, lapply(1:(length(xDiff) - 1), function(j) {
      x_p <- (xDiff[j] + 1):xDiff[j + 1]
      order_A <- order(IdentificationScore[x_p], decreasing = TRUE)
      A <- annotated_molf[x_p[order_A], ]
      A[, 15] <- seq(1, length(x_p), 1)
      return(A)
    }))
    annotated_molf_updated <- data.frame(annotated_molf_updated)
    colnames(annotated_molf_updated) <- colnames_annotated_molf
    rownames(annotated_molf_updated) <- NULL
    ##
    MolecularFormulaAnnotationTable <- annotated_molf_updated
    save(MolecularFormulaAnnotationTable, file = paste0(output_annotated_molf_address, "/", file_names[i], ".Rdata"))
    write.csv(MolecularFormulaAnnotationTable, file = paste0(output_annotated_molf_address, "/", file_names[i], ".csv"), row.names = TRUE)
    ##
    return()
  }
  ##
  if (number_processing_threads == 1) {
    ##
    null_var <- do.call(rbind, lapply(1:length(file_names), function(i) {
      call_UFA_score_coefficients_corrector(i)
    }))
    ##
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, setdiff(ls(), c("clust")), envir = environment())
      ##
      null_var <- do.call(rbind, parLapply(clust, 1:length(file_names), function(i) {
        call_UFA_score_coefficients_corrector(i)
      }))
      ##
      stopCluster(clust)
      ##
    } else {
      ##
      null_var <- do.call(rbind, mclapply(1:length(file_names), function(i) {
        call_UFA_score_coefficients_corrector(i)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
    }
  }
}
