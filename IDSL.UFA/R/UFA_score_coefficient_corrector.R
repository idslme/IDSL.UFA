UFA_score_coefficient_corrector <- function(input_annotated_molf_address, output_annotated_molf_address, IPDB_address, maxNEME, Score_coeff, number_processing_threads = 1) {
  ## To re-rank individual annotated molecular formula tables
  IPDB <- IDSL.IPA::loadRdata(IPDB_address)
  size_IPDB <- IPDB[[6]]
  IPDB <- 0
  ann_peaklist <- dir(input_annotated_molf_address, pattern = '.Rdata')
  file_names <- gsub(".Rdata", "", ann_peaklist)
  ##
  if (!dir.exists(output_annotated_molf_address)) {
    dir.create(output_annotated_molf_address)
  }
  ##
  null_var_call <- function(i) {
    ##
    annotated_molf <- IDSL.IPA::loadRdata(paste0(input_annotated_molf_address, "/", ann_peaklist[i]))
    colnames_annotated_molf <- colnames(annotated_molf)
    ##
    PCS <- as.numeric(annotated_molf[, 9])
    RCS <- as.numeric(annotated_molf[, 13])
    NEME <- as.numeric(annotated_molf[, 8])
    R13C_PL <- as.numeric(annotated_molf[, 10])
    R13C_IP <- as.numeric(annotated_molf[, 11])
    size_IP <- size_IPDB[as.numeric(annotated_molf[, 2])]
    ##
    IdentificationScore <- identification_score(Score_coeff, size_IP, PCS, RCS, NEME, maxNEME, R13C_PL, R13C_IP)
    ##
    peaks <- unique(annotated_molf[, 1])
    annotated_molf_updated <- do.call(rbind, lapply(1:length(peaks), function(j) {
      x_p <- which(annotated_molf[, 1] == peaks[j])
      order_A <- order(IdentificationScore[x_p], decreasing = TRUE)
      A <- annotated_molf[x_p[order_A], ]
      A[, 14] <- seq(1:length(x_p))
      return(A)
    }))
    annotated_molf_updated <- data.frame(annotated_molf_updated)
    colnames(annotated_molf_updated) <- colnames_annotated_molf
    rownames(annotated_molf_updated) <- c()
    ##
    MolecularFormulaAnnotationTable <- annotated_molf_updated
    save(MolecularFormulaAnnotationTable, file = paste0(output_annotated_molf_address, "/", file_names[i], ".Rdata"))
    write.csv(MolecularFormulaAnnotationTable, file = paste0(output_annotated_molf_address, "/", file_names[i], ".csv"))
    c()
  }
  ##
  osType <- Sys.info()[['sysname']]
  ##
  if (osType == "Windows") {
    clust <- makeCluster(number_processing_threads)
    registerDoParallel(clust)
    ##
    null_var <- foreach(k = 1:length(file_names), .combine = 'rbind', .verbose = FALSE) %dopar% {
      null_var_call(k)
    }
    ##
    stopCluster(clust)
    ##
  } else if (osType == "Linux") {
    ##
    null_var <- do.call(rbind, mclapply(1:length(file_names), function(k) {
      null_var_call(k)
    }, mc.cores = number_processing_threads))
    ##
    closeAllConnections()
  }
}
