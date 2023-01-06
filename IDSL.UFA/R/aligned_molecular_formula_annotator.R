aligned_molecular_formula_annotator <- function(PARAM) {
  ##
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Initiated generating data for the aligned molecular formula annotated table!")
  ##
  number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2])
  peak_alignment_folder <- PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0014'), 2]
  output_path_annotated_mf_tables <- paste0(output_path, "/annotated_mf_tables")
  mf_table_list <- dir(path = output_path_annotated_mf_tables, pattern = ".Rdata")
  ##
  peakXcol <- IDSL.IPA::loadRdata(paste0(peak_alignment_folder, "/peakXcol.Rdata"))
  ColPL <- colnames(peakXcol)
  L_peaks <- dim(peakXcol)[1]
  Lsamples3 <- dim(peakXcol)[2]
  Lsamples <- Lsamples3 - 3
  ColPL <- ColPL[4:Lsamples3]
  ##
  ##############################################################################
  ##############################################################################
  ##
  seqXcolSample <- do.call(rbind, lapply(1:length(ColPL), function(i) {
    patternSampleName <- paste0("MolecularFormulaAnnotationTable_", ColPL[i], ".Rdata")
    ##
    xPatternCheck <- grep(patternSampleName, mf_table_list)
    if (length(xPatternCheck) > 0) {
      c(i, xPatternCheck)
    } else {
      c(i, 0)
    }
  }))
  ##
  MissedPL <- which(seqXcolSample[, 2] == 0)
  if (length(MissedPL) > 0) {
    IPA_logRecorder("WARNING!!! MolecularFormulaAnnotationTables are not avialable for the following HRMS files:")
    for (i in MissedPL) {
      IPA_logRecorder(ColPL[i])
    }
    ##
    seqXcolSample <- matrix(seqXcolSample[-MissedPL, ], ncol = 2)
  }
  ##
  seqXcolSample <- matrix(seqXcolSample[order(seqXcolSample[, 2], decreasing = FALSE), ], ncol = 2)
  orderSeqSample <- seqXcolSample[, 1]
  seqSample <- seqXcolSample[, 2]
  peakXcol <- peakXcol[, c(seq(1, 3, 1), (orderSeqSample + 3))]
  ##
  ##############################################################################
  ##############################################################################
  ##
  maxRankSample <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0028'), 2])
  Ncandidate <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0029'), 2])
  ##
  adjustFreqRankCheck <- if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0030'), 2]) == "yes") {TRUE} else {FALSE}
  ##
  ##############################################################################
  ##
  call_MF_Zcol <- function(i) {
    peak_table_id <- peakXcol[, (i + 3)]
    MolecularFormulaAnnotationTable <- IDSL.IPA::loadRdata(paste0(output_path_annotated_mf_tables, "/", mf_table_list[i]))
    matched_peak_ids <- as.numeric(MolecularFormulaAnnotationTable$'IDSL.IPA_PeakID')
    x_peak_ids <- which(peak_table_id %in% unique(matched_peak_ids))
    ##
    if (length(x_peak_ids) > 0) {
      matchedMolecularFormulas <- MolecularFormulaAnnotationTable$'IonFormula'
      ##
      do.call(rbind, lapply(x_peak_ids, function(j) {
        x_j <- which(matched_peak_ids == peak_table_id[j])
        max_k <- min(c(maxRankSample, length(x_j)))
        ##
        cbind(rep(j, max_k), seq(1, max_k, 1), matchedMolecularFormulas[x_j[1:max_k]])
      }))
    }
  }
  ##
  ##############################################################################
  ##
  repNA00Ncandidate <- rep(c(NA, 0, 0), Ncandidate)
  ##
  call_calculating_median_ranks <- function(i) {
    molecularFormulaFreqRank <- repNA00Ncandidate
    ##
    if (xZcol[i, 1] != 0) {
      rankMolecularFormula <- MF_Zcol[xZcol[i, 1]:xZcol[i, 2], 2:3]
      ##
      uMolecularFormula <- unique(rankMolecularFormula$MolecularFormula)
      ##
      molecularFormulaRank <- do.call(rbind, lapply(uMolecularFormula, function(j) {
        jRow <- subset(rankMolecularFormula, MolecularFormula == j)
        jFreq <- dim(jRow)[1]
        jMed <- median(jRow$Rank) # A median is calculated for the rank of candidate compounds across samples
        ##
        c(j, jFreq, jMed)
      }))
      ##
      if (adjustFreqRankCheck) {    # To adjust ranking and frequencies
        oderAdjustFreqRank <- order(sqrt(as.numeric(molecularFormulaRank[, 2]))/as.numeric(molecularFormulaRank[, 3]), decreasing = TRUE)
        molecularFormulaRank <- matrix(molecularFormulaRank[oderAdjustFreqRank, ], ncol = 3)
      } else {
        molecularFormulaRank <- matrix(molecularFormulaRank[order(as.numeric(molecularFormulaRank[, 3]), decreasing = FALSE), ], ncol = 3)
        molecularFormulaRank <- matrix(molecularFormulaRank[order(as.numeric(molecularFormulaRank[, 2]), decreasing = TRUE), ], ncol = 3)
      }
      ##
      minNcandidate <- min(Ncandidate, dim(molecularFormulaRank)[1])
      for (k in 1:minNcandidate) {
        molecularFormulaFreqRank[3*k - 2] <- molecularFormulaRank[k, 1]
        molecularFormulaFreqRank[3*k - 1] <- molecularFormulaRank[k, 2]
        molecularFormulaFreqRank[3*k] <- molecularFormulaRank[k, 3]
      }
    }
    return(molecularFormulaFreqRank)
  }
  ##
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    IPA_logRecorder("Initiated matching peak IDs!")
    progressBARboundaries <- txtProgressBar(min = 0, max = Lsamples, initial = 0, style = 3)
    #
    MF_Zcol <- do.call(rbind, lapply(seqSample, function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      call_MF_Zcol(i)
    }))
    close(progressBARboundaries)
    #
    MF_Zcol <- data.frame(MF_Zcol)
    colnames(MF_Zcol) <- c("XID", "Rank", "MolecularFormula")
    MF_Zcol$XID <- as.numeric(MF_Zcol$XID)
    MF_Zcol$Rank <- as.numeric(MF_Zcol$Rank)
    MF_Zcol <- MF_Zcol[order(MF_Zcol[, 1], decreasing = FALSE), ]
    rownames(MF_Zcol) <- NULL
    xDiff <- which(diff(MF_Zcol[, 1]) > 0)
    #
    xZcol <- matrix(rep(0, 2*L_peaks), ncol = 2)
    #
    u_peakid <- unique(MF_Zcol[, 1])
    xZcol[u_peakid, 1] <- c(1, (xDiff + 1))
    xZcol[u_peakid, 2] <- c(xDiff, dim(MF_Zcol)[1])
    #
    IPA_logRecorder("Completed matching peak IDs!")
    ##
    IPA_logRecorder("Initiated calculating median ranks!")
    progressBARboundaries <- txtProgressBar(min = 0, max = L_peaks, initial = 0, style = 3)
    #
    aligned_molecular_formula <- do.call(rbind, lapply(1:L_peaks, function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      call_calculating_median_ranks(i)
    }))
    MF_Zcol <- NULL
    close(progressBARboundaries)
    IPA_logRecorder("Completed calculating median ranks!")
    ##
    title_mat <- do.call(c, lapply(1:Ncandidate, function(i) {
      c(paste0("IonFormula_", i), paste0("Frequency_", i), paste0("MedianRank_", i))
    }))
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      IPA_logRecorder("Initiated matching peak IDs!")
      MF_Zcol <- do.call(rbind, mclapply(seqSample, function(i) {
        call_MF_Zcol(i)
      }, mc.cores = number_processing_threads))
      #
      MF_Zcol <- data.frame(MF_Zcol)
      colnames(MF_Zcol) <- c("XID", "Rank", "MolecularFormula")
      MF_Zcol$XID <- as.numeric(MF_Zcol$XID)
      MF_Zcol$Rank <- as.numeric(MF_Zcol$Rank)
      MF_Zcol <- MF_Zcol[order(MF_Zcol[, 1], decreasing = FALSE), ]
      rownames(MF_Zcol) <- NULL
      xDiff <- which(diff(MF_Zcol[, 1]) > 0)
      #
      xZcol <- matrix(rep(0, 2*L_peaks), ncol = 2)
      #
      u_peakid <- unique(MF_Zcol[, 1])
      xZcol[u_peakid, 1] <- c(1, (xDiff + 1))
      xZcol[u_peakid, 2] <- c(xDiff, dim(MF_Zcol)[1])
      #
      IPA_logRecorder("Completed matching peak IDs!")
      ##
      IPA_logRecorder("Initiated calculating median ranks!")
      aligned_molecular_formula <- do.call(rbind, mclapply(1:L_peaks, function(i) {
        call_calculating_median_ranks(i)
      }, mc.cores = number_processing_threads))
      MF_Zcol <- NULL
      IPA_logRecorder("Completed calculating median ranks!")
      ##
      title_mat <- do.call(c, mclapply(1:Ncandidate, function(i) {
        c(paste0("IonFormula_", i), paste0("Frequency_", i), paste0("MedianRank_", i))
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      IPA_logRecorder("Initiated matching peak IDs!")
      MF_Zcol <- foreach(i = seqSample, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_MF_Zcol(i)
      }
      #
      MF_Zcol <- data.frame(MF_Zcol)
      colnames(MF_Zcol) <- c("XID", "Rank", "MolecularFormula")
      MF_Zcol$XID <- as.numeric(MF_Zcol$XID)
      MF_Zcol$Rank <- as.numeric(MF_Zcol$Rank)
      MF_Zcol <- MF_Zcol[order(MF_Zcol[, 1], decreasing = FALSE), ]
      rownames(MF_Zcol) <- NULL
      xDiff <- which(diff(MF_Zcol[, 1]) > 0)
      #
      xZcol <- matrix(rep(0, 2*L_peaks), ncol = 2)
      #
      u_peakid <- unique(MF_Zcol[, 1])
      xZcol[u_peakid, 1] <- c(1, (xDiff + 1))
      xZcol[u_peakid, 2] <- c(xDiff, dim(MF_Zcol)[1])
      #
      IPA_logRecorder("Completed matching peak IDs!")
      ##
      IPA_logRecorder("Initiated calculating median ranks!")
      aligned_molecular_formula <- foreach(i = 1:L_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_calculating_median_ranks(i)
      }
      MF_Zcol <- NULL
      IPA_logRecorder("Completed calculating median ranks!")
      ##
      title_mat <-  foreach(i = 1:Ncandidate, .combine = 'c', .verbose = FALSE) %dopar% {
        c(paste0("IonFormula_", i), paste0("Frequency_", i), paste0("MedianRank_", i))
      }
      ##
      stopCluster(clust)
    }
  }
  ##
  medianPeakHeight <- IDSL.IPA::loadRdata(paste0(peak_alignment_folder, "/peak_height.Rdata"))[, 4:5]
  medianPeakArea <- IDSL.IPA::loadRdata(paste0(peak_alignment_folder, "/peak_area.Rdata"))[, 4]
  medianR13C <- IDSL.IPA::loadRdata(paste0(peak_alignment_folder, "/peak_R13C.Rdata"))[, 4]
  ##
  aligned_molecular_formula <- cbind(peakXcol[, 1:3], medianPeakHeight[, 2], medianPeakHeight[, 1], medianPeakArea, medianR13C, aligned_molecular_formula)
  rownames(aligned_molecular_formula) <- NULL
  colnames(aligned_molecular_formula) <- c("mz", "RT", "frequencyPeakXcol", "Flag", "medianPeakHeight", "medianPeakArea", "medianR13C", title_mat)
  ##
  IPA_logRecorder("Completed processing of the peak property table!")
  output_path_aligned_table <- paste0(output_path, "/aligned_molecular_formula_table")
  if (!dir.exists(output_path_aligned_table)) {
    dir.create(output_path_aligned_table, recursive = TRUE)
  }
  IPA_logRecorder("Initiated saving the aligned molecular formula annotated table!")
  save(aligned_molecular_formula, file = paste0(output_path_aligned_table, "/aligned_molecular_formula.Rdata"))
  write.csv(aligned_molecular_formula, file = paste0(output_path_aligned_table, "/aligned_molecular_formula.csv"), row.names = TRUE)
  ##
  ##############################################################################
  ##
  IPA_logRecorder("Stored annotated aligned table as 'aligned_molecular_formula' in the `.Rdata` and `.csv` formats in the `aligned_molecular_formula_table` folder!")
  IPA_logRecorder("Completed generating data for the aligned molecular formula annotated table!")
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ##############################################################################
  ##
  return()
}
