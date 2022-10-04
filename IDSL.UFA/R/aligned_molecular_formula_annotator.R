aligned_molecular_formula_annotator <- function(PARAM) {
  ##
  UFA_logRecorder(paste0(rep("", 100), collapse = "-"))
  UFA_logRecorder("Initiated generating data for the aligned molecular formula annotated table!")
  ##
  number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == "PARAM0009"), 2])
  input_path_peakXcol <- PARAM[which(PARAM[, 1] == 'PARAM0025'), 2]
  peakXcol <- IDSL.IPA::loadRdata(input_path_peakXcol)
  L_peaks <- dim(peakXcol)[1]
  L_samples2 <- dim(peakXcol)[2]
  L_samples <- L_samples2 - 2
  ##
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0014'), 2]
  output_path_annotated_mf_tables <- paste0(output_path, "/annotated_mf_tables")
  mf_table_list <- dir(path = output_path_annotated_mf_tables, pattern = ".Rdata")
  ##
  ##############################################################################
  ##
  ColPL <- colnames(peakXcol)[3:L_samples2]
  ##
  seqSample <- do.call(c, lapply(1:L_samples, function(i) {
    patternSampleName <- paste0("MolecularFormulaAnnotationTable_", ColPL[i], ".Rdata")
    ##
    patternCheck <- grep(patternSampleName, mf_table_list)
    if (length(patternCheck) > 0) {
      i
    }
  }))
  ##
  MissedPL <- setdiff((1:L_samples), seqSample)
  if (length(MissedPL) > 0) {
    UFA_logRecorder("WARNING!!! The following MolecularFormulaAnnotationTables are not avialable:")
    for (i in MissedPL) {
      UFA_logRecorder(ColPL[i])
    }
  }
  ##
  ##############################################################################
  ##
  input_path_peak_property <- PARAM[which(PARAM[, 1] == 'PARAM0026'), 2]
  peak_property <- IDSL.IPA::loadRdata(input_path_peak_property)
  if (dim(peak_property)[1] != L_peaks | dim(peak_property)[2] != (L_samples2)) {
    stop(UFA_logRecorder("Aligned peak property and peak indexed tables are not in the same size!"))
  }
  peak_property <- peak_property[, c(1, 2, (seqSample + 2))]
  ##
  maxRankSample <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0027'), 2])
  Ncandidate <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0028'), 2])
  ##
  adjustFreqRankCheck <- ifelse((tolower(PARAM[which(PARAM[, 1] == 'PARAM0029'), 2]) == "yes"), TRUE, FALSE)
  ##
  call_MF_Zcol <- function(i) {
    peak_table_id <- peakXcol[, (i + 2)]
    MolecularFormulaAnnotationTable <- IDSL.IPA::loadRdata(paste0(output_path_annotated_mf_tables, "/", mf_table_list[i]))
    matched_peak_ids <- as.numeric(MolecularFormulaAnnotationTable[, 1])
    x_peak_ids <- which(peak_table_id %in% unique(matched_peak_ids) == TRUE)
    ##
    if (length(x_peak_ids) > 0) {
      matchedMolecularFormulas <- MolecularFormulaAnnotationTable$IonFormula
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
  call_calculating_median_ranks <- function(i) {
    molecularFormulaFreqRank <- rep(c(NA, 0, 0), Ncandidate)
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
  call_mz_rt_freq_median_peak_property <- function(i) {
    x_h <- which(peak_property[i, 3:L_samples2] != 0)
    freq_h <- length(x_h)
    if (freq_h > 0) {
      m_h <- median(peak_property[i, (x_h + 2)])
    } else {
      m_h <- 0
    }
    c(freq_h, m_h)
  }
  ##
  ##############################################################################
  if (number_processing_threads == 1) {
    ##
    UFA_logRecorder("Initiated matching peak IDs!")
    progressBARboundaries <- txtProgressBar(min = 0, max = L_samples, initial = 0, style = 3)
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
    UFA_logRecorder("Completed matching peak IDs!")
    ##
    UFA_logRecorder("Initiated calculating median ranks!")
    progressBARboundaries <- txtProgressBar(min = 0, max = L_peaks, initial = 0, style = 3)
    #
    aligned_molecular_formula <- do.call(rbind, lapply(1:L_peaks, function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      call_calculating_median_ranks(i)
    }))
    MF_Zcol <- NULL
    close(progressBARboundaries)
    UFA_logRecorder("Completed calculating median ranks!")
    ##
    UFA_logRecorder("Initiated processing the peak property table!")
    progressBARboundaries <- txtProgressBar(min = 0, max = L_peaks, initial = 0, style = 3)
    IPA_Xcol <- do.call(rbind, lapply(1:L_peaks, function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      c(peakXcol[i, 1:2], length(which(peakXcol[i, 3:L_samples2] > 0)))
    }))
    close(progressBARboundaries)
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = L_peaks, initial = 0, style = 3)
    mz_rt_freq_median_peak_property <- do.call(rbind, lapply(1:L_peaks, function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      call_mz_rt_freq_median_peak_property(i)
    }))
    close(progressBARboundaries)
    ##
    title_mat <- do.call(cbind, lapply(1:Ncandidate, function(i) {
      cbind(paste0("IonFormula_", i), paste0("Frequency_", i), paste0("MedianRank_", i))
    }))
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      UFA_logRecorder("Initiated matching peak IDs!")
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
      UFA_logRecorder("Completed matching peak IDs!")
      ##
      UFA_logRecorder("Initiated calculating median ranks!")
      aligned_molecular_formula <- do.call(rbind, mclapply(1:L_peaks, function(i) {
        call_calculating_median_ranks(i)
      }, mc.cores = number_processing_threads))
      MF_Zcol <- NULL
      UFA_logRecorder("Completed calculating median ranks!")
      ##
      UFA_logRecorder("Initiated processing the peak property table!")
      ##
      IPA_Xcol <- do.call(rbind, mclapply(1:L_peaks, function(i) {
        c(peakXcol[i, 1:2], length(which(peakXcol[i, 3:L_samples2] > 0)))
      }, mc.cores = number_processing_threads))
      ##
      mz_rt_freq_median_peak_property <- do.call(rbind, mclapply(1:L_peaks, function(i) {
        call_mz_rt_freq_median_peak_property(i)
      }, mc.cores = number_processing_threads))
      ##
      title_mat <- do.call(cbind, mclapply(1:Ncandidate, function(i) {
        cbind(paste0("IonFormula_", i), paste0("Frequency_", i), paste0("MedianRank_", i))
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      UFA_logRecorder("Initiated matching peak IDs!")
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
      UFA_logRecorder("Completed matching peak IDs!")
      ##
      UFA_logRecorder("Initiated calculating median ranks!")
      aligned_molecular_formula <- foreach(i = 1:L_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_calculating_median_ranks(i)
      }
      MF_Zcol <- NULL
      UFA_logRecorder("Completed calculating median ranks!")
      ##
      UFA_logRecorder("Initiated processing the peak property table!")
      ##
      IPA_Xcol <- foreach(i = 1:L_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
        c(peakXcol[i, 1:2], length(which(peakXcol[i, 3:L_samples2] > 0)))
      }
      ##
      mz_rt_freq_median_peak_property <- foreach(i = 1:L_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_mz_rt_freq_median_peak_property(i)
      }
      ##
      title_mat <-  foreach(i = 1:Ncandidate, .combine = 'cbind', .verbose = FALSE) %dopar% {
        cbind(paste0("IonFormula_", i), paste0("Frequency_", i), paste0("MedianRank_", i))
      }
      ##
      stopCluster(clust)
    }
  }
  ##
  aligned_molecular_formula <- data.frame(cbind(IPA_Xcol, mz_rt_freq_median_peak_property, aligned_molecular_formula))
  ##
  ppn1 <- strsplit(input_path_peak_property, "/")[[1]]
  ppn <- ppn1[length(ppn1)]
  peak_property_name <- gsub(".Rdata$", "", ppn)
  #
  title_mat <- cbind("m/z", "RT", "IPA detection frequency", paste0(peak_property_name, " frequency"), paste0("median ", peak_property_name), title_mat)
  colnames(aligned_molecular_formula) <- title_mat
  rownames(aligned_molecular_formula) <- NULL
  UFA_logRecorder("Completed processing of the peak property table!")
  output_path_aligned_table <- paste0(output_path, "/aligned_molecular_formula_table")
  if (!dir.exists(output_path_aligned_table)) {
    dir.create(output_path_aligned_table, recursive = TRUE)
  }
  UFA_logRecorder("Initiated saving the aligned molecular formula annotated table!")
  save(aligned_molecular_formula, file = paste0(output_path_aligned_table, "/aligned_molecular_formula.Rdata"))
  write.csv(aligned_molecular_formula, file = paste0(output_path_aligned_table, "/aligned_molecular_formula.csv"))
  UFA_logRecorder("Stored annotated aligned table as `aligned_molecular_formula` in the `.Rdata` and `.csv` formats in the `aligned_molecular_formula_table` folder!")
  UFA_logRecorder("Completed generating data for the aligned molecular formula annotated table!")
  UFA_logRecorder(paste0(rep("", 100), collapse = "-"))
}
