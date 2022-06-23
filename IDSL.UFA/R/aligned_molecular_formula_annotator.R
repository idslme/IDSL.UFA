aligned_molecular_formula_annotator <- function(PARAM) {
  print("Initiated creating the aligned molecular formula annotated table!")
  number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == "PARAM0009"), 2])
  ##
  input_path_peak_Xcol <- PARAM[which(PARAM[, 1] == 'PARAM0025'), 2]
  peak_Xcol <- loadRdata(input_path_peak_Xcol)
  L_peaks <- dim(peak_Xcol)[1]
  L_samples <- dim(peak_Xcol)[2] - 2
  ##
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0014'), 2]
  output_path_annotated_mf_tables <- paste0(output_path, "/annotated_mf_tables")
  mf_table_list <- dir(path = output_path_annotated_mf_tables, pattern = ".Rdata")
  if (length(mf_table_list) != L_samples) {
    AnPL <- gsub("MolecularFormulaAnnotationTable_", "", mf_table_list)
    ColPL <- colnames(peak_Xcol)[3:(L_samples + 2)]
    MissedPL <- setdiff(ColPL, AnPL)
    print("Error!!! The following MolecularFormulaAnnotationTables are not avialable:")
    for (i in 1:length(MissedPL)) {
      print(MissedPL[i])
    }
    stop()
  }
  ##
  input_path_peak_property <- PARAM[which(PARAM[, 1] == 'PARAM0026'), 2]
  peak_property <- loadRdata(input_path_peak_property)
  if (dim(peak_property)[1] != L_peaks | dim(peak_property)[2] != (L_samples + 2)) {
    stop("Error!!! aligned peak property and peak indexed tables are not in the same size!!!")
  }
  ##
  maxRankSample <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0027'), 2])
  N_candidate <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0028'), 2])
  ##
  adjust_freq_rank <- PARAM[which(PARAM[, 1] == 'PARAM0029'), 2]
  N_candidate3 <- N_candidate*3
  ##
  call_MF_Zcol <- function(i) {
    peak_table_id <- peak_Xcol[, (i + 2)]
    MolecularFormulaAnnotationTable <- loadRdata(paste0(output_path_annotated_mf_tables, "/", mf_table_list[i]))
    matched_peak_ids <- as.numeric(MolecularFormulaAnnotationTable[, 1])
    x_peak_ids <- which(peak_table_id %in% unique(matched_peak_ids) == TRUE)
    ##
    if (length(x_peak_ids) > 0) {
      matched_mf_ids <- as.numeric(MolecularFormulaAnnotationTable[, 2])
      ##
      do.call(rbind, lapply(x_peak_ids, function(j) {
        x_j <- which(matched_peak_ids == peak_table_id[j])
        max_k <- min(c(maxRankSample, length(x_j)))
        ##
        cbind(rep(j, max_k), matched_mf_ids[x_j[1:max_k]], seq(1, max_k, 1))
      }))
    }
  }
  ##
  call_calculating_median_ranks <- function(j) {
    ID_freq_Rank <- rep(0, N_candidate3)
    ##
    if (xZcol[j, 1] != 0) {
      ID_rank <- matrix(MF_Zcol[xZcol[j, 1]:xZcol[j, 2], 2:3], ncol = 2)
      ##
      t_IDs <- sort(table(ID_rank[, 1]), decreasing = TRUE)
      max_k <- min(c(N_candidate, length(t_IDs)))
      t_freq <- as.numeric(t_IDs[1:max_k])
      t_id <- as.numeric(names(t_IDs[1:max_k]))
      ID_freq_Rank3 <- do.call(rbind, lapply(1:max_k, function(i) {
        x_id_t <- which(ID_rank[, 1] == t_id[i])
        med_rank <- median(ID_rank[x_id_t, 2]) # A median is calculated for the rank of candidate compounds across samples
        c(t_id[i], t_freq[i], med_rank)
      }))
      ID_freq_Rank3 <- matrix(ID_freq_Rank3[order(ID_freq_Rank3[, 3], decreasing = FALSE), ], ncol = 3)
      ID_freq_Rank3 <- matrix(ID_freq_Rank3[order(ID_freq_Rank3[, 2], decreasing = TRUE), ], ncol = 3)
      ##
      for (i in 1:max_k) {
        ID_freq_Rank[3*i - 2] <- ID_freq_Rank3[i, 1]
        ID_freq_Rank[3*i - 1] <- ID_freq_Rank3[i, 2]
        ID_freq_Rank[3*i] <- ID_freq_Rank3[i, 3]
      }
    }
    return(ID_freq_Rank)
  }
  ##
  call_creating_aligned_table <- function(i) {
    molf_IDs <- M_IDs[, (3*i - 2)]
    x_non0 <- which(molf_IDs != 0)
    molf <- rep(NA, L_peaks)
    if (length(x_non0) > 0) {
      matched_IDs_vec_db <- matrix(MolVecList_DB[molf_IDs[x_non0], ], ncol = L_Elements)
      molf_hill <- hill_molecular_formula_printer(Elements, matched_IDs_vec_db)
      molf[x_non0] <- molf_hill
    }
    sub_table <- cbind(molf, M_IDs[, (3*i - 1)], M_IDs[, 3*i])
    return(sub_table)
  }
  ##
  call_mz_rt_freq_median_peak_property <- function(i) {
    m_h <- 0
    x_h <- which(peak_property[i, 3:(L_samples + 2)] != 0)
    freq_h <- length(x_h)
    if (freq_h > 0) {
      m_h <- median(peak_property[i, (x_h + 2)])
    }
    c(freq_h, m_h)
  }
  ##
  if (tolower(adjust_freq_rank) == "yes") {
    ##
    call_freq_rank_table <- function(k) {
      freq <- M_IDs[, (3*k - 1)]
      rank <- M_IDs[, (3*k)]
      sqrt(freq)/rank           # To adjust ranking and frequencies
    }
    ##
    call_M_IDs2 <- function(k) {
      x_0 <- which(freq_rank_table[k, ] > 0)
      L_x_0 <- length(x_0)
      if (L_x_0 > 1) {
        order_rank <- order(freq_rank_table[k, x_0], decreasing = TRUE)
        M_ID_ordered <- do.call(cbind, lapply(order_rank, function (i) {
          cbind(M_IDs[k, (3*i - 2)], M_IDs[k, (3*i - 1)], M_IDs[k, 3*i])
        }))
        if (L_x_0 < N_candidate) {
          M_ID_ordered <- c(M_ID_ordered, rep(0, (N_candidate3 - 3*L_x_0)))
        }
      } else {
        M_ID_ordered <- M_IDs[k, ]
      }
      return(M_ID_ordered)
    }
  }
  ####
  if (number_processing_threads == 1) {
    ##
    print("Initiated matching peak IDs!")
    progressBARboundaries <- txtProgressBar(min = 0, max = L_samples, initial = 0, style = 3)
    #
    MF_Zcol <- do.call(rbind, lapply(1:L_samples, function(k) {
      setTxtProgressBar(progressBARboundaries, k)
      ##
      call_MF_Zcol(k)
    }))
    close(progressBARboundaries)
    #
    MF_Zcol <- MF_Zcol[order(MF_Zcol[, 1], decreasing = FALSE), ]
    xDiff <- which(diff(MF_Zcol[, 1]) > 0)
    #
    xZcol <- matrix(rep(0, 2*L_peaks), ncol = 2)
    #
    u_peakid <- unique(MF_Zcol[, 1])
    xZcol[u_peakid, 1] <- c(1, (xDiff + 1))
    xZcol[u_peakid, 2] <- c(xDiff, dim(MF_Zcol)[1])
    #
    print("Completed matching peak IDs!")
    ##
    print("Initiated calculating median ranks!")
    progressBARboundaries <- txtProgressBar(min = 0, max = L_peaks, initial = 0, style = 3)
    #
    M_IDs <- do.call(rbind, lapply(1:L_peaks, function(k) {
      setTxtProgressBar(progressBARboundaries, k)
      ##
      call_calculating_median_ranks(k)
    }))
    MF_Zcol <- NULL
    close(progressBARboundaries)
    print("Completed calculating median ranks!")
    ##
    if (tolower(adjust_freq_rank) == "yes") {
      ##
      print("Initiated adjusting frequencies and ranks!")
      progressBARboundaries <- txtProgressBar(min = 0, max = N_candidate, initial = 0, style = 3)
      #
      freq_rank_table <- do.call(cbind, lapply(1:N_candidate, function(k) {
        setTxtProgressBar(progressBARboundaries, k)
        ##
        call_freq_rank_table(k)
      }))
      close(progressBARboundaries)
      #
      freq_rank_table[is.nan(freq_rank_table)] <- 0
      ##
      progressBARboundaries <- txtProgressBar(min = 0, max = L_peaks, initial = 0, style = 3)
      #
      M_IDs <- do.call(rbind, lapply(1:L_peaks, function(k) {
        setTxtProgressBar(progressBARboundaries, k)
        ##
        call_M_IDs2(k)
      }))
      close(progressBARboundaries)
      print("Completed adjusting frequencies and ranks!")
    }
    ##
    address_sav_IPDB <- PARAM[which(PARAM[, 1] == "PARAM0004"), 2]
    print("Loading the isotopic profiles database!")
    IPDB <- loadRdata(address_sav_IPDB)
    MolVecList0 <- IPDB[[2]]
    IPDB <- NULL
    Elements <- MolVecList0[[1]]
    MolVecList_DB <- MolVecList0[[2]]
    L_Elements <- length(Elements)
    print("Initiated creating the aligned table!")
    progressBARboundaries <- txtProgressBar(min = 0, max = N_candidate, initial = 0, style = 3)
    #
    aligned_molecular_formula <- do.call(cbind, lapply(1:N_candidate, function(k) {
      setTxtProgressBar(progressBARboundaries, k)
      ##
      call_creating_aligned_table(k)
    }))
    close(progressBARboundaries)
    print("Completed creating the aligned table!")
    ##
    print("Initiated processing the peak property table!")
    progressBARboundaries <- txtProgressBar(min = 0, max = L_peaks, initial = 0, style = 3)
    IPA_Xcol <- do.call(rbind, lapply(1:L_peaks, function(k) {
      setTxtProgressBar(progressBARboundaries, k)
      ##
      c(peak_Xcol[k, 1:2], (length(which(peak_Xcol[k, ] > 0)) - 2))
    }))
    close(progressBARboundaries)
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = L_peaks, initial = 0, style = 3)
    mz_rt_freq_median_peak_property <- do.call(rbind, lapply(1:L_peaks, function(k) {
      setTxtProgressBar(progressBARboundaries, k)
      ##
      call_mz_rt_freq_median_peak_property(k)
    }))
    close(progressBARboundaries)
    ##
    title_mat <- do.call(cbind, lapply(1:N_candidate, function(k) {
      cbind(paste0("IonFormula_", k), paste0("Frequency_", k), paste0("MedianRank_", k))
    }))
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      print("Initiated matching peak IDs!")
      MF_Zcol <- do.call(rbind, mclapply(1:L_samples, function(k) {
        call_MF_Zcol(k)
      }, mc.cores = number_processing_threads))
      #
      MF_Zcol <- MF_Zcol[order(MF_Zcol[, 1], decreasing = FALSE), ]
      xDiff <- which(diff(MF_Zcol[, 1]) > 0)
      #
      xZcol <- matrix(rep(0, 2*L_peaks), ncol = 2)
      #
      u_peakid <- unique(MF_Zcol[, 1])
      xZcol[u_peakid, 1] <- c(1, (xDiff + 1))
      xZcol[u_peakid, 2] <- c(xDiff, dim(MF_Zcol)[1])
      #
      print("Completed matching peak IDs!")
      ##
      print("Initiated calculating median ranks!")
      M_IDs <- do.call(rbind, mclapply(1:L_peaks, function(k) {
        call_calculating_median_ranks(k)
      }, mc.cores = number_processing_threads))
      MF_Zcol <- NULL
      print("Completed calculating median ranks!")
      ##
      if (tolower(adjust_freq_rank) == "yes") {
        print("Initiated adjusting frequencies and ranks!")
        ##
        freq_rank_table <- do.call(cbind, mclapply(1:N_candidate, function(k) {
          call_freq_rank_table(k)
        }, mc.cores = number_processing_threads))
        #
        freq_rank_table[is.nan(freq_rank_table)] <- 0
        ##
        M_IDs <- do.call(rbind, mclapply(1:L_peaks, function(k) {
          call_M_IDs2(k)
        }, mc.cores = number_processing_threads))
        print("Completed adjusting frequencies and ranks!")
      }
      ##
      address_sav_IPDB <- PARAM[which(PARAM[, 1] == "PARAM0004"), 2]
      print("Loading the isotopic profiles database!")
      IPDB <- loadRdata(address_sav_IPDB)
      MolVecList0 <- IPDB[[2]]
      IPDB <- NULL
      Elements <- MolVecList0[[1]]
      MolVecList_DB <- MolVecList0[[2]]
      L_Elements <- length(Elements)
      print("Initiated creating the aligned table!")
      aligned_molecular_formula <- do.call(cbind, mclapply(1:N_candidate, function(k) {
        call_creating_aligned_table(k)
      }, mc.cores = number_processing_threads))
      print("Completed creating the aligned table!")
      ##
      print("Initiated processing the peak property table!")
      ##
      IPA_Xcol <- do.call(rbind, mclapply(1:L_peaks, function(k) {
        c(peak_Xcol[k, 1:2], (length(which(peak_Xcol[k, ] > 0)) - 2))
      }, mc.cores = number_processing_threads))
      ##
      mz_rt_freq_median_peak_property <- do.call(rbind, mclapply(1:L_peaks, function(k) {
        call_mz_rt_freq_median_peak_property(k)
      }, mc.cores = number_processing_threads))
      ##
      title_mat <- do.call(cbind, mclapply(1:N_candidate, function(k) {
        cbind(paste0("IonFormula_", k), paste0("Frequency_", k), paste0("MedianRank_", k))
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      print("Initiated matching peak IDs!")
      MF_Zcol <- foreach(k = 1:L_samples, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_MF_Zcol(k)
      }
      #
      MF_Zcol <- MF_Zcol[order(MF_Zcol[, 1], decreasing = FALSE), ]
      xDiff <- which(diff(MF_Zcol[, 1]) > 0)
      #
      xZcol <- matrix(rep(0, 2*L_peaks), ncol = 2)
      #
      u_peakid <- unique(MF_Zcol[, 1])
      xZcol[u_peakid, 1] <- c(1, (xDiff + 1))
      xZcol[u_peakid, 2] <- c(xDiff, dim(MF_Zcol)[1])
      #
      print("Completed matching peak IDs!")
      ##
      print("Initiated calculating median ranks!")
      M_IDs <- foreach(k = 1:L_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_calculating_median_ranks(k)
      }
      MF_Zcol <- NULL
      print("Completed calculating median ranks!")
      ##
      if (tolower(adjust_freq_rank) == "yes") {
        print("Initiated adjusting frequencies and ranks!")
        ##
        freq_rank_table <- foreach(k = 1:N_candidate, .combine = 'cbind', .verbose = FALSE) %dopar% {
          call_freq_rank_table(k)
        }
        #
        freq_rank_table[is.nan(freq_rank_table)] <- 0
        ##
        M_IDs <- foreach(k = 1:L_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
          call_M_IDs2(k)
        }
        print("Completed adjusting frequencies and ranks!")
      }
      ##
      address_sav_IPDB <- PARAM[which(PARAM[, 1] == "PARAM0004"), 2]
      print("Loading the isotopic profiles database!")
      IPDB <- loadRdata(address_sav_IPDB)
      MolVecList0 <- IPDB[[2]]
      IPDB <- NULL
      Elements <- MolVecList0[[1]]
      MolVecList_DB <- MolVecList0[[2]]
      L_Elements <- length(Elements)
      print("Initiated creating the aligned table!")
      aligned_molecular_formula <- foreach(k = 1:N_candidate, .combine = 'cbind', .verbose = FALSE) %dopar% {
        call_creating_aligned_table(k)
      }
      print("Completed creating the aligned table!")
      ##
      print("Initiated processing the peak property table!")
      ##
      IPA_Xcol <- foreach(k = 1:L_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
        c(peak_Xcol[k, 1:2], (length(which(peak_Xcol[k, ] > 0)) - 2))
      }
      ##
      mz_rt_freq_median_peak_property <- foreach(k = 1:L_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_mz_rt_freq_median_peak_property(k)
      }
      ##
      title_mat <-  foreach(k = 1:N_candidate, .combine = 'cbind', .verbose = FALSE) %dopar% {
        cbind(paste0("IonFormula_", k), paste0("Frequency_", k), paste0("MedianRank_", k))
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
  peak_property_name <- gsub(".Rdata", "", ppn)
  #
  title_mat <- cbind("m/z", "RT", "IPA detection frequency", paste0(peak_property_name, " frequency"), paste0("median ", peak_property_name), title_mat)
  colnames(aligned_molecular_formula) <- title_mat
  rownames(aligned_molecular_formula) <- c()
  print("Completed processing of the peak property table!")
  output_path_aligned_table <- paste0(output_path, "/aligned_molecular_formula_table")
  if (!dir.exists(output_path_aligned_table)) {
    dir.create(output_path_aligned_table)
  }
  print("Initiated saving the aligned molecular formula table!")
  save(aligned_molecular_formula, file = paste0(output_path_aligned_table, "/aligned_molecular_formula.Rdata"))
  write.csv(aligned_molecular_formula, file = paste0(output_path_aligned_table, "/aligned_molecular_formula.csv"))
  print("Completed saving the aligned molecular formula table!")
}
