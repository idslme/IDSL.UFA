score_coefficient_evaluation <- function(PARAM_SFT) {
  ##
  UFA_logRecorder(paste0(rep("", 100), collapse = "-"))
  UFA_logRecorder("Initiated evaluating score coefficients optimization!")
  ##
  output_path <- PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0010"), 2]
  output_path_score_function_calculations <- paste0(output_path, "/score_function_calculations")
  Entire_final_list_unoptimized <- IDSL.IPA::loadRdata(paste0(output_path_score_function_calculations, "/Entire_final_list_unoptimized.Rdata"))
  GA_score <- IDSL.IPA::loadRdata(paste0(output_path_score_function_calculations, "/GA_score.Rdata"))
  GA_score <- summary(GA_score)
  Score_coeff <- GA_score$solution
  Score_coeff <- Score_coeff[1, ]
  maxNEME <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0014"), 2])
  PCS <- as.numeric(Entire_final_list_unoptimized[, 11])
  RCS <- as.numeric(Entire_final_list_unoptimized[, 15])
  NEME <- as.numeric(Entire_final_list_unoptimized[, 10])
  R13C_PL <- as.numeric(Entire_final_list_unoptimized[, 12])
  R13C_IP <- as.numeric(Entire_final_list_unoptimized[, 13])
  size_IP <- as.numeric(Entire_final_list_unoptimized[, 9])
  N_compounds <- max(as.numeric(Entire_final_list_unoptimized$CompoundID))
  x_c <- lapply(1:N_compounds, function(i) {
    which(Entire_final_list_unoptimized$CompoundID == i)
  })
  N_candidate <- as.numeric(sapply(1:N_compounds, function(i) {
    Entire_final_list_unoptimized$CandidateCount[x_c[[i]][1]]
  }))
  IdentificationScore <- identification_score(Score_coeff, size_IP, PCS, RCS, NEME, maxNEME, R13C_PL, R13C_IP)
  Entire_final_list_unoptimized <- cbind(Entire_final_list_unoptimized, IdentificationScore)
  progressBARboundaries <- txtProgressBar(min = 0, max = N_compounds, initial = 0, style = 3)
  Entire_final_list_optimized <- do.call(rbind, lapply(1:N_compounds, function(i) {
    setTxtProgressBar(progressBARboundaries, i)
    A <- Entire_final_list_unoptimized[x_c[[i]], ]
    A <- A[order(A[, 20], decreasing = TRUE), ]
    A$Rank <- seq(1, N_candidate[i])
    A[, -20]
  }))
  close(progressBARboundaries)
  names(Entire_final_list_optimized) <- c("FileName", "PeakID",
                                          "ID_IonFormula", "IonFormula", "m/z Isotopic Profile",
                                          "m/z peaklist", "RT(min)", "PeakHeight", "size IP",
                                          "NEME(mDa)", "PCS", "R13C peakList", "R13C Isotopic Profile",
                                          "NDCS", "RCS(%)", "Rank", "CandidateCount", "CompoundID", "MolFMatch")
  rownames(Entire_final_list_optimized) <- NULL
  save(Entire_final_list_optimized, file = paste0(output_path_score_function_calculations, "/Entire_final_list_optimized.Rdata"))
  obj_function <- gsub(" ", "", tolower(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0018"), 2]))
  if (obj_function == "toprank") {
    max_rank <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0019"), 2])
    x <- which(as.numeric(Entire_final_list_unoptimized$MolFMatch) == 1)
    UFA_logRecorder(paste0("There detected totally ", length(x), " compounds for score coefficients optimization!!!"))
    r_unop <- length(which(as.numeric(Entire_final_list_unoptimized$Rank[x]) <= max_rank))
    UFA_logRecorder(paste0("There met ", r_unop, " peaks the <=", max_rank, " ranking with score coefficients of 1!!!"))
    x <- which(as.numeric(Entire_final_list_optimized$MolFMatch) == 1)
    r_op <- length(which(as.numeric(Entire_final_list_optimized$Rank[x]) <= max_rank))
    UFA_logRecorder(paste0("There met ", r_op, " peaks the <=", max_rank, " ranking after score coefficients optimization!!!"))
    R <- round((r_unop - r_op)/(r_unop - length(x)) * 100, 2)
  }
  if (obj_function == "overalrank") {
    x <- which(as.numeric(Entire_final_list_unoptimized$MolFMatch) == 1)
    NC <- as.numeric(Entire_final_list_optimized$CandidateCount[x])
    F_min <- sum(1/NC)
    UFA_logRecorder(paste0("The minimum achievable value of objective function is `", round(F_min, 2), "`!"))
    r_unop <- as.numeric(Entire_final_list_unoptimized$Rank[x])
    F_unop <- sum(r_unop/NC)
    UFA_logRecorder(paste0("The objective function was ", round(F_unop, 2), " with score coefficients of 1!!!"))
    x <- which(as.numeric(Entire_final_list_optimized$MolFMatch) == 1)
    r_op <- as.numeric(Entire_final_list_unoptimized$Rank[x])
    F_op <- sum(r_op/NC)
    UFA_logRecorder(paste0("The objective function became ", round(F_op, 2), " after score coefficients optimization!!!"))
    R <- round((F_unop - F_op)/(F_unop - F_min) * 100, 2)
  }
  UFA_logRecorder(paste0("The score coefficient optimization was ", R, "% successful with respect to score coefficients of 1 !!!"))
}
