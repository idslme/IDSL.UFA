scoreCoefficientsEvaluation <- function(PARAM_ScoreFunc) {
  ##
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Initiated evaluating score coefficients optimization!")
  ##
  output_path <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0008'), 2]
  output_path_score_function_calculations <- paste0(output_path, "/score_function_calculations")
  unoptimized_molecular_formula_annotation <- IDSL.IPA::loadRdata(paste0(output_path_score_function_calculations, "/unoptimized_molecular_formula_annotation.Rdata"))
  scoreCoefficients <- data.frame(read.csv(paste0(output_path_score_function_calculations, "/scoreCoefficients.csv"), header = TRUE))
  scoreCoefficients <- as.numeric(scoreCoefficients[1, ])
  ##
  PCS <- as.numeric(unoptimized_molecular_formula_annotation$'PCS')
  RCS <- as.numeric(unoptimized_molecular_formula_annotation$'ratioChromScan')
  NEME <- as.numeric(unoptimized_molecular_formula_annotation$'normEucMassError')
  R13C_PL <- as.numeric(unoptimized_molecular_formula_annotation$'r13cPeaklist')
  R13C_IP <- as.numeric(unoptimized_molecular_formula_annotation$'r13cTheoretical')
  size_IP <- as.numeric(unoptimized_molecular_formula_annotation$'isotopologueCount')
  ##
  unoptimized_molecular_formula_annotation$'CompoundID' <- as.numeric(unoptimized_molecular_formula_annotation$'CompoundID')
  N_compounds <- max(unoptimized_molecular_formula_annotation$'CompoundID')
  updatedIdentificationScore <- identificationScoreCalculator(scoreCoefficients, size_IP, PCS, RCS, NEME, R13C_PL, R13C_IP)
  xDiff <- c(0, which(abs(diff(unoptimized_molecular_formula_annotation$'CompoundID')) > 0), nrow(unoptimized_molecular_formula_annotation))
  ##
  progressBARboundaries <- txtProgressBar(min = 0, max = N_compounds, initial = 0, style = 3)
  ##
  optimized_molecular_formula_annotation <- do.call(rbind, lapply(1:(length(xDiff) - 1), function(i) {
    setTxtProgressBar(progressBARboundaries, i)
    ##
    xCompound <- seq((xDiff[i] + 1), xDiff[i + 1], 1)
    nCandidates <- as.numeric(unoptimized_molecular_formula_annotation$'CandidateCount'[xCompound[1]])
    A <- unoptimized_molecular_formula_annotation[xCompound, ]
    if (nCandidates > 1) {
      orderA <- order(updatedIdentificationScore[xCompound], decreasing = TRUE)
      A <- A[orderA, ]
      A$Rank <- seq(1, nCandidates, 1)
    }
    ##
    return(A)
  }))
  close(progressBARboundaries)
  ##
  optimized_molecular_formula_annotation <- data.frame(optimized_molecular_formula_annotation)
  rownames(optimized_molecular_formula_annotation) <- NULL
  save(optimized_molecular_formula_annotation, file = paste0(output_path_score_function_calculations, "/optimized_molecular_formula_annotation.Rdata"))
  ##
  ##############################################################################
  ##
  obj_function <- gsub(" ", "", tolower(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0016'), 2]))
  if (obj_function == "toprank") {
    maxRank <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0017'), 2])
    x <- which(as.numeric(unoptimized_molecular_formula_annotation$'MolFMatch') == 1)
    IPA_logRecorder(paste0("There detected totally `", N_compounds, "` compounds for score coefficients optimization!!!"))
    r_unop <- length(which(as.numeric(unoptimized_molecular_formula_annotation$'rank'[x]) <= maxRank))
    IPA_logRecorder(paste0("There met ", r_unop, " peaks the `<= ", maxRank, "` ranking with score coefficients of 1!!!"))
    x <- which(as.numeric(optimized_molecular_formula_annotation$MolFMatch) == 1)
    r_op <- length(which(as.numeric(optimized_molecular_formula_annotation$Rank[x]) <= maxRank))
    IPA_logRecorder(paste0("There met ", r_op, " peaks the `<= ", maxRank, "` ranking after score coefficients optimization!!!"))
    R <- round((r_unop - r_op)/(r_unop - N_compounds) * 100, 2)
    ##
  } else if (obj_function == "overalrank") {
    ##
    x <- which(as.numeric(unoptimized_molecular_formula_annotation$'MolFMatch') == 1)
    NC <- as.numeric(optimized_molecular_formula_annotation$'CandidateCount'[x])
    F_min <- sum(1/NC)
    IPA_logRecorder(paste0("The minimum achievable value of objective function is `", round(F_min, 2), "`!"))
    r_unop <- as.numeric(unoptimized_molecular_formula_annotation$'rank'[x])
    F_unop <- sum(r_unop/NC)
    IPA_logRecorder(paste0("The objective function was `", round(F_unop, 2), "` with score coefficients of 1!!!"))
    x <- which(as.numeric(optimized_molecular_formula_annotation$'MolFMatch') == 1)
    r_op <- as.numeric(unoptimized_molecular_formula_annotation$'rank'[x])
    F_op <- sum(r_op/NC)
    IPA_logRecorder(paste0("The objective function became `", round(F_op, 2), "` after score coefficients optimization!!!"))
    R <- round((F_unop - F_op)/(F_unop - F_min) * 100, 2)
  }
  IPA_logRecorder(paste0("The score coefficient optimization was `", R, "%` successful with respect to score coefficients of 1 !!!"))
  ##
  return()
}
