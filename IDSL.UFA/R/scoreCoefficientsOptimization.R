scoreCoefficientsOptimization <- function(PARAM_ScoreFunc) {
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  output_path <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0008'), 2]
  output_path_score_function_calculations <- paste0(output_path, "/score_function_calculations")
  unoptimized_molecular_formula_annotation <- IDSL.IPA::loadRdata(paste0(output_path_score_function_calculations, "/unoptimized_molecular_formula_annotation.Rdata"))
  PCS <- as.numeric(unoptimized_molecular_formula_annotation$'PCS')
  RCS <- as.numeric(unoptimized_molecular_formula_annotation$'ratioChromScan')
  NEME <- as.numeric(unoptimized_molecular_formula_annotation$'normEucMassError')
  R13C_PL <- as.numeric(unoptimized_molecular_formula_annotation$'r13cPeaklist')
  R13C_IP <- as.numeric(unoptimized_molecular_formula_annotation$'r13cTheoretical')
  size_IP <- as.numeric(unoptimized_molecular_formula_annotation$'isotopologueCount')
  ##
  N_compounds <- max(as.numeric(unoptimized_molecular_formula_annotation$'CompoundID'))
  if (N_compounds > 0) {
    x_c <- base::tapply(seq(1, nrow(unoptimized_molecular_formula_annotation), 1), unoptimized_molecular_formula_annotation$'CompoundID', 'c', simplify = FALSE)
    ##
    x_true_matchedcompound <- do.call(rbind, lapply(1:N_compounds, function(i) {
      x_true <- 0
      x_t <- which(unoptimized_molecular_formula_annotation$'MolFMatch'[x_c[[i]]] == 1)
      if (length(x_t) > 0) {
        x_true <- x_t[1]
      }
      c(i, x_true)
    }))
    x_true_list <- which(x_true_matchedcompound[, 2] > 0)
    L_x_true_list <- length(x_true_list)
    if (L_x_true_list > 0) {
      IPA_logRecorder(paste0("There are totally `", L_x_true_list, "` compounds matched for the score coefficients optimization!"))
      matchedcompounds <- x_true_matchedcompound[x_true_list, 1]
      x_true_false <- x_true_matchedcompound[, 2]
      ##
      N_candidate <- as.numeric(sapply(1:N_compounds, function(i) {
        unoptimized_molecular_formula_annotation$'CandidateCount'[x_c[[i]][1]]
      }))
      ##
      obj_function <- gsub(" ", "", tolower(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0016'), 2]))
      IPA_logRecorder(paste0("The objective function method of `", obj_function, "` is used for the score coefficients optimization!"))
      #
      if (obj_function == "toprank") {
        maxRank <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0017'), 2])
        obj_ga <- function(scoreCoefficients) {
          IdentificationScores <- identificationScoreCalculator(scoreCoefficients, size_IP, PCS, RCS, NEME, R13C_PL, R13C_IP)
          ObjF <- do.call(sum, lapply(matchedcompounds, function(i) {
            R2N <- 0
            x_true <- x_true_false[i]
            x_candidate <- rep(0, N_candidate[i])
            IdS <- IdentificationScores[x_c[[i]]]
            x_candidate[x_true] <- 1
            x_candidate <- x_candidate[order(IdS, decreasing = TRUE)]
            x_1 <- which(x_candidate == 1)[1]
            if (x_1 <= maxRank) {
              R2N <- 1
            }
            R2N
          }))
          return(ObjF)
        }
      }
      #
      if (obj_function == "overalrank") {
        obj_ga <- function(scoreCoefficients) {
          IdentificationScores <- identificationScoreCalculator(scoreCoefficients, size_IP, PCS, RCS, NEME, R13C_PL, R13C_IP)
          ObjF <- do.call(sum, lapply(matchedcompounds, function(i) {
            x_true <- x_true_false[i]
            x_candidate <- rep(0, N_candidate[i])
            IdS <- IdentificationScores[x_c[[i]]]
            x_candidate[x_true] <- 1
            x_candidate <- x_candidate[order(IdS, decreasing = TRUE)]
            x_1 <- which(x_candidate == 1)[1]
            x_1/N_candidate[i]
          }))
          return(-ObjF)
        }
      }
      ######## GA ########
      number_processing_threads <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0018'), 2])
      lower_limit <- eval(parse(text = PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0019'), 2]))
      upper_limit <- eval(parse(text = PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0020'), 2]))
      population_size <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0021'), 2])
      max_iteration <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0022'), 2])
      ## To clear cache memory
      k <- 6
      while (k > 0) {
        Sys.sleep(k)
        closeAllConnections()
        gc()
        k <- k - 1
      }
      ##
      GApackageCheck <- tryCatch(requireNamespace('GA', quietly = TRUE), error = function(e) {FALSE})
      if (!GApackageCheck) {
        IPA_logRecorder("IDSL.UFA requires the 'GA' package of R for the score coefficients optimization workflow!")
        stop(IPA_logRecorder(" <<< install.packages('GA') >>> "))
      }
      ##
      IPA_logRecorder("Initiated the genetic algorithm optimization! This step may take from a minute to several hours!")
      GA_score <- GA::ga(type = "real-valued", fitness = obj_ga, lower = lower_limit, upper = upper_limit, popSize = population_size, maxiter = max_iteration, parallel = number_processing_threads)
      save(GA_score, file = paste0(output_path_score_function_calculations, "/GA_score.Rdata"))
      GA_score <- summary(GA_score)
      scoreCoefficients <- GA_score$solution
      write.csv(scoreCoefficients, file = paste0(output_path_score_function_calculations, "/scoreCoefficients.csv"), row.names = FALSE)
      IPA_logRecorder("Stored the genetic algorithm results in `GA_score.Rdata`")
      IPA_logRecorder("Stored score coefficients as `scoreCoefficients.csv` in the `score_function_calculations` folder!")
    } else {
      stop(IPA_logRecorder("Error!!! No molecular formula in the reference list have been matched to optimize score coefficients!"))
    }
  } else {
    stop(IPA_logRecorder("Error!!! No molecular formula in the reference list have been matched to optimize score coefficients!"))
  }
  ##
  return()
}
