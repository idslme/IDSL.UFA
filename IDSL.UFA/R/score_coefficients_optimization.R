score_coefficients_optimization <- function(PARAM_SFT) {
  UFA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  number_processing_threads <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0005"), 2])
  output_path <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0010'), 2]
  output_path_score_function_calculations <- paste0(output_path, "/score_function_calculations")
  Entire_final_list_unoptimized <- IDSL.IPA::loadRdata(paste0(output_path_score_function_calculations, "/Entire_final_list_unoptimized.Rdata"))
  maxNEME <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0014"), 2])
  PCS <- as.numeric(Entire_final_list_unoptimized[, 11])
  RCS <- as.numeric(Entire_final_list_unoptimized[, 15])
  NEME <- as.numeric(Entire_final_list_unoptimized[, 10])
  R13C_PL <- as.numeric(Entire_final_list_unoptimized[, 12])
  R13C_IP <- as.numeric(Entire_final_list_unoptimized[, 13])
  size_IP <- as.numeric(Entire_final_list_unoptimized[, 9])
  ##
  N_compounds <- max(as.numeric(Entire_final_list_unoptimized$CompoundID))
  x_c <- lapply(1:N_compounds, function(i) {
    which(Entire_final_list_unoptimized$CompoundID == i)
  })
  ##
  x_true_matchedcompound <- do.call(rbind, lapply(1:N_compounds, function(i) {
    x_true <- 0
    x_t <- which(Entire_final_list_unoptimized$MolFMatch[x_c[[i]]] == 1)
    if (length(x_t) > 0) {
      x_true <- x_t[1]
    }
    c(i, x_true)
  }))
  x_true_list <- which(x_true_matchedcompound[, 2] > 0)
  L_x_true_list <- length(x_true_list)
  if (L_x_true_list > 0) {
    UFA_logRecorder(paste0("There are totally ", L_x_true_list, " compounds matched for the score coefficients optimization!"))
    matchedcompounds <- x_true_matchedcompound[x_true_list, 1]
    x_true_false <- x_true_matchedcompound[, 2]
    ##
    N_candidate <- as.numeric(sapply(1:N_compounds, function(i) {
      Entire_final_list_unoptimized$CandidateCount[x_c[[i]][1]]
    }))
    ##
    obj_function <- gsub(" ", "", tolower(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0018'), 2]))
    UFA_logRecorder(paste0("The objective function method `", obj_function, "` is used for the score coefficients optimization!"))
    #
    if (obj_function == "toprank") {
      max_rank <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0019"), 2])
      obj_ga <- function(Score_coeff) {
        IdentificationScore <- identification_score(Score_coeff, size_IP, PCS, RCS, NEME, maxNEME, R13C_PL, R13C_IP)
        ObjF <- sum(sapply(matchedcompounds, function(i) {
          R2N <- 0
          x_true <- x_true_false[i]
          x_candidate <- rep(0, N_candidate[i])
          IdS <- IdentificationScore[x_c[[i]]]
          x_candidate[x_true] <- 1
          sd <- matrix(cbind(IdS, x_candidate), ncol = 2)
          sd <- matrix(sd[order(sd[, 1], decreasing = TRUE), ], ncol = 2)
          x_1 <- which(sd[, 2] == 1)[1]
          if (x_1 <= max_rank) {
            R2N <- 1
          }
          R2N
        }))
        return(ObjF)
      }
    }
    #
    if (obj_function == "overalrank") {
      obj_ga <- function(Score_coeff) {
        IdentificationScore <- identification_score(Score_coeff, size_IP, PCS, RCS, NEME, maxNEME, R13C_PL, R13C_IP)
        ObjF <- sum(sapply(matchedcompounds, function(i) {
          x_true <- x_true_false[i]
          x_candidate <- rep(0, N_candidate[i])
          IdS <- IdentificationScore[x_c[[i]]]
          x_candidate[x_true] <- 1
          sd <- matrix(cbind(IdS, x_candidate), ncol = 2)
          sd <- matrix(sd[order(sd[, 1], decreasing = TRUE), ], ncol = 2)
          x_1 <- which(sd[, 2] == 1)[1]
          x_1/N_candidate[i]
        }))
        return(-ObjF)
      }
    }
    ######## GA ########
    lower_limit <- eval(parse(text = PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0020"), 2]))
    upper_limit <- eval(parse(text = PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0021"), 2]))
    population_size <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0022"), 2])
    max_iteration <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0023"), 2])
    ## To clear cache memory
    k <- 6
    while (k > 0) {
      Sys.sleep(k)
      closeAllConnections()
      gc()
      k <- k - 1
    }
    ##
    UFA_logRecorder("Initiated the genetic algorithm optimization! This step may take from a minute to several hours!")
    GA_score <- GA::ga(type = "real-valued", fitness = obj_ga, lower = lower_limit, upper = upper_limit, popSize = population_size, maxiter = max_iteration, parallel = number_processing_threads)
    save(GA_score, file = paste0(output_path_score_function_calculations, "/GA_score.Rdata"))
    GA_score <- summary(GA_score)
    Score_coeff <- GA_score$solution
    write.csv(Score_coeff, file = paste0(output_path_score_function_calculations, "/score_coefficients.csv"))
    UFA_logRecorder("Stored the genetic algorithm results in 'GA_score.Rdata'")
    UFA_logRecorder("Stored score coefficients as `score_coefficients.csv` in the `score_function_calculations` folder!")
  } else {
    stop(UFA_logRecorder("Error!!! No molecular formula in the reference list have been matched to optimize score coefficients!"))
  }
}
