molecularFormula2IPdb <- function(molecularFormulaDatabase, retentionTime = NULL, peak_spacing = 0, intensity_cutoff_str = 1,
                                  IonPathways = "[M]+", number_processing_threads = 1, UFA_IP_memeory_variables = c(1e30, 1e-12, 100),
                                  allowedMustRunCalculation = FALSE, allowedVerbose = TRUE) {
  ##
  ##############################################################################
  ##############################################################################
  ##
  IsoProfInf100 <- matrix(c(Inf, 100), ncol = 2)
  ##
  failed_IP_calculator <- function(i) {
    if (allowedVerbose) {
      failedFormula <- hill_molecular_formula_printer(Elements, MoleFormVecMat[i, ], number_processing_threads = 1)
      IPA_logRecorder(paste0("Failed to calculate isotopic profile for `", failedFormula,"`!"), allowedPrinting = FALSE)
    }
    ##
    return(IsoProfInf100)
  }
  ##
  if (number_processing_threads > 1) {
    osType <- Sys.info()[['sysname']]
  }
  ##
  if (is.null(retentionTime)) {
    retentionTimeCheck <- FALSE
  } else {
    retentionTimeCheck <- TRUE
    retentionTime <- suppressWarnings(as.numeric(retentionTime))
    retentionTime[is.na(retentionTime)] <- Inf
  }
  ##
  if (typeof(molecularFormulaDatabase) == "list") {
    ##
    Elements <- molecularFormulaDatabase[[1]]
    MoleFormVecMat <- molecularFormulaDatabase[[2]]
    molVecCharge <- molecularFormulaDatabase[[3]] ## Charge of the molecular ions to divide the mass of profile by
    nMoleFormVecMat <- dim(MoleFormVecMat)[1]
    molecularFormulaDatabase <- NULL
    formulaDeconvolutionCheck <- FALSE
  } else {
    formulaDeconvolutionCheck <- TRUE
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  IP_calculator <- "IP_calculator <- function(i) {
    c <- MoleFormVecMat[i, x_el_c]
    if (elementBcheck) {
      b <- MoleFormVecMat[i, x_el_b]
    } else {
      b <- 0
    }
    if (elementBrCheck) {
      br <- MoleFormVecMat[i, x_el_br]
    } else {
      br <- 0
    }
    if (elementClcheck) {
      cl <- MoleFormVecMat[i, x_el_cl]
    } else {
      cl <- 0
    }
    if (elementKcheck) {
      k <- MoleFormVecMat[i, x_el_k]
    } else {
      k <- 0
    }
    if (elementScheck) {
      s <- MoleFormVecMat[i, x_el_s]
    } else {
      s <- 0
    }
    if (elementSeCheck) {
      se <- MoleFormVecMat[i, x_el_se]
    } else {
      se <- 0
    }
    if (elementSiCheck) {
      si <- MoleFormVecMat[i, x_el_si]
    } else {
      si <- 0
    }
    ##
    intensity_cutoff <- intensity_cutoff_str
    ##
    IPP <- tryCatch(isotopic_profile_calculator(MoleFormVecMat[i, ], massAbundanceList, peak_spacing, intensity_cutoff, UFA_IP_memeory_variables),
                    error = function(e) {failed_IP_calculator(i)},
                    warning = function(w) {failed_IP_calculator(i)})
    ##
    IPP[, 1] <- round(IPP[, 1]/molVecCharge[i], 6)
    IPP[, 2] <- round(IPP[, 2], 3)
    ##
    return(IPP)
  }"
  IP_calculator <- gsub("intensity_cutoff_str", intensity_cutoff_str, IP_calculator)
  eval(parse(text = IP_calculator))
  ##
  ##############################################################################
  ##############################################################################
  ##
  parametersIPDBcalculator <- function(i) {
    IPP <- IsotopicProfileList[[i]]
    x100 <- which.max(IPP[, 2])
    LIP <- length(IPP[, 2])
    ##
    IP_R13C <- 0
    if (LIP > x100) {
      M13C <- abs(IPP[, 1] - IPP[x100, 1] - 1.003354835336/molVecCharge[i])
      M13C <- M13C[(x100 + 1):LIP]
      x101 <- which.min(M13C)[1]
      if (M13C[x101] <= 0.015) {
        x101 <- x101 + x100
        IP_R13C <- IPP[x101, 2]/IPP[x100, 2]*100
      }
    }
    c(IPP[x100, 1], IP_R13C, x100, LIP)
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  if (formulaDeconvolutionCheck) {
    ##
    Elements <- element_sorter()[["Elements"]]
    x_el_c <- which(Elements == "C")
    LElements <- length(Elements)
    ##
    IonPW_DC <- ionization_pathway_deconvoluter(IonPathways, Elements)
    for (i in 1:length(IonPW_DC)) {
      ionStrCharge <- gsub("[+]|-", "", IonPW_DC[[i]][[3]])
      if (ionStrCharge == "") {
        IonPW_DC[[i]][[3]] <- 1
      } else {
        IonPW_DC[[i]][[3]] <- tryCatch(as.numeric(ionStrCharge), warning = function(w) {1}, error = function(e) {1})
      }
    }
    ##
    molf_deconvoluter <- function(i) {
      FormulaVector <- formula_vector_generator(molecularFormulaDatabase[i], Elements, LElements, allowedRedundantElements = TRUE)
      if (FormulaVector[x_el_c] > 0) {
        do.call(rbind, lapply(IonPW_DC, function(pathway) {
          Ion_coeff <- pathway[[1]]
          Ion_adduct <- pathway[[2]]
          Ion_charge <- pathway[[3]]
          MoleFormVec <- Ion_coeff*FormulaVector + Ion_adduct
          xNeg <- which(MoleFormVec < 0)
          if (length(xNeg) == 0) {
            if (retentionTimeCheck) {
              c(Ion_charge, retentionTime[i], MoleFormVec)
            } else {
              c(Ion_charge, MoleFormVec)
            }
          }
        }))
      }
    }
    ##
    ############################################################################
    ##
    if (allowedVerbose) {IPA_logRecorder("Initiated deconvoluting molecular formulas!")}
    ##
    if (number_processing_threads == 1) {
      MoleFormVecMat <- do.call(rbind, lapply(1:length(molecularFormulaDatabase), function(i) {
        molf_deconvoluter(i)
      }))
      ##
      ##########################################################################
      ##
    } else {
      ##
      ##########################################################################
      ##
      if (osType == "Windows") {
        ##
        clust <- makeCluster(number_processing_threads)
        clusterExport(clust, setdiff(ls(), c("clust")), envir = environment())
        ##
        MoleFormVecMat <- do.call(rbind, parLapply(clust, 1:length(molecularFormulaDatabase), function(i) {
          molf_deconvoluter(i)
        }))
        ##
        stopCluster(clust)
        ##
        ########################################################################
        ##
      } else {
        ##
        MoleFormVecMat <- do.call(rbind, mclapply(1:length(molecularFormulaDatabase), function(i) {
          molf_deconvoluter(i)
        }, mc.cores = number_processing_threads))
        ##
        closeAllConnections()
        ##
        ########################################################################
        ##
      }
    }
    ##
    molecularFormulaDatabase <- NULL
    ##
    if (allowedVerbose) {IPA_logRecorder("Completed deconvoluting molecular formulas!")}
    ##
    if (is.null(MoleFormVecMat)) {
      stop(IPA_logRecorder("Molecular formulas are not consistent with the ionization pathways!"))
    }
    ##
    if (retentionTimeCheck) {
      MoleFormVecMat <- matrix(MoleFormVecMat, ncol = (LElements + 2))
    } else {
      MoleFormVecMat <- matrix(MoleFormVecMat, ncol = (LElements + 1))
    }
    ##
    nMoleFormVecMat <- dim(MoleFormVecMat)[1]
    if (nMoleFormVecMat > 1) {
      MoleFormVecMat <- unique(as.matrix(MoleFormVecMat)) # To remove redundant rows
      nMoleFormVecMat <- dim(MoleFormVecMat)[1]
    }
    ##
    molVecCharge <- MoleFormVecMat[, 1]
    indexCol2Remove <- -1
    ##
    if (retentionTimeCheck) {
      retentionTime <- MoleFormVecMat[, 2]
      indexCol2Remove <- c(-1, -2)
    }
    ##
    MoleFormVecMat <- MoleFormVecMat[, indexCol2Remove]
    ##
    if (nMoleFormVecMat == 1) {
      MoleFormVecMat <- matrix(MoleFormVecMat, nrow = 1)
    }
    ##
    if (allowedVerbose) {IPA_logRecorder(paste0("There are `", nMoleFormVecMat, "` molecular formula ions for isotopic profile calculations!"))}
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  molecularFormulaMatrixElementSorterList <- molecular_formula_elements_filter(MoleFormVecMat, Elements)
  MoleFormVecMat <- molecularFormulaMatrixElementSorterList[["molecularFormulaMatrix"]]
  ##
  Elements <- molecularFormulaMatrixElementSorterList[["elementSorterList"]][["Elements"]]
  massAbundanceList <- molecularFormulaMatrixElementSorterList[["elementSorterList"]][["massAbundanceList"]]
  molecularFormulaMatrixElementSorterList <- NULL
  ##
  x_el_c <- which(Elements == "C")
  x_el_b <- which(Elements == "B")
  elementBcheck <- if (length(x_el_b) == 0) {FALSE} else {TRUE}
  x_el_br <- which(Elements == "Br")
  elementBrCheck <- if (length(x_el_br) == 0) {FALSE} else {TRUE}
  x_el_cl <- which(Elements == "Cl")
  elementClcheck <- if (length(x_el_cl) == 0) {FALSE} else {TRUE}
  x_el_k <- which(Elements == "K")
  elementKcheck <- if (length(x_el_k) == 0) {FALSE} else {TRUE}
  x_el_s <- which(Elements == "S")
  elementScheck <- if (length(x_el_s) == 0) {FALSE} else {TRUE}
  x_el_se <- which(Elements == "Se")
  elementSeCheck <- if (length(x_el_se) == 0) {FALSE} else {TRUE}
  x_el_si <- which(Elements == "Si")
  elementSiCheck <- if (length(x_el_si) == 0) {FALSE} else {TRUE}
  ##
  ##############################################################################
  ##############################################################################
  ##
  if ((number_processing_threads == 1) | allowedMustRunCalculation) {
    ##
    if (allowedVerbose) {IPA_logRecorder("Initiated calculating isotopic profiles!")}
    if (allowedVerbose) {progressBARboundaries <- txtProgressBar(min = 0, max = nMoleFormVecMat, initial = 0, style = 3)}
    ##
    IsotopicProfileList <- vector(mode = "list", nMoleFormVecMat)
    for (i in 1:nMoleFormVecMat) {
      IsotopicProfileList[[i]] <- tryCatch(IP_calculator(i),  error = function(e) {failed_IP_calculator(i)}, warning = function(w) {failed_IP_calculator(i)})
      ##
      if (allowedVerbose) {setTxtProgressBar(progressBARboundaries, i)}
    }
    ##
    if (allowedVerbose) {close(progressBARboundaries)}
    if (allowedVerbose) {IPA_logRecorder("Completed calculating isotopic profiles!")}
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    if (allowedVerbose) {IPA_logRecorder("Initiated calculating the database parameters!")}
    parametersIPDB <- do.call(rbind, lapply(1:nMoleFormVecMat, function(i) {
      parametersIPDBcalculator(i)
    }))
    ##
    ############################################################################
    ##
  } else {
    ##
    ############################################################################
    ##
    if (osType == "Windows") {
      ##
      ##########################################################################
      ##
      if (!allowedMustRunCalculation) {
        if (allowedVerbose) {IPA_logRecorder("Initiated calculating isotopic profiles!")}
        clust <- makeCluster(number_processing_threads)
        clusterExport(clust, setdiff(ls(), c("clust", "nMoleFormVecMat")), envir = environment())
        IsotopicProfileList <- tryCatch(parLapply(clust, 1:nMoleFormVecMat, function(i) {
          IP_calculator(i)
        }), warning = function(w) {stop("Isotopic profile calculations ran out of memory! Update parameter `FS0009` or apply a brute-force method `FS0005`!")})
        stopCluster(clust)
        if (allowedVerbose) {IPA_logRecorder("Completed calculating isotopic profiles!")}
      }
      ##
      if (allowedVerbose) {IPA_logRecorder("Initiated calculating the database parameters!")}
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, c("parametersIPDBcalculator", "IsotopicProfileList", "molVecCharge"), envir = environment())
      parametersIPDB <- do.call(rbind, parLapply(clust, 1:nMoleFormVecMat, function(i) {
        parametersIPDBcalculator(i)
      }))
      stopCluster(clust)
      ##
      ##########################################################################
      ##
    } else {
      ##
      if (!allowedMustRunCalculation) {
        if (allowedVerbose) {IPA_logRecorder("Initiated calculating isotopic profiles!")}
        IsotopicProfileList <- tryCatch(mclapply(1:nMoleFormVecMat, function(i) {
          IP_calculator(i)
        }, mc.cores = number_processing_threads), warning = function(w) {stop("Isotopic profile calculations ran out of memory! Update parameter `FS0009` or apply a brute-force method `FS0005`!")})
        if (allowedVerbose) {IPA_logRecorder("Completed calculating isotopic profiles!")}
      }
      ##
      if (allowedVerbose) {IPA_logRecorder("Initiated calculating the database parameters!")}
      parametersIPDB <- do.call(rbind, mclapply(1:nMoleFormVecMat, function(i) {
        parametersIPDBcalculator(i)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
      ##########################################################################
      ##
    }
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  if (allowedVerbose) {IPA_logRecorder("Completed calculating the database parameters!")}
  ##
  if (nMoleFormVecMat == 1) {
    parametersIPDB <- matrix(parametersIPDB, ncol = 4)
  }
  ##
  xNonInf <- which(!is.infinite(parametersIPDB[, 1]))
  LxNonInf <- length(xNonInf)
  if (LxNonInf > 0) {
    MassMAIso <- parametersIPDB[xNonInf, 1]
    IP_R13C <- parametersIPDB[xNonInf, 2]
    IndexMAIso <- parametersIPDB[xNonInf, 3]
    IPsize <- parametersIPDB[xNonInf, 4]
    parametersIPDB <- NULL
    ##
    MoleFormVecMat <- MoleFormVecMat[xNonInf, ]
    if (LxNonInf == 1) {
      MoleFormVecMat <- matrix(MoleFormVecMat, nrow = 1)
    }
    ##
    if (retentionTimeCheck) {
      retentionTime <- retentionTime[xNonInf]
    }
    ##
    if (allowedVerbose) {IPA_logRecorder("Initiated generating molecular formulas!")}
    MolecularFormulaVec <- hill_molecular_formula_printer(Elements, MoleFormVecMat, number_processing_threads)
    MoleFormVecMat <- NULL
    IsotopicProfileList <- IsotopicProfileList[xNonInf]
    if (allowedVerbose) {IPA_logRecorder("Completed generating molecular formulas!")}
  } else {
    MassMAIso <- Inf
    IP_R13C <- 0
    IndexMAIso <- 1
    IPsize <- 1
    MolecularFormulaVec <- ""
    IsotopicProfileList <- IsoProfInf100
    retentionTime <- -Inf
  }
  ##
  ##############################################################################
  ##
  AggregatedList <- aggregatedIPdbListGenerator(MassMAIso)
  ##
  ##############################################################################
  ##
  if (retentionTimeCheck) {
    IPDB <- list(AggregatedList, MassMAIso, MolecularFormulaVec, IsotopicProfileList, IP_R13C, IndexMAIso, IPsize, retentionTime)
    names(IPDB) <- c("AggregatedList", "MassMAIso", "MolecularFormula", "IsotopicProfile", "R13C", "IndexMAIso", "IPsize", "Retention Time")
  } else {
    IPDB <- list(AggregatedList, MassMAIso, MolecularFormulaVec, IsotopicProfileList, IP_R13C, IndexMAIso, IPsize)
    names(IPDB) <- c("AggregatedList", "MassMAIso", "MolecularFormula", "IsotopicProfile", "R13C", "IndexMAIso", "IPsize")
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  return(IPDB)
}
