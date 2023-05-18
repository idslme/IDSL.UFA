scoreCoefficientsReplicate <- function(PARAM_ScoreFunc) {
  ##
  address_IPDB <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0002'), 2]
  if (!is.na(address_IPDB)) {
    if (!file.exists(address_IPDB)) {
      stop(IPA_logRecorder("ERROR!!! Problem with SFT0002! The isotopic profile database file is not available!"))
    }
  } else {
    stop(IPA_logRecorder("ERROR!!! Problem with SFT0002! The isotopic profile database file is not available!"))
  }
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder(paste0("Loading the Isotopic Profiles DataBase (IPDB) from `", address_IPDB, "`!"))
  IPDB <- IDSL.IPA::loadRdata(address_IPDB)
  ##
  input_path_hrms <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0005'), 2]
  ##
  IPA_logRecorder("Deconvoluting the reference spreadsheet file!")
  excelfile_address <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == "SFT0006"), 2]
  excelfile <- readxl::read_xlsx(excelfile_address, col_names = TRUE)
  ##
  listIndexHRMSfileNames <- base::tapply(seq(1, length(excelfile$'FileName'), 1), excelfile$'FileName', 'c', simplify = FALSE)
  HRMSfileNames <- names(listIndexHRMSfileNames)
  LHRMS <- length(HRMSfileNames)
  ##
  inputPathPeaklist <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0007'), 2]
  ##
  peaklistFileNames <- dir(path = inputPathPeaklist, pattern = ".Rdata$", ignore.case = TRUE)
  peaklistFileNames <- peaklistFileNames[grep("^peaklist_", peaklistFileNames, ignore.case = TRUE)]
  L_PL <- length(peaklistFileNames)
  ##
  if (LHRMS > L_PL) {
    peaklistHRMSfileNames <- paste0("peaklist_", HRMSfileNames, ".Rdata")
    ndPeaklists <- setdiff(peaklistHRMSfileNames, peaklistFileNames)
    ndPeaklists <- gsub("^peaklist_|.Rdata$", "", ndPeaklists, ignore.case = TRUE)
    IPA_logRecorder("Error!!! peaklist files are not available for the following HRMS file(s):")
    for (i in ndPeaklists) {
      IPA_logRecorder(i)
    }
    stop()
  }
  ##
  output_path <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0008'), 2]
  output_path_score_function_calculations <- paste0(output_path, "/score_function_calculations")
  if (!dir.exists(output_path_score_function_calculations)) {
    dir.create(output_path_score_function_calculations, recursive = TRUE)
  }
  ##
  RTf <- as.numeric(excelfile$'RetentionTime(min)')
  ##
  Molf <- excelfile$'MolcularFormula'
  Molf <- gsub(" ", "", Molf, fixed = TRUE)
  Molf <- gsub("[+]", "", Molf, fixed = TRUE)
  Molf <- gsub("-", "", Molf, fixed = TRUE)
  LMolf <- length(Molf)
  ##
  Elements <- element_sorter()[["Elements"]]
  LElements <- length(Elements)
  ##
  Ionf <- excelfile$'IonizationPathway'
  Ionf <- gsub(" ", "", Ionf, fixed = TRUE)
  Ionf <- gsub("*", "", Ionf, fixed = TRUE)
  uniqueIonf <- unique(Ionf)
  uniqueIonPathways <- ionization_pathway_deconvoluter(uniqueIonf, Elements)
  ##
  ##############################################################################
  ##
  call_MoleFormVecMat <- function(i) {
    mol_v <- formula_vector_generator(Molf[i], Elements, LElements, allowedRedundantElements = TRUE)
    #
    ion_dc <- uniqueIonPathways[[Ionf[i]]]
    ion_coeff <- ion_dc[[1]]
    ion_adduct <- ion_dc[[2]]
    #
    molvec <- mol_v*ion_coeff + ion_adduct
    xNeg <- which(molvec < 0)
    if (length(xNeg) == 0) {
      molvec
    }
  }
  ##
  NPT <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0009'), 2])
  RTtolerance <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0010'), 2])
  massAccuracy <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0011'), 2])
  maxNEME <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0012'), 2])
  minPCS <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0013'), 2])
  minNDCS <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0014'), 2])
  minRCS <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0015'), 2])
  scoreCoefficients <- rep(1, 5)
  ##
  ##############################################################################
  ##
  call_scoreCoefficientsReplicate <- function(i) {
    ##
    x_ss <- listIndexHRMSfileNames[[HRMSfileNames[i]]]
    ##
    Ltargeted <- length(x_ss)
    if (Ltargeted > 0) {
      mz_ss <- mzf[x_ss]
      rt_ss <- RTf[x_ss]
      ##
      peaklist <- IDSL.IPA::loadRdata(paste0(inputPathPeaklist, "/peaklist_", HRMSfileNames[i], ".Rdata"))
      x_ss_i_j <- do.call(rbind, lapply(1:Ltargeted, function(j) {
        x <- mzRTindexer(peaklist[, 8], peaklist[, 3], mz_ss[j], rt_ss[j], massAccuracy, RTtolerance)
        ##
        if (!is.null(x)) {
          c(x, j)
        }
      }))
      L_x_ss_i <- length(x_ss_i_j)/2
      if (L_x_ss_i > 0) {
        x_ss_i <- x_ss_i_j[, 1]
        ##
        outputer <- IDSL.IPA::IPA_MSdeconvoluter(input_path_hrms, HRMSfileNames[i])
        spectraList <- outputer[["spectraList"]]
        outputer <- NULL
        ##
        FinalList <- molecular_formula_annotator(IPDB, spectraList, peaklist, selectedIPApeaks = x_ss_i, massAccuracy, maxNEME, minPCS,
                                                 minNDCS, minRCS, scoreCoefficients, RTtolerance = NA, correctedRTpeaklist = NULL,
                                                 exportSpectraParameters = NULL, number_processing_threads = 1)
        ##
        if (!is.null(FinalList)) {
          FinalList[, 1] <- as.numeric(FinalList[, 1])
          Molf_ss <- Molf_product[x_ss[x_ss_i_j[, 2]]]
          x_fl_compounds <- do.call(c, lapply(1:L_x_ss_i, function(id) {
            x_id <- which(FinalList[, 1] == x_ss_i[id])
            if (length(x_id) > 0) {
              x_compound <- which(FinalList[x_id, 4] == Molf_ss[id])
              if (length(x_compound) > 0) {
                x_id[x_compound]
              }
            }
          }))
          LFinalList <- dim(FinalList)[1]
          true_false <- rep(0, LFinalList)
          true_false[x_fl_compounds] <- 1
          ##
          xDiff <- c(0, which(abs(diff(FinalList[, 1])) > 0), LFinalList)
          p_count <- do.call(c, lapply(1:(length(xDiff) - 1), function(id) {
            Lid <- xDiff[id + 1] - xDiff[id]
            rep(Lid, Lid)
          }))
          ##
          cbind(rep(HRMSfileNames[i], LFinalList), FinalList, p_count, true_false)
        }
      }
    }
  }
  ##
  osType <- Sys.info()[['sysname']]
  ##
  if (osType == "Windows") {
    clust <- makeCluster(NPT)
    clusterExport(clust, setdiff(ls(), c("clust", "LMolf")), envir = environment())
    ##
    MoleFormVecMat <- do.call(rbind, parLapply(clust, 1:LMolf, function(i) {
      call_MoleFormVecMat(i)
    }))
    ##
    stopCluster(clust)
    ##
  } else {
    ##
    MoleFormVecMat <- do.call(rbind, mclapply(1:LMolf, function(i) {
      call_MoleFormVecMat(i)
    }, mc.cores = NPT))
    ##
    closeAllConnections()
  }
  ##
  ##############################################################################
  ##
  molecularFormulaMatrixElementSorterList <- molecular_formula_elements_filter(MoleFormVecMat, Elements)
  MoleFormVecMat <- molecularFormulaMatrixElementSorterList[["molecularFormulaMatrix"]]
  Elements <- molecularFormulaMatrixElementSorterList[["elementSorterList"]][["Elements"]]
  molecularFormulaMatrixElementSorterList <- NULL
  ##
  ##############################################################################
  ##
  Molf_product <- hill_molecular_formula_printer(Elements, MoleFormVecMat, NPT)
  MoleFormVecMat <- NULL
  ##
  xIPDB <- which(IPDB[["MolecularFormula"]] %in% Molf_product)
  mzf <- rep(NA, LMolf)
  for (i in xIPDB) {
    x_product <- which(Molf_product == IPDB[["MolecularFormula"]][i])
    if (length(x_product) > 0) {
      mzf[x_product] <- IPDB[["MassMAIso"]][i]
    }
  }
  ##
  xNA <- which(is.na(mzf))
  LxNA <- length(xNA)
  ##
  if (LxNA == LMolf) {
    stop("No `ion molecular formula` from reference spreadsheet was detected in the IPDB!")
  }
  ##
  if (LxNA > 0) {
    NA_molf <- unique(Molf[xNA])
    IPA_logRecorder("WARNING!!! The following `ion molecular formulas` were not included in the isotopic profile database (IPDB):")
    for (i in NA_molf) {
      IPA_logRecorder(i)
    }
  }
  ##
  IPA_logRecorder("Initiated producing the unoptimized list of candidate molecular formulas!")
  if (osType == "Windows") {
    ##
    clust <- makeCluster(NPT)
    clusterExport(clust, setdiff(ls(), c("clust", "LHRMS")), envir = environment())
    ##
    unoptimized_molecular_formula_annotation <- do.call(rbind, parLapply(clust, 1:LHRMS, function(i) {
      tryCatch(call_scoreCoefficientsReplicate(i),
               error = function(e) {IPA_logRecorder(paste0("Problem with `", HRMSfileNames[i],"`!"))})
    }))
    ##
    stopCluster(clust)
    ##
    ############################################################################
    ##
  } else {
    unoptimized_molecular_formula_annotation <- do.call(rbind, mclapply(1:LHRMS, function(i) {
      tryCatch(call_scoreCoefficientsReplicate(i),
               error = function(e) {IPA_logRecorder(paste0("Problem with `", HRMSfileNames[i],"`!"))})
    }, mc.cores = NPT))
    ##
    closeAllConnections()
    ##
  }
  ##############################################################################
  counterCompound <- 1
  dimunoptimized_molecular_formula_annotation <- dim(unoptimized_molecular_formula_annotation)
  L_unoptimized_molecular_formula_annotation <- dimunoptimized_molecular_formula_annotation[1]
  if (L_unoptimized_molecular_formula_annotation > 0) {
    ##
    CompoundID <- rep(1, L_unoptimized_molecular_formula_annotation)
    if (L_unoptimized_molecular_formula_annotation > 1) {
      for (i in 2:L_unoptimized_molecular_formula_annotation) {
        if ((unoptimized_molecular_formula_annotation[(i - 1), 1] != unoptimized_molecular_formula_annotation[i, 1]) |
            (unoptimized_molecular_formula_annotation[(i - 1), 2] != unoptimized_molecular_formula_annotation[i, 2])) {
          counterCompound <- counterCompound + 1
        }
        CompoundID[i] <- counterCompound
      }
    }
    ##
    noProperties <- dimunoptimized_molecular_formula_annotation[2]
    unoptimized_molecular_formula_annotation <- cbind(unoptimized_molecular_formula_annotation[, 1:(noProperties - 1)], CompoundID, unoptimized_molecular_formula_annotation[, noProperties])
    unoptimized_molecular_formula_annotation <- data.frame(unoptimized_molecular_formula_annotation)
    colnames(unoptimized_molecular_formula_annotation) <- c('FileName', 'IDSL.IPA_PeakID', 'IonFormulaID', 'isotopologueCount', 'IonFormula', 'mzTheoretical',
                                                            'mzPeaklist', 'retentionTime', 'PeakHeight', 'normEucMassError', 'PCS', 'r13cPeaklist',
                                                            'r13cTheoretical', 'NDCS', 'ratioChromScan', 'rank', 'CandidateCount', 'CompoundID', 'MolFMatch')
    rownames(unoptimized_molecular_formula_annotation) <- NULL
    save(unoptimized_molecular_formula_annotation, file = paste0(output_path_score_function_calculations, "/unoptimized_molecular_formula_annotation.Rdata"))
    IPA_logRecorder("Completed producing the unoptimized list of candidate molecular formulas!")
    ##
  } else {
    stop(IPA_logRecorder("Production of the unoptimized list of candidate molecular formulas was not successful!!!"))
  }
  ##
  gc()
  closeAllConnections()
  ##
  return()
}
