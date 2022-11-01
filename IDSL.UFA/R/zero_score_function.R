zero_score_function <- function(PARAM_SFT) {
  ##
  address_IPDB <- PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0002"), 2]
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
  NPT <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0005"), 2])
  ##
  input_path_hrms <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0006'), 2]
  if (tolower(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0007'), 2]) == "all") {
    file_name_hrms <- dir(path = input_path_hrms)
    file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0008'), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
  } else {
    samples_string <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0007'), 2]
    file_name_hrms <- strsplit(samples_string, ";")[[1]]
  }
  ##
  input_path_pl <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0009'), 2]
  file_names_peaklist1 <- dir(path = input_path_pl, pattern = ".Rdata")
  file_names_peaklist2 <- dir(path = input_path_pl, pattern = "peaklist_")
  file_names_peaklist <- file_names_peaklist1[file_names_peaklist1 %in% file_names_peaklist2]
  L_PL <- length(file_names_peaklist)
  #
  file_names_peaklist_hrms1 <- gsub(".Rdata", "", file_names_peaklist)
  file_names_peaklist_hrms2 <- gsub("peaklist_", "", file_names_peaklist_hrms1)
  file_names_peaklist_hrms <- file_name_hrms %in% file_names_peaklist_hrms2
  if (length(which(file_names_peaklist_hrms == TRUE)) != L_PL) {
    stop(IPA_logRecorder("Error!!! peaklist files are not available for the entire selected HRMS files!"))
  }
  ##
  output_path <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0010'), 2]
  output_path_score_function_calculations <- paste0(output_path, "/score_function_calculations")
  if (!dir.exists(output_path_score_function_calculations)) {
    dir.create(output_path_score_function_calculations, recursive = TRUE)
  }
  opendir(output_path_score_function_calculations)
  ##
  IPA_logRecorder("Deconvoluting the reference spreadsheet file!")
  excelfile_address <- PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0011"), 2]
  excelfile <- readxl::read_xlsx(excelfile_address)
  FileName <- excelfile$FileName
  Molf <- excelfile$MolcularFormula
  Ionf <- excelfile$IonizationPathway
  RTf <- as.numeric(excelfile$`RetentionTime(min)`)
  ##
  Molf <- gsub(" ", "", Molf, fixed = TRUE)
  Molf <- gsub("[+]", "", Molf, fixed = TRUE)
  Molf <- gsub("-", "", Molf, fixed = TRUE)
  ##############################################################################
  EL <- element_sorter()
  EL_alpha <- EL[[1]]
  Mass_Abundance <- EL[[2]]
  L_EL <- length(EL_alpha)
  L_Molf <- length(Molf)
  ##
  u_Ionf <- unique(Ionf)
  index_ion_dc <- rep(0, L_Molf)
  for (i in 1:length(u_Ionf)) {
    x_grep <- which(Ionf == u_Ionf[i])
    index_ion_dc[x_grep] <- i
  }
  u_ion_dc <- ionization_pathway_deconvoluter(u_Ionf, EL_alpha)
  ##
  MolVecMatList_call <- function(k) {
    mol_v <- formula_vector_generator(Molf[k], EL_alpha, L_EL)
    #
    ion_dc <- u_ion_dc[[index_ion_dc[k]]]
    ion_coeff <- ion_dc[[1]]
    ion_adduct <- ion_dc[[2]]
    #
    mol_v*ion_coeff + ion_adduct
  }
  ##
  mass_accuracy <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0012"), 2])
  rt_error <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0013"), 2])
  maxNEME <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0014"), 2])
  minPCS <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0015"), 2])
  minNDCS <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0016"), 2])
  minRCS <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0017"), 2])
  Score_coeff <- rep(1, 5)
  ##
  Entire_final_list_unoptimized_call <- function(i) {
    e_u <- NULL
    ##
    x_ss <- which(FileName == file_name_hrms[i])
    L_targeted <- length(x_ss)
    if (L_targeted > 0) {
      mz_ss <- mzf[x_ss]
      rt_ss <- RTf[x_ss]
      ##
      peaklist <- IDSL.IPA::loadRdata(paste0(input_path_pl, "/peaklist_", file_name_hrms[i], ".Rdata"))
      x_ss_i_j <- do.call(rbind, lapply(1:L_targeted, function(j) {
        x <-mzRTindexer(peaklist[, 8], peaklist[, 3], mz_ss[j], rt_ss[j], mass_accuracy, rt_error)
        ##
        if (!is.null(x)) {
          c(x, j)
        }
      }))
      x_ss_i <- x_ss_i_j[, 1]
      L_x_ss_i <- length(x_ss_i)
      if (L_x_ss_i > 0) {
        peaklist <- matrix(peaklist[x_ss_i, ], nrow = L_x_ss_i)
        ##
        outputer <- IDSL.IPA::IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
        spectraList <- outputer[["spectraList"]]
        ##
        FinalList <- molecular_formula_annotator(IPDB, spectraList, peaklist, mass_accuracy, maxNEME, minPCS, minNDCS, minRCS, Score_coeff, number_processing_threads = 1)
        ##
        if (length(FinalList) > 0) {
          Molf_ss <- Molf_product[x_ss[x_ss_i_j[, 2]]]
          x_fl_compounds <- do.call(rbind, lapply(1:L_x_ss_i, function(id) {
            id_compound <- NULL
            x_id <- which(FinalList[, 1] == id)
            x_compound <- which(FinalList[x_id, 3] == Molf_ss[id])
            if (length(x_compound) > 0) {
              id_compound <- x_id[1] + x_compound - 1
            }
            id_compound
          }))
          true_false <- rep(0, dim(FinalList)[1])
          true_false[x_fl_compounds] <- 1
          ##
          p_id <- as.numeric(unique(FinalList[, 1]))
          p_count <- do.call(rbind, lapply(1:length(p_id), function(k) {
            x_id <- length(which(FinalList[, 1] == p_id[k]))
            matrix(rep(x_id, x_id), ncol = 1)
          }))
          ##
          e_u <- cbind(rep(file_name_hrms[i], dim(FinalList)[1]), FinalList, p_count, true_false)
        }
      }
    }
    e_u
  }
  ##
  osType <- Sys.info()[['sysname']]
  ##
  if (osType == "Windows") {
    clust <- makeCluster(NPT)
    registerDoParallel(clust)
    ##
    MolVecMatList <- foreach(k = 1:L_Molf, .combine = 'rbind', .verbose = FALSE) %dopar% {
      MolVecMatList_call(k)
    }
    ##
    stopCluster(clust)
    ##
  } else if (osType == "Linux") {
    ##
    MolVecMatList <- do.call(rbind, mclapply(1:L_Molf, function(k) {
      MolVecMatList_call(k)
    }, mc.cores = NPT))
    ##
    closeAllConnections()
  }
  ##
  Molf_product <- hill_molecular_formula_printer(EL_alpha, MolVecMatList, NPT)
  MolVecMatList <- 0
  ##
  Molf_IPDB <- hill_molecular_formula_printer(IPDB[["MolecularFormulaDB"]][["Elements"]], IPDB[["MolecularFormulaDB"]][["MolecularFormulaMatrix"]], NPT)
  x_IPDB <- which((Molf_IPDB %in% Molf_product) == TRUE)
  mzf <- rep(NA, L_Molf)
  for (k in x_IPDB) {
    x_product <- which(Molf_product == Molf_IPDB[k])
    if (length(x_product) > 0) {
      mzf[x_product] <- IPDB[["MassMAIso"]][k]
    }
  }
  Molf_IPDB <- NULL
  ##
  x_NA <- which(is.na(mzf))
  if (length(x_NA) > 0) {
    NA_molf <- unique(Molf[x_NA])
    IPA_logRecorder("WARNING!!! The following molecular formulas were not included in the isotopic profile database (IPDB):")
    for (i in 1:length(NA_molf)) {
      IPA_logRecorder(NA_molf[i])
    }
  }
  ##
  IPA_logRecorder("Initiated producing the unoptimized list of candidate molecular formulas!")
  if (osType == "Windows") {
    clust <- makeCluster(NPT)
    registerDoParallel(clust)
    ##
    Entire_final_list_unoptimized <- foreach(i = 1:L_PL, .combine = 'rbind', .verbose = FALSE) %dopar% {
      Entire_final_list_unoptimized_call(i)
    }
    ##
    stopCluster(clust)
    ##
  } else if (osType == "Linux") {
    Entire_final_list_unoptimized <- do.call(rbind, mclapply(1:L_PL, function (i) {
      Entire_final_list_unoptimized_call(i)
    }, mc.cores = NPT))
    ##
    closeAllConnections()
  }
  ##############################################################################
  counter_c <- 1
  L_Entire_final_list_unoptimized <- dim(Entire_final_list_unoptimized)[1]
  if (L_Entire_final_list_unoptimized > 0) {
    ##
    CompoundID <- rep(1, L_Entire_final_list_unoptimized)
    if (L_Entire_final_list_unoptimized > 1) {
      for (i in 2:L_Entire_final_list_unoptimized) {
        if ((Entire_final_list_unoptimized[(i - 1), 1] != Entire_final_list_unoptimized[i, 1]) |
            (Entire_final_list_unoptimized[(i - 1), 2] != Entire_final_list_unoptimized[i, 2])) {
          counter_c <- counter_c + 1
        }
        CompoundID[i] <- counter_c
      }
    }
    ##
    SizeIP_IsotopicProfile_DataBase <- IPDB[["IPsize"]]
    x_ip <- SizeIP_IsotopicProfile_DataBase[as.numeric(Entire_final_list_unoptimized[, 3])]
    ##
    Entire_final_list_unoptimized <- cbind(Entire_final_list_unoptimized[, 1:8], x_ip, Entire_final_list_unoptimized[, 9:16], CompoundID, Entire_final_list_unoptimized[, dim(Entire_final_list_unoptimized)[2]])
    Entire_final_list_unoptimized <- data.frame(Entire_final_list_unoptimized)
    colnames(Entire_final_list_unoptimized) <- c("FileName", "PeakID", "ID_IonFormula",
                                                 "IonFormula", "m/z Isotopic Profile", "m/z peaklist",
                                                 "RT(min)", "PeakHeight", "size IP", "NEME(mDa)", "PCS",
                                                 "R13C peakList", "R13C Isotopic Profile", "NDCS", "RCS(%)",
                                                 "Rank", "CandidateCount", "CompoundID", "MolFMatch")
    rownames(Entire_final_list_unoptimized) <- NULL
    save(Entire_final_list_unoptimized, file = paste0(output_path_score_function_calculations, "/Entire_final_list_unoptimized.Rdata"))
    IPA_logRecorder("Completed producing the unoptimized list of candidate molecular formulas!")
    ##
  } else {
    stop(IPA_logRecorder("Production of the unoptimized list of candidate molecular formulas was not successful!!!"))
  }
  ##
  gc()
  closeAllConnections()
  ##
}
