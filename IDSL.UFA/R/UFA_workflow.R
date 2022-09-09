UFA_workflow <- function(spreadsheet) {
  ##
  gc()
  closeAllConnections()
  initiation_time <- Sys.time()
  print("Initiated testing the spreadsheet consistency!")
  PARAM <- UFA_xlsxAnalyzer(spreadsheet)
  if (length(PARAM) > 1) {
    print("The spreadsheet is consistent with the IDSL.UFA workflow!")
    ############################################################################
    x0001 <- PARAM[which(PARAM[, 1] == 'PARAM0001'), 2]
    if (tolower(x0001) == "yes") {
      x0002 <- PARAM[which(PARAM[, 1] == 'PARAM0002'), 2]
      if (tolower(x0002) == "yes") {
        PARAM_MF <- readxl::read_xlsx(spreadsheet, sheet = "enumerated_chemical_space")
        UFA_enumerated_chemical_space(PARAM_MF)
      }
      x0003 <- PARAM[which(PARAM[, 1] == 'PARAM0003'), 2]
      if (tolower(x0003) == "yes") {
        PARAM_SF <- readxl::read_xlsx(spreadsheet, sheet = "formula_source")
        PARAM_SF <- cbind(PARAM_SF[, 2], PARAM_SF[, 4])
        molecular_formulas_source_IPDB(PARAM_SF)
      }
    }
    ##
    x0005 <- PARAM[which(PARAM[, 1] == 'PARAM0005'), 2]
    if (tolower(x0005) == "yes") {
      ##########################################################################
      address_IPDB <- PARAM[which(PARAM[, 1] == "PARAM0004"), 2]
      print("Loading the isotopic profiles database!")
      IPDB <- loadRdata(address_IPDB)
      ## Temp
      if (length(IPDB) != 8) {
        stop("The selected IPDB is not consistent with the IDSL.UFA version 1.5!")
      }
      ##
      NPT <- as.numeric(PARAM[which(PARAM[, 1] == "PARAM0009"), 2])
      ##
      input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
      if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0011'), 2]) == "all") {
        file_name_hrms <- dir(path = input_path_hrms)
        file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
      } else {
        samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0011'), 2]
        file_name_hrms <- strsplit(samples_string, ";")[[1]]
      }
      ##
      input_path_pl <- PARAM[which(PARAM[, 1] == 'PARAM0013'), 2]
      file_names_peaklist1 <- dir(path = input_path_pl, pattern = ".Rdata")
      file_names_peaklist2 <- dir(path = input_path_pl, pattern = "peaklist_")
      file_names_peaklist <- file_names_peaklist1[file_names_peaklist1 %in% file_names_peaklist2]
      #
      file_names_peaklist_hrms1 <- gsub(".Rdata", "", file_names_peaklist)
      file_names_peaklist_hrms2 <- gsub("peaklist_", "", file_names_peaklist_hrms1)
      file_names_peaklist_hrms <- file_name_hrms %in% file_names_peaklist_hrms2
      L_PL <- length(which(file_names_peaklist_hrms == TRUE))
      if (length(file_name_hrms) != L_PL) {
        stop("Error!!! peaklist files are not available for the entire selected HRMS files!")
      }
      ##
      output_path <- PARAM[which(PARAM[, 1] == 'PARAM0014'), 2]
      output_path_annotated_mf_tables <- paste0(output_path, "/annotated_mf_tables")
      if (!dir.exists(output_path_annotated_mf_tables)) {
        dir.create(output_path_annotated_mf_tables, recursive = TRUE)
      }
      opendir(output_path_annotated_mf_tables)
      ##
      para_mode <- PARAM[which(PARAM[, 1] == 'PARAM0015'), 2]
      mass_accuracy <- as.numeric(PARAM[which(PARAM[, 1] == "PARAM0016"), 2])
      maxNEME <- as.numeric(PARAM[which(PARAM[, 1] == "PARAM0017"), 2])
      minPCS <- as.numeric(PARAM[which(PARAM[, 1] == "PARAM0018"), 2])
      minNDCS <- as.numeric(PARAM[which(PARAM[, 1] == "PARAM0019"), 2])
      minRCS <- as.numeric(PARAM[which(PARAM[, 1] == "PARAM0020"), 2])
      Score_coeff <- tryCatch(eval(parse(text = PARAM[which(PARAM[, 1] == 'PARAM0021'), 2])), error = function(e){rep(1, 5)})
      ##
      MF_library_search_TRUE <- FALSE
      x0022 <- PARAM[which(PARAM[, 1] == 'PARAM0022'), 2]
      if (tolower(x0022) == "yes") {
        MF_library_search_TRUE <- TRUE
        IonPathways <- tryCatch(eval(parse(text = PARAM[which(PARAM[, 1] == 'PARAM0023'), 2])), error = function(e){c("[M]")})
        ##
        x0024 <- PARAM[which(PARAM[, 1] == 'PARAM0024'), 2]
        MF_library <- loadRdata(x0024)
      }
      ##
      print("Initiated molecular formula annotation on individual peaklists!")
      call_molecular_formula_annotation <- function(i) {
        peaklist <- loadRdata(paste0(input_path_pl, "/peaklist_", file_name_hrms[i], ".Rdata"))
        ##
        outputer <- IDSL.IPA::IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
        spectraList <- outputer[["spectraList"]]
        ##
        MolecularFormulaAnnotationTable <- molecular_formula_annotator(IPDB, spectraList, peaklist, mass_accuracy, maxNEME, minPCS, minNDCS, minRCS, Score_coeff, number_processing_threads = NPT)
        ##
        if (MF_library_search_TRUE) {
          MolecularFormulaAnnotationTable <- molecular_formula_library_search(MolecularFormulaAnnotationTable, IPDB, MF_library, IonPathways, number_processing_threads = NPT)
        }
        ##
        save(MolecularFormulaAnnotationTable, file = paste0(output_path_annotated_mf_tables, "/MolecularFormulaAnnotationTable_", file_name_hrms[i], ".Rdata"))
        write.csv(MolecularFormulaAnnotationTable, file = paste0(output_path_annotated_mf_tables, "/MolecularFormulaAnnotationTable_", file_name_hrms[i], ".csv"))
        ##
        return()
      }
      if (para_mode == "peakmode" | NPT == 1) {
        progressBARboundaries <- txtProgressBar(min = 0, max = L_PL, initial = 0, style = 3)
        for (i in 1:L_PL) {
          setTxtProgressBar(progressBARboundaries, i)
          ##
          call_molecular_formula_annotation(i)
        }
        close(progressBARboundaries)
      } else if (para_mode == "samplemode") {
        NPT0 <- NPT
        NPT <- 1
        ##
        osType <- Sys.info()[['sysname']]
        ##
        if (osType == "Linux") {
          ##
          null_variable <- mclapply(1:L_PL, function (k) {
            call_molecular_formula_annotation(k)
          }, mc.cores = NPT0)
          ##
          closeAllConnections()
          ##
        } else if (osType == "Windows") {
          ##
          clust <- makeCluster(NPT0)
          registerDoParallel(clust)
          ##
          null_variable <- foreach(k = 1:L_PL, .verbose = FALSE) %dopar% {
            call_molecular_formula_annotation(k)
          }
          ##
          stopCluster(clust)
          ##
        }
      }
      print("Sucessfully completed molecular formula annotation on individual peaklists!")
      ##
      gc()
      closeAllConnections()
    }
    ############################################################################
    x0006 <- PARAM[which(PARAM[, 1] == 'PARAM0006'), 2]
    if (tolower(x0006) == "yes") {
      aligned_molecular_formula_annotator(PARAM)
    }
    ############################################################################
    x0007 <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
    if (tolower(x0007) == "yes") {
      UFA_score_coefficient_workflow(spreadsheet)
    }
    ############################################################################
    x0008 <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
    if (tolower(x0008) == "yes") {
      PARAM_SA <- UFA_profile_visualizer_xlsxAnalyzer(spreadsheet)
      AnnotatedSpectraTable <- UFA_profile_visualizer(PARAM_SA)
      ##
      exportedAnnotatedSpectraTable <- ifelse((tolower(PARAM_SA[which(PARAM_SA[, 1] == 'SA0011'), 2]) == "yes"), TRUE, FALSE)
      if (exportedAnnotatedSpectraTable) {
        output_path <- PARAM_SA[which(PARAM_SA[, 1] == 'PARAM0014'), 2]
        save(AnnotatedSpectraTable, file = paste0(output_path, "/AnnotatedSpectraTable.Rdata"))
        write.csv(AnnotatedSpectraTable, file = paste0(output_path, "/AnnotatedSpectraTable.csv"))
      }
    }
    ############################################################################
    required_time <- Sys.time() - initiation_time
    print(required_time)
  }
  return()
}
