UFA_xlsxAnalyzer <- function(spreadsheet) {
  ##
  print("Initiated testing the spreadsheet consistency!")
  ##
  PARAM_ECS <- NULL
  PARAM_FormSource <- NULL
  PARAM_ScoreFunc <- NULL
  ##
  checkpoint_parameter <- FALSE
  if (length(spreadsheet) >= 4) {
    if (typeof(spreadsheet) == "list") {
      PARAM <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
    } else {
      print("The UFA spreadsheet was not produced properly!")
    }
  } else if (length(spreadsheet) == 1) {
    if (typeof(spreadsheet) == "character") {
      if (file.exists(spreadsheet)) {
        spreadsheet_UFA <- readxl::read_xlsx(spreadsheet, sheet = "parameters")
        PARAM <- cbind(spreadsheet_UFA[, 2], spreadsheet_UFA[, 4])
        checkpoint_parameter <- TRUE
      } else {
        print("The UFA spreadsheet not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      print("The UFA spreadsheet was not produced properly!")
    }
  } else {
    print("The UFA spreadsheet was not produced properly!")
  }
  if (checkpoint_parameter) {
    ############################################################################
    ########################### Global parameters ##############################
    ############################################################################
    x0001 <- which(PARAM[, 1] == 'PARAM0001')
    if (length(x0001) == 0) {
      print("ERROR!!! Problem with PARAM0001!")
      checkpoint_parameter <- FALSE
      PARAM0001 <- "no"
    } else {
      PARAM0001 <- tolower(PARAM[x0001, 2])
      if (!(PARAM0001 == "yes" | PARAM0001 == "no")) {
        print("ERROR!!! Problem with PARAM0001!")
        checkpoint_parameter <- FALSE
      } else {
        PARAM[x0001, 2] <- PARAM0001
      }
    }
    ##
    if (PARAM0001 == "yes") {
      x0002 <- which(PARAM[, 1] == 'PARAM0002')
      if (length(x0002) == 0) {
        print("ERROR!!! Problem with PARAM0002!")
        checkpoint_parameter <- FALSE
        PARAM0002 <- "no"
      } else {
        PARAM0002 <- tolower(PARAM[x0002, 2])
        if (!(PARAM0002 == "yes" | PARAM0002 == "no")) {
          print("ERROR!!! Problem with PARAM0002!")
          checkpoint_parameter <- FALSE
        } else {
          PARAM[x0002, 2] <- PARAM0002
        }
      }
      ##
      x0003 <- which(PARAM[, 1] == 'PARAM0003')
      if (length(x0003) == 0) {
        print("ERROR!!! Problem with PARAM0003!")
        checkpoint_parameter <- FALSE
        PARAM0003 <- "no"
      } else {
        PARAM0003 <- tolower(PARAM[x0003, 2])
        if (!(PARAM0003 == "yes" | PARAM0003 == "no")) {
          print("ERROR!!! Problem with PARAM0003!")
          checkpoint_parameter <- FALSE
        } else {
          PARAM[x0003, 2] <- PARAM0003
        }
      }
      ##
      if ((PARAM0002 == "yes" & PARAM0003 == "yes") | (PARAM0002 == "no" & PARAM0003 == "no")) {
        PARAM0002 <- "no"
        PARAM0003 <- "no"
        print("ERROR!!! Problem with PARAM0002 & PARAM0003!")
        print("ERROR!!! Both PARAM0002 & PARAM0003 cannot be 'YES' and 'NO' at the same time! You may choose only one method to generate isotopic profiles database (IPDB)!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0005 <- which(PARAM[, 1] == 'PARAM0005')
    if (length(x0005) == 0) {
      print("ERROR!!! Problem with PARAM0005!")
      checkpoint_parameter <- FALSE
      PARAM0005 <- "no"
    } else {
      PARAM0005 <- tolower(PARAM[x0005, 2])
      if (!(PARAM0005 == "yes" | PARAM0005 == "no")) {
        print("ERROR!!! Problem with PARAM0005!")
        checkpoint_parameter <- FALSE
      } else {
        PARAM[x0005, 2] <- PARAM0005
      }
    }
    ##
    x0004 <- which(PARAM[, 1] == 'PARAM0004') # This must be here
    ##
    if (PARAM0005 == "yes") {
      PARAM0004 <- PARAM[x0004, 2]
      PARAM0004 <- gsub("\\", "/", PARAM0004, fixed = TRUE)
      PARAM[x0004, 2] <- PARAM0004
      if (!is.na(PARAM0004)) {
        if (!file.exists(PARAM0004)) {
          print("ERROR!!! Problem with PARAM0004! The isotopic profile database (IPDB) file is not available!")
          checkpoint_parameter <- FALSE
        }
      } else {
        print("ERROR!!! Problem with PARAM0004!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0006 <- which(PARAM[, 1] == 'PARAM0006')
    if (length(x0006) == 0) {
      print("ERROR!!! Problem with PARAM0006!")
      checkpoint_parameter <- FALSE
      PARAM0006 <- "no"
    } else {
      PARAM0006 <- tolower(PARAM[x0006, 2])
      if (!(PARAM0006 == "yes" | PARAM0006 == "no")) {
        print("ERROR!!! Problem with PARAM0006!")
        checkpoint_parameter <- FALSE
      } else {
        PARAM[x0006, 2] <- PARAM0006
      }
    }
    ##
    x0007 <- which(PARAM[, 1] == 'PARAM0007')
    if (length(x0007) == 0) {
      print("ERROR!!! Problem with PARAM0007!")
      checkpoint_parameter <- FALSE
      PARAM0007 <- "no"
    } else {
      PARAM0007 <- tolower(PARAM[x0007, 2])
      if (!(PARAM0007 == "yes" | PARAM0007 == "no")) {
        print("ERROR!!! Problem with PARAM0007!")
        checkpoint_parameter <- FALSE
      } else {
        PARAM[x0007, 2] <- PARAM0007
      }
    }
    ##
    NPT <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2])
    if (length(NPT) == 0) {
      print("ERROR!!! Problem with PARAM0008! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (NPT >= 1) {
        if ((NPT %% 1) != 0) {
          print("ERROR!!! Problem with PARAM0008! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        print("ERROR!!! Problem with PARAM0008! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ############################################################################
    ######################## Data import and export ############################
    ############################################################################
    if (PARAM0005 == "yes") {
      x0009 <- which(PARAM[, 1] == 'PARAM0009')
      if (length(x0009) == 0) {
        print("ERROR!!! Problem with PARAM0009!")
        checkpoint_parameter <- FALSE
      } else {
        input_path_hrms <- PARAM[x0009, 2]
        input_path_hrms <- gsub("\\", "/", input_path_hrms, fixed = TRUE)
        PARAM[x0009, 2] <- input_path_hrms
        if (!dir.exists(input_path_hrms)) {
          print("ERROR!!! Problem with PARAM0009! Please make sure the full path is provided!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      LHRMS <- 0
      x0010 <- which(PARAM[, 1] == 'PARAM0010')
      if (length(x0010) == 0) {
        print("ERROR!!! Problem with PARAM0010!")
        checkpoint_parameter <- FALSE
      } else {
        if (tolower(PARAM[x0010, 2]) == "all") {
          file_name_hrms <- dir(path = input_path_hrms)
          file_name_hrms <- file_name_hrms[grep(pattern = ".mzML$|.mzXML$|.CDF$", file_name_hrms, ignore.case = TRUE)]
          LHRMS <- length(file_name_hrms)
          if (LHRMS == 0) {
            print("ERROR!!! Problem with PARAM0010! No mzML/mzXML/CDF file was detected in the folder!")
          }
        } else {
          samples_string <- PARAM[x0010, 2]
          file_name_hrms <- strsplit(samples_string, ";")[[1]]
          LHRMS <- length(file_name_hrms)
          ndHRMS <- do.call(c, lapply(file_name_hrms, function(i) {
            if (!file.exists(paste0(input_path_hrms, "/", i))) {
              i
            }
          }))
          ##
          if (!is.null(ndHRMS)) {
            print("ERROR!!! Problem with PARAM0010! not detected the following file(s) (case sensitive even for file extensions):")
            for (i in ndHRMS) {
              print(i)
            }
            checkpoint_parameter <- FALSE
          }
        }
      }
      ##
      x0011 <- which(PARAM[, 1] == 'PARAM0011')
      if (length(x0011) == 0) {
        print("ERROR!!! Problem with PARAM0011!")
        checkpoint_parameter <- FALSE
      } else {
        inputPathPeaklist <- PARAM[x0011, 2]
        inputPathPeaklist <- gsub("\\", "/", inputPathPeaklist, fixed = TRUE)
        PARAM[x0011, 2] <- inputPathPeaklist
        if (!dir.exists(inputPathPeaklist)) {
          print("ERROR!!! Problem with PARAM0011! Please make sure the full path is provided!")
          checkpoint_parameter <- FALSE
        } else {
          ######################################################################
          ## To see if the entire peaklists were generated for all HRMS files ##
          ######################################################################
          if (LHRMS > 0) {
            peaklistFileNames <- dir(path = inputPathPeaklist, pattern = ".Rdata$")
            peaklistFileNames <- peaklistFileNames[grep("^peaklist_", peaklistFileNames)]
            L_PL <- length(peaklistFileNames)
            ##
            if (LHRMS > L_PL) {
              checkpoint_parameter <- FALSE
              peaklistHRMSfileNames <- paste0("peaklist_", file_name_hrms, ".Rdata")
              ndPeaklists <- setdiff(peaklistHRMSfileNames, peaklistFileNames)
              ndPeaklists <- gsub("^peaklist_|.Rdata$", "", ndPeaklists)
              print("Error!!! peaklist files are not available for the following HRMS file(s):")
              for (i in ndPeaklists) {
                print(i)
              }
            }
          }
        }
      }
    }
    ##
    peak_alignment_folder_check <- TRUE
    if (PARAM0006 == "yes") {
      listAlignmentFolderCheck <- IPA_peak_alignment_folder_xlsxAnalyzer(PARAM, PARAM_ID = 'PARAM0012', checkpoint_parameter, correctedRTcheck = FALSE, CSAcheck = FALSE, allowedVerbose = TRUE)
      PARAM <- listAlignmentFolderCheck[[1]]
      checkpoint_parameter <- listAlignmentFolderCheck[[2]]
      peak_alignment_folder_check <- FALSE
    }
    ##
    if (PARAM0005 == "yes") {
      ##
      if (LHRMS == 1) {
        PARAM0013 <- PARAM[which(PARAM[, 1] == "PARAM0013"), 2]
        if (is.na(PARAM0013)) {
          checkpoint_parameter <- FALSE
          print("ERROR!!! Problem with PARAM0013! This parameter should be 'All' or a vector of indices!")
        } else if (gsub(" ", "", tolower(PARAM0013)) == "all") {
          print("The enitre 12C m/z values in the peaklist were placed in the processing row!")
        } else {
          peaklist <- IDSL.IPA::loadRdata(paste0(inputPathPeaklist, "/peaklist_", file_name_hrms, ".Rdata"))
          n_peaks <- dim(peaklist)[1]
          ##
          selectedIPApeaks <- tryCatch(eval(parse(text = paste0("c(", PARAM0013, ")"))), error = function(e){NULL})
          if (is.null(selectedIPApeaks) | (max(selectedIPApeaks) > n_peaks)) {
            checkpoint_parameter <- FALSE
            print("ERROR!!! Problem with PARAM0013! The range of indices are out of the peaklist dimension!")
          } else {
            print("The following peak IDs were selected for processing: ")
            for (id in 1:length(selectedIPApeaks)) {
              print(paste0(selectedIPApeaks[id], " - ", peaklist[selectedIPApeaks[id], 3],  " - ", peaklist[selectedIPApeaks[id], 8]))
            }
          }
        }
      }
    }
    if (PARAM0005 == "yes" | PARAM0006 == "yes") {
      ##
      x0014 <- which(PARAM[, 1] == 'PARAM0014')
      if (length(x0014) == 0) {
        print("ERROR!!! Problem with PARAM0014!")
        checkpoint_parameter <- FALSE
      } else {
        output_path <- gsub("\\", "/", PARAM[x0014, 2], fixed = TRUE)
        PARAM[x0014, 2] <- output_path
        if (!dir.exists(output_path)) {
          tryCatch(dir.create(output_path, recursive = TRUE), warning = function(w){warning("Problem with PARAM0014! R cannot create the folder!")})
          if (!dir.exists(output_path)) {
            checkpoint_parameter <- FALSE
          }
        }
      }
    }
    ############################################################################
    ################## Molecular formula annotation criteria ###################
    ############################################################################
    if (PARAM0005 == "yes") {
      ##
      if (NPT > 1) {
        x0015 <- which(PARAM[, 1] == 'PARAM0015')
        parallelizationMode <- PARAM[x0015, 2]
        if (is.na(parallelizationMode)) {
          print("ERROR!!! Problem with PARAM0015!")
          checkpoint_parameter <- FALSE
        } else {
          parallelizationMode <- gsub(" ", "", tolower(parallelizationMode))
          if (parallelizationMode == "samplemode" | parallelizationMode == "peakmode") {
            PARAM[x0015, 2] <- parallelizationMode
          } else {
            print("ERROR!!! Problem with PARAM0015!")
            checkpoint_parameter <- FALSE
          }
        }
      }
      ##
      x0016 <- which(PARAM[, 1] == 'PARAM0016')
      RTtolerance <- tryCatch(as.numeric(PARAM[x0016, 2]), warning = function(w){NA})
      PARAM[x0016, 2] <- RTtolerance
      ##
      if (!is.na(RTtolerance)) {
        x0017 <- which(PARAM[, 1] == 'PARAM0017')
        if (length(x0017) == 0) {
          print("ERROR!!! Problem with PARAM0017!")
          checkpoint_parameter <- FALSE
          PARAM0017 <- "no"
        } else {
          PARAM0017 <- tolower(PARAM[x0017, 2])
          if (!(PARAM0017 == "yes" | PARAM0017 == "no")) {
            print("ERROR!!! Problem with PARAM0017!")
            checkpoint_parameter <- FALSE
          } else {
            PARAM[x0017, 2] <- PARAM0017
          }
        }
        ##
        if (!peak_alignment_folder_check) {
          if (PARAM0017 == "yes") {
            listAlignmentFolderCheck <- IPA_peak_alignment_folder_xlsxAnalyzer(PARAM, PARAM_ID = 'PARAM0012', checkpoint_parameter, correctedRTcheck = TRUE, CSAcheck = FALSE, allowedVerbose = TRUE)
            PARAM <- listAlignmentFolderCheck[[1]]
            checkpoint_parameter <- listAlignmentFolderCheck[[2]]
          }
        }
      }
      ##
      massAccuracy <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0018'), 2])
      if (is.na(massAccuracy)) {
        print("ERROR!!! Problem with PARAM0018!")
        checkpoint_parameter <- FALSE
      } else {
        if (massAccuracy > 0.01) {
          print("ERROR!!! Problem with PARAM0018! Mass accuracy suggested to be below `0.01 Da`")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      maxNEME <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0019'), 2])
      if (is.na(maxNEME)) {
        print("ERROR!!! Problem with PARAM0019!")
        checkpoint_parameter <- FALSE
      } else {
        if (maxNEME < 0) {
          print("ERROR!!! Problem with PARAM0019!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      minPCS <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0020'), 2])
      if (is.na(minPCS)) {
        print("ERROR!!! Problem with PARAM0020!")
        checkpoint_parameter <- FALSE
      } else {
        if (minPCS <= 0) {
          print("ERROR!!! Problem with PARAM0020!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      minNDCS <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0021'), 2])
      if (is.na(minNDCS)) {
        print("ERROR!!! Problem with PARAM0021!")
        checkpoint_parameter <- FALSE
      } else {
        if (minNDCS < 0) {
          print("ERROR!!! Problem with PARAM0021!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      minRCS <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0022'), 2])
      if (is.na(minRCS)) {
        print("ERROR!!! Problem with PARAM0021! This parameter should be between 0-100!")
        checkpoint_parameter <- FALSE
      } else {
        if ((minRCS < 0) | (minRCS > 100)) {
          print("ERROR!!! Problem with PARAM0021! This parameter should be between 0-100!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      scoreCoefficients <- tryCatch(eval(parse(text = PARAM[which(PARAM[, 1] == 'PARAM0023'), 2])), error = function(e){NULL})
      if (is.null(scoreCoefficients)) {
        print("ERROR!!! Problem with PARAM0023! This parameter should be a vector of five positive numbers!")
        checkpoint_parameter <- FALSE
      } else {
        if (length(scoreCoefficients) != 5) {
          print("ERROR!!! Problem with PARAM0023! This parameter should be a vector of five positive numbers!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      maxAllowedNumberHits <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0024'), 2])
      if (is.na(maxAllowedNumberHits)) {
        print("ERROR!!! Problem with PARAM0024!")
        checkpoint_parameter <- FALSE
      } else {
        if (maxAllowedNumberHits < 0) {
          print("ERROR!!! Problem with PARAM0024!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0025 <- which(PARAM[, 1] == 'PARAM0025')
      if (length(x0025) == 0) {
        print("ERROR!!! Problem with PARAM0025!")
        checkpoint_parameter <- FALSE
        PARAM0025 <- "no"
      } else {
        PARAM0025 <- tolower(PARAM[x0025, 2])
        if (!(PARAM0025 == "yes" | PARAM0025 == "no")) {
          print("ERROR!!! Problem with PARAM0025!")
          checkpoint_parameter <- FALSE
        } else {
          PARAM[x0025, 2] <- PARAM0025
        }
      }
      ##
      if (PARAM0025 == "yes") {
        IonPathways <- tryCatch(eval(parse(text = paste0("c(", PARAM[which(PARAM[, 1] == 'PARAM0026'), 2], ")"))), error = function(e){NULL})
        if (is.null(IonPathways)) {
          print("ERROR!!! Problem with PARAM0026!")
          checkpoint_parameter <- FALSE
        }
        ##
        x0027 <- which(PARAM[, 1] == 'PARAM0027')
        if (length(x0027) == 0) {
          print("ERROR!!! Problem with PARAM0027! PubChem library data is not available! You should use the 'molecular_formula_library_generator' module to produce the molecular formula library!")
          checkpoint_parameter <- FALSE
        } else {
          MFlibraryPath <- gsub("\\", "/", PARAM[x0027, 2], fixed = TRUE)
          PARAM[x0027, 2] <- MFlibraryPath
          if (!file.exists(MFlibraryPath)) {
            print("ERROR!!! Problem with PARAM0027! PubChem library data is not available! You should use the 'molecular_formula_library_generator' module to produce the molecular formula library!")
            checkpoint_parameter <- FALSE
          }
        }
      }
    }
    ############################################################################
    ############### Aligned table molecular formula annotation #################
    ############################################################################
    if (PARAM0006 == "yes") {
      maxRankSample <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0028'), 2])
      if (is.na(maxRankSample)) {
        print("ERROR!!! Problem with PARAM0028!")
        checkpoint_parameter <- FALSE
      } else {
        if (maxRankSample <= 0) {
          print("ERROR!!! Problem with PARAM0028!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      Ncandidate <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0029'), 2])
      if (is.na(Ncandidate)) {
        print("ERROR!!! Problem with PARAM0029!")
        checkpoint_parameter <- FALSE
      } else {
        if (Ncandidate <= 0) {
          print("ERROR!!! Problem with PARAM0029!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0030 <- which(PARAM[, 1] == 'PARAM0030')
      if (length(x0030) == 0) {
        print("ERROR!!! Problem with PARAM0030!")
        checkpoint_parameter <- FALSE
        adjustFreqRank <- "no"
      } else {
        adjustFreqRank <- tolower(PARAM[x0030, 2])
        if (!(adjustFreqRank == "yes" | adjustFreqRank == "no")) {
          print("ERROR!!! Problem with PARAM0030!")
          checkpoint_parameter <- FALSE
        } else {
          PARAM[x0030, 2] <- adjustFreqRank
        }
      }
    }
    ############################################################################
    ############################ IPDB production ###############################
    ############################################################################
    if (PARAM0001 == "yes") {
      if (PARAM0002 == "yes") {
        ##
        PARAM_ECS <- UFA_enumerated_chemical_space_xlsxAnalyzer(spreadsheet)
        if (!is.null(PARAM_ECS)) {
          addressIPDB <- PARAM_ECS$`User input 2`[which(PARAM_ECS$Parameter == "IPDB output address")]
          nameIPDB <- PARAM_ECS$`User input 2`[which(PARAM_ECS$Parameter == "IPDB file name")]
          addressIPDB <- paste0(addressIPDB, "/", nameIPDB, ".Rdata")
          addressIPDB <- gsub("\\", "/", addressIPDB, fixed = TRUE)
          PARAM[x0004, 2] <- addressIPDB
        } else {
          checkpoint_parameter <- FALSE
        }
      }
      if (PARAM0003 == "yes") {
        ##
        PARAM_FormSource <- UFA_formula_source_xlsxAnalyzer(spreadsheet)
        if (!is.null(PARAM_FormSource)) {
          addressIPDB <- PARAM_FormSource[which(PARAM_FormSource[, 1] == "FS0002"), 2]
          nameIPDB <- PARAM_FormSource[which(PARAM_FormSource[, 1] == "FS0003"), 2]
          addressIPDB <- paste0(addressIPDB, "/", nameIPDB, ".Rdata")
          addressIPDB <- gsub("\\", "/", addressIPDB, fixed = TRUE)
          PARAM[x0004, 2] <- addressIPDB
        } else {
          checkpoint_parameter <- FALSE
        }
      }
    }
    ############################################################################
    ###################### Score Function Optimization #########################
    ############################################################################
    if (PARAM0007 == "yes") {
      PARAM_ScoreFunc <- UFA_score_function_optimization_xlsxAnalyzer(spreadsheet)
      if (is.null(PARAM_ScoreFunc)) {
        checkpoint_parameter <- FALSE
      }
    }
  }
  ##
  ##############################################################################
  ##
  if (!checkpoint_parameter) {
    listPARAM <- NULL
  } else {
    print("The spreadsheet is consistent with the IDSL.UFA workflow!")
    listPARAM <- list(PARAM, PARAM_ECS, PARAM_FormSource, PARAM_ScoreFunc)
    names(listPARAM) <- c("PARAM", "PARAM_ECS", "PARAM_FormSource", "PARAM_ScoreFunc")
  }
  ##
  return(listPARAM)
}
