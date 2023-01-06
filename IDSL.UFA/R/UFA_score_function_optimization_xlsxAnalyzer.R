UFA_score_function_optimization_xlsxAnalyzer <- function(spreadsheet) {
  print("Initiated testing the `score_function_optimization` tab consistency!")
  ##
  checkpoint_parameter <- FALSE
  if (length(spreadsheet) >= 4) {
    if (typeof(spreadsheet) == "list") {
      PARAM_ScoreFunc <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
    } else {
      print("The UFA spreadsheet was not produced properly!")
    }
  } else if (length(spreadsheet) == 1) {
    if (typeof(spreadsheet) == "character") {
      if (file.exists(spreadsheet)) {
        spreadsheet_UFA <- readxl::read_xlsx(spreadsheet, sheet = "score_function_optimization")
        PARAM_ScoreFunc <- cbind(spreadsheet_UFA[, 2], spreadsheet_UFA[, 4])
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
    ########################## Global parameters ###############################
    ############################################################################
    x0001 <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0001'), 2]
    if (length(x0001) == 0) {
      print("ERROR!!! Problem with SFT0001!")
      checkpoint_parameter <- FALSE
      x0001 <- 0
    } else {
      if (!(tolower(x0001) == "yes" | tolower(x0001) == "no")) {
        print("ERROR!!! Problem with SFT0001!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (tolower(x0001) == "yes") {
      x0002 <- which(PARAM_ScoreFunc[, 1] == 'SFT0002')
      if (length(x0002) > 0) {
        SFT0002 <- PARAM_ScoreFunc[x0002, 2]
        SFT0002 <- gsub("\\", "/", SFT0002, fixed = TRUE)
        PARAM_ScoreFunc[x0002, 2] <- SFT0002
        if (is.na(SFT0002)) {
          print("ERROR!!! Problem with SFT0002! The isotopic profile database file is not available!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    x0003 <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0003'), 2]
    if (length(x0003) == 0) {
      print("ERROR!!! Problem with SFT0003!")
      checkpoint_parameter <- FALSE
      x0003 <- 0
    } else {
      if (!(tolower(x0003) == "yes" | tolower(x0003) == "no")) {
        print("ERROR!!! Problem with SFT0003!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0004 <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0004'), 2]
    if (length(x0004) == 0) {
      print("ERROR!!! Problem with SFT0004!")
      checkpoint_parameter <- FALSE
      x0004 <- 0
    } else {
      if (!(tolower(x0004) == "yes" | tolower(x0004) == "no")) {
        print("ERROR!!! Problem with SFT0004!")
        checkpoint_parameter <- FALSE
      }
    }
    ############################################################################
    ######################## Data import and export ############################
    ############################################################################
    if (tolower(x0001) == "yes") {
      x0005 <- which(PARAM_ScoreFunc[, 1] == 'SFT0005')
      if (length(x0005) == 0) {
        print("ERROR!!! Problem with SFT0005!")
        checkpoint_parameter <- FALSE
      } else {
        input_path_hrms <- PARAM_ScoreFunc[x0005, 2]
        input_path_hrms <- gsub("\\", "/", input_path_hrms, fixed = TRUE)
        PARAM_ScoreFunc[x0005, 2] <- input_path_hrms
        if (!dir.exists(input_path_hrms)) {
          print("ERROR!!! Problem with SFT0005! Please make sure the full path is provided!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      LHRMS <- 0
      excelfile_address <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0006'), 2]
      if (is.na(excelfile_address)) {
        print("Error!!! SFT0006 is empty. Please also check SFT0005!")
        checkpoint_parameter <- FALSE
      } else {
        excelfile_address <- gsub("\\", "/", excelfile_address, fixed = TRUE)
        PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0006'), 2] <- excelfile_address
        if (file.exists(excelfile_address)) {
          PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0006'), 2] <- excelfile_address
          excelfile <- readxl::read_xlsx(excelfile_address)
          col <- colnames(excelfile)
          x_fn <- which(col == 'FileName')
          x_mf <- which(col == 'MolcularFormula')
          x_ipw <- which(col == 'IonizationPathway')
          x_RT <- which(col == 'RetentionTime(min)')
          if (!(length(x_fn) > 0 & length(x_mf) > 0 & length(x_ipw) > 0 & length(x_RT) > 0)) {
            print("ERROR!!! Problem with SFT0006! Incorrect column headers in the reference spreadsheet -> The following columns should be detected in the spreadsheet : 'FileName', 'MolcularFormula', 'IonizationPathway', `RetentionTime(min)` - case sensitive")
            checkpoint_parameter <- FALSE
          } else {
            FileNames <- excelfile$'FileName'
            HRMSfileNames <- unique(FileNames)
            LHRMS <- length(HRMSfileNames)
            ##
            if (LHRMS > 0) {
              xHRMSfileNames <- do.call(c, lapply(HRMSfileNames, function(i) {
                if (!file.exists(paste0(input_path_hrms, "/", i))) {
                  i
                }
              }))
              ##
              if (length(xHRMSfileNames) > 0) {
                print("ERROR!!! Problem with SFT0006! not detected the following HRMS file(s) (case sensitive even for file extensions) in the reference spreadsheet:")
                for (i in xHRMSfileNames) {
                  print(i)
                }
                checkpoint_parameter <- FALSE
              }
            } else {
              print("ERROR!!! Problem with SFT0006! No selected mzML/mzXML/CDF file was detected in the folder!")
              checkpoint_parameter <- FALSE
            }
          }
        } else {
          print("ERROR!!! Problem with SFT0006! The reference spreadsheet not found!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0007 <- which(PARAM_ScoreFunc[, 1] == 'SFT0007')
      if (length(x0007) == 0) {
        print("ERROR!!! Problem with SFT0007!")
        checkpoint_parameter <- FALSE
      } else {
        inputPathPeaklist <- PARAM_ScoreFunc[x0007, 2]
        inputPathPeaklist <- gsub("\\", "/", inputPathPeaklist, fixed = TRUE)
        PARAM_ScoreFunc[x0007, 2] <- inputPathPeaklist
        if (!dir.exists(inputPathPeaklist)) {
          print("ERROR!!! Problem with SFT0007! Please make sure the full path is provided!")
          checkpoint_parameter <- FALSE
        } else {
          ######################################################################
          ## To see if the entire peaklists were generated for all HRMS files ##
          ######################################################################
          if (LHRMS > 0) {
            peaklistFileNames <- dir(path = inputPathPeaklist, pattern = ".Rdata$", ignore.case = TRUE)
            peaklistFileNames <- peaklistFileNames[grep("^peaklist_", peaklistFileNames)]
            L_PL <- length(peaklistFileNames)
            ##
            if (LHRMS > L_PL) {
              checkpoint_parameter <- FALSE
              peaklistHRMSfileNames <- paste0("peaklist_", HRMSfileNames, ".Rdata")
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
    x0008 <- which(PARAM_ScoreFunc[, 1] == 'SFT0008')
    if (length(x0008) == 0) {
      print("ERROR!!! Problem with SFT0008!")
      checkpoint_parameter <- FALSE
    } else {
      output_path <- gsub("\\", "/", PARAM_ScoreFunc[x0008, 2], fixed = TRUE)
      PARAM_ScoreFunc[x0008, 2] <- output_path
      if (!dir.exists(output_path)) {
        tryCatch(dir.create(output_path, recursive = TRUE), warning = function(w){warning("Problem with SFT0008! R cannot create the folder!")})
        if (!dir.exists(output_path)) {
          checkpoint_parameter <- FALSE
        }
      }
    }
    ############################################################################
    ################# Molecular formula annotation criteria ####################
    ############################################################################
    if (tolower(x0003) == "yes") {
      ##
      NPT <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0009'), 2])
      if (length(NPT) == 0) {
        print("ERROR!!! Problem with SFT0009! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      } else {
        if (NPT >= 1) {
          if ((NPT %% 1) != 0) {
            print("ERROR!!! Problem with SFT0009! This parameter should be a positive integer!")
            checkpoint_parameter <- FALSE
          }
        } else {
          print("ERROR!!! Problem with SFT0009! This parameter should be at least 1 !")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      RTtolerance <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0010'), 2])
      if (is.na(RTtolerance)) {
        print("ERROR!!! Problem with SFT0010!")
        checkpoint_parameter <- FALSE
      } else {
        if (RTtolerance < 0) {
          print("ERROR!!! Problem with SFT0010!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      massAccuracy <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0011'), 2])
      if (is.na(massAccuracy)) {
        print("ERROR!!! Problem with SFT0011!")
        checkpoint_parameter <- FALSE
      } else {
        if (massAccuracy > 0.01) {
          print("ERROR!!! Problem with SFT0011! Mass accuracy suggested to be below `0.01 Da`")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      maxNEME <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0012'), 2])
      if (is.na(maxNEME)) {
        print("ERROR!!! Problem with SFT0012!")
        checkpoint_parameter <- FALSE
      } else {
        if (maxNEME < 0) {
          print("ERROR!!! Problem with SFT0012!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      minPCS <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0013'), 2])
      if (is.na(minPCS)) {
        print("ERROR!!! Problem with SFT0013!")
        checkpoint_parameter <- FALSE
      } else {
        if (minPCS <= 0) {
          print("ERROR!!! Problem with SFT0013!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      minNDCS <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0014'), 2])
      if (is.na(minNDCS)) {
        print("ERROR!!! Problem with SFT0014!")
        checkpoint_parameter <- FALSE
      } else {
        if (minNDCS < 0) {
          print("ERROR!!! Problem with SFT0014!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      minRCS <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0015'), 2])
      if (is.na(minRCS)) {
        print("ERROR!!! Problem with SFT0015! This parameter should be between 0-100!")
        checkpoint_parameter <- FALSE
      } else {
        if ((minRCS < 0) | (minRCS > 100)) {
          print("ERROR!!! Problem with SFT0015! This parameter should be between 0-100!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ############################################################################
    ########################### Genetics Algorithm #############################
    ############################################################################
    if ((tolower(x0003) == "yes") | (tolower(x0004) == "yes")) {
      x0016 <- which(PARAM_ScoreFunc[, 1] == 'SFT0016')
      SFT0016 <- PARAM_ScoreFunc[x0016, 2]
      if (is.na(SFT0016)) {
        print("ERROR!!! Problem with SFT0016!")
        checkpoint_parameter <- FALSE
      } else {
        SFT0016 <- gsub(" ", "", tolower(SFT0016))
        if (SFT0016 == "toprank" | SFT0016 == "overalrank") {
          PARAM_ScoreFunc[x0016, 2] <- SFT0016
        } else {
          print("ERROR!!! Problem with SFT0016!")
          checkpoint_parameter <- FALSE
        }
      }
      if (!is.na(SFT0016)) {
        if (SFT0016 == "toprank") {
          maxRank <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0017'), 2])
          if (is.na(maxRank)) {
            print("ERROR!!! Problem with SFT0017! This parameter should be a positive integer greater than or equal to 1 !")
            checkpoint_parameter <- FALSE
          } else {
            if (maxRank <= 0) {
              print("ERROR!!! Problem with SFT0017! This parameter should be a positive integer greater than or equal to 1 !")
              checkpoint_parameter <- FALSE
            } else {
              if ((maxRank %% 1) != 0) {
                print("ERROR!!! Problem with SFT0017! This parameter should be a positive integer greater than or equal to 1 !")
                checkpoint_parameter <- FALSE
              }
            }
          }
        }
      }
    }
    ##
    if (tolower(x0003) == "yes") {
      ##
      GApackageCheck <- tryCatch(requireNamespace('GA', quietly = TRUE), error = function(e) {FALSE})
      if (!GApackageCheck) {
        print("IDSL.UFA requires the 'GA' package of R for the score coefficients optimization workflow!")
        print(" <<< install.packages('GA') >>> ")
        checkpoint_parameter <- FALSE
      }
      ##
      number_processing_threads <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0018'), 2])
      if (length(number_processing_threads) == 0) {
        print("ERROR!!! Problem with SFT0018! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      } else {
        if (number_processing_threads >= 1) {
          if ((number_processing_threads %% 1) != 0) {
            print("ERROR!!! Problem with SFT0018! This parameter should be a positive integer!")
            checkpoint_parameter <- FALSE
          }
        } else {
          print("ERROR!!! Problem with SFT0018! This parameter should be at least 1 !")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      lower_limit <- tryCatch(eval(parse(text = PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0019'), 2])), error = function(e){NULL})
      if (is.null(lower_limit)) {
        print("ERROR!!! Problem with SFT0019! This parameter should be a vector of five positive numbers!")
        checkpoint_parameter <- FALSE
      } else {
        if (length(lower_limit) != 5) {
          print("ERROR!!! Problem with SFT0019! This parameter should be a vector of five positive numbers!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      upper_limit <- tryCatch(eval(parse(text = PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0020'), 2])), error = function(e){NULL})
      if (is.null(upper_limit)) {
        print("ERROR!!! Problem with SFT0020! This parameter should be a vector of five positive numbers!")
        checkpoint_parameter <- FALSE
      } else {
        if (length(upper_limit) != 5) {
          print("ERROR!!! Problem with SFT0020! This parameter should be a vector of five positive numbers!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      diff_limits <- upper_limit - lower_limit
      x_neg <- which(diff_limits < 0)
      if (length(x_neg) > 0) {
        print("ERROR!!! Visit SFT0020 and SFT0019! Upper limits must be greater than lower limits!")
        checkpoint_parameter <- FALSE
      }
      ##
      population_size <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0021'), 2])
      if (is.na(population_size)) {
        print("ERROR!!! Problem with SFT0021!")
        checkpoint_parameter <- FALSE
      } else {
        if (population_size <= 0) {
          print("ERROR!!! Problem with SFT0021!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      max_iteration <- as.numeric(PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0022'), 2])
      if (is.na(max_iteration)) {
        print("ERROR!!! Problem with SFT0022!")
        checkpoint_parameter <- FALSE
      } else {
        if (max_iteration <= 0) {
          print("ERROR!!! Problem with SFT0022!")
          checkpoint_parameter <- FALSE
        }
      }
    }
  }
  if (!checkpoint_parameter) {
    PARAM_ScoreFunc <- NULL
  } else {
    print("The `score_function_optimization` tab is consistent with the score coefficients optimization workflow!")
  }
  ##
  return(PARAM_ScoreFunc)
}
