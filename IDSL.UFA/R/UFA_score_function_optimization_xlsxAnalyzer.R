UFA_score_function_optimization_xlsxAnalyzer <- function(spreadsheet) {
  checkpoint_parameter <- 0
  if (length(spreadsheet) >= 4) {
    if (typeof(spreadsheet) == "list") {
      PARAM_SFT <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
    } else {
      print("The UFA spreadsheet was not produced properly!")
    }
  } else if (length(spreadsheet) == 1) {
    if (typeof(spreadsheet) == "character") {
      if (file.exists(spreadsheet)) {
        spreadsheet_UFA <- readxl::read_xlsx(spreadsheet, sheet = "score_function_optimization")
        PARAM_SFT <- cbind(spreadsheet_UFA[, 2], spreadsheet_UFA[, 4])
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
  if (checkpoint_parameter == TRUE) {
    ########################## Global parameters #################################
    x0001 <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0001'), 2]
    if (length(x0001) == 0) {
      print("ERROR!!! Problem with SFT0001!")
      checkpoint_parameter <- 0
      x0001 <- 0
    } else {
      if (!(tolower(x0001) == "yes" | tolower(x0001) == "no")) {
        print("ERROR!!! Problem with SFT0001!")
        checkpoint_parameter <- 0
      }
    }
    ##
    if (tolower(x0001) == "yes") {
      x0002 <- which(PARAM_SFT[, 1] == 'SFT0002')
      if (length(x0002) > 0) {
        SFT0002 <- PARAM_SFT[x0002, 2]
        SFT0002 <- gsub("\\", "/", SFT0002, fixed = TRUE)
        PARAM_SFT[x0002, 2] <- SFT0002
        if (is.na(SFT0002)) {
          print("ERROR!!! Problem with SFT0002! The isotopic profile database file is not available!")
          checkpoint_parameter <- 0
        }
      }
    }
    ##
    x0003 <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0003'), 2]
    if (length(x0003) == 0) {
      print("ERROR!!! Problem with SFT0003!")
      checkpoint_parameter <- 0
      x0003 <- 0
    } else {
      if (!(tolower(x0003) == "yes" | tolower(x0003) == "no")) {
        print("ERROR!!! Problem with SFT0003!")
        checkpoint_parameter <- 0
      }
    }
    ##
    x0004 <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0004'), 2]
    if (length(x0004) == 0) {
      print("ERROR!!! Problem with SFT0004!")
      checkpoint_parameter <- 0
      x0004 <- 0
    } else {
      if (!(tolower(x0004) == "yes" | tolower(x0004) == "no")) {
        print("ERROR!!! Problem with SFT0004!")
        checkpoint_parameter <- 0
      }
    }
    ##
    x0005 <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0005'), 2])
    if (length(x0005) == 0) {
      print("ERROR!!! Problem with SFT0005! This parameter should be a positive integer!")
      checkpoint_parameter <- 0
    } else {
      if (x0005 >= 1) {
        if ((x0005 %% 1) != 0) {
          print("ERROR!!! Problem with SFT0005! This parameter should be a positive integer!")
          checkpoint_parameter <- 0
        }
      } else {
        print("ERROR!!! Problem with SFT0005! This parameter should be at least 1 !")
        checkpoint_parameter <- 0
      }
    }
    ################################## Data ####################################
    if (!(tolower(x0001) == "no"  & tolower(x0003) == "no")) {
      x0006 <- which(PARAM_SFT[, 1] == 'SFT0006')
      if (length(x0006) == 0) {
        print("ERROR!!! Problem with SFT0006!")
        checkpoint_parameter <- 0
      } else {
        input_path_hrms <- PARAM_SFT[x0006, 2]
        input_path_hrms <- gsub("\\", "/", input_path_hrms, fixed = TRUE)
        PARAM_SFT[x0006, 2] <- input_path_hrms
        if (!dir.exists(input_path_hrms)) {
          print("ERROR!!! Problem with SFT0006! Please make sure the full path is provided!")
          checkpoint_parameter <- 0
        }
      }
      ##
      x0007 <- which(PARAM_SFT[, 1] == 'SFT0007')
      if (is.null(PARAM_SFT[x0007, 2])) {
        print("ERROR!!! Problem with SFT0007!")
        checkpoint_parameter <- 0
      } else {
        if (tolower(PARAM_SFT[x0007, 2]) != "all") {
          samples_string <- PARAM_SFT[x0007, 2]
          name <- strsplit(samples_string, ";")[[1]]
          ID <- sapply(1:length(name), function(i) {
            ID_name <- paste0(input_path_hrms, "/", name[i])
            as.numeric(file.exists(ID_name))
          })
          x_ID <- which(ID == 0)
          if (length(x_ID) > 0) {
            print("ERROR!!! Problem with SFT0007! not detected the following file(s) (case sensitive even for file extensions):")
            for (i in 1:length(x_ID)) {
              print(name[x_ID[i]])
            }
            checkpoint_parameter <- 0
          }
        }
        ##
        if (tolower(PARAM_SFT[x0007, 2]) == "all") {
          x0008 <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0008'), 2]
          if (is.na(x0008)) {
            print("ERROR!!! Problem with SFT0008!")
            checkpoint_parameter <- 0
          } else {
            if (!(tolower(x0008) == "mzml" | tolower(x0008) == "mzxml" | tolower(x0008) == "cdf")) {
              print("ERROR!!! Problem with SFT0008! HRMS data are incompatible!")
              checkpoint_parameter <- 0
            }
          }
        }
        ##
        x0009 <- which(PARAM_SFT[, 1] == 'SFT0009')
        if (length(x0009) == 0) {
          print("ERROR!!! Problem with SFT0009!")
          checkpoint_parameter <- 0
        } else {
          input_path_pl <- PARAM_SFT[x0009, 2]
          input_path_pl <- gsub("\\", "/", input_path_pl, fixed = TRUE)
          PARAM_SFT[x0009, 2] <- input_path_pl
          if (!dir.exists(input_path_pl)) {
            print("ERROR!!! Problem with SFT0009! Please make sure the full path is provided!")
            checkpoint_parameter <- 0
          }
        }
      }
      ##### To see if the entire peaklists were generated for all HRMS files #####
      if (dir.exists(input_path_pl)) {
        file_names_peaklist1 <- dir(path = input_path_pl, pattern = ".Rdata")
        file_names_peaklist2 <- dir(path = input_path_pl, pattern = "peaklist_")
        file_names_peaklist <- file_names_peaklist1[file_names_peaklist1%in%file_names_peaklist2]
        L_PL <- length(file_names_peaklist)
        ##
        if (dir.exists(input_path_hrms)) {
          if (tolower(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0007'), 2]) == "all") {
            file_name_hrms <- dir(path = input_path_hrms)
            file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0008'), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
          } else {
            samples_string <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0007'), 2]
            file_name_hrms <- strsplit(samples_string, ";")[[1]]
          }
          ##
          file_names_peaklist_hrms1 <- gsub(".Rdata", "", file_names_peaklist)
          file_names_peaklist_hrms2 <- gsub("peaklist_", "", file_names_peaklist_hrms1)
          file_names_peaklist_hrms <- file_name_hrms%in%file_names_peaklist_hrms2
          if (length(which(file_names_peaklist_hrms == TRUE)) != L_PL) {
            checkpoint_parameter <- 0
            print("Error!!! peaklist files are not available for the selected HRMS files!")
          }
        }
      }
      ##
      x0010 <- which(PARAM_SFT[, 1] == 'SFT0010')
      if (length(x0010) == 0) {
        print("ERROR!!! Problem with SFT0010!")
        checkpoint_parameter <- 0
      } else {
        output_path <- gsub("\\", "/", PARAM_SFT[x0010, 2], fixed = TRUE)
        PARAM_SFT[x0010, 2] <- output_path
        if (!dir.exists(output_path)) {
          tryCatch(dir.create(output_path, recursive = TRUE), warning = function(w){warning("Problem with SFT0010! R cannot create the folder!")})
          if (!dir.exists(output_path)) {
            checkpoint_parameter <- 0
          }
        }
      }
    }
    ################# Molecular formula annotation criteria ######################
    ##############################################################################
    if (tolower(x0003) == "yes") {
      ref_xlsx_file <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0011'), 2]
      if (is.na(ref_xlsx_file)) {
        print("Error!!! SFT0011 is empty. Please also check SFT0006!")
        checkpoint_parameter <- 0
      } else {
        ref_xlsx_file <- gsub("\\", "/", ref_xlsx_file, fixed = TRUE)
        PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0011'), 2] <- ref_xlsx_file
        if (file.exists(ref_xlsx_file)) {
          PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0011'), 2] <- ref_xlsx_file
          ref_table <- readxl::read_xlsx(ref_xlsx_file)
          col <- colnames(ref_table)
          x_fn <- which(col == 'FileName')
          x_mf <- which(col == 'MolcularFormula')
          x_ipw <- which(col == 'IonizationPathway')
          x_RT <- which(col == 'RetentionTime(min)')
          if (!(length(x_fn) > 0 & length(x_mf) > 0 & length(x_ipw) > 0 & length(x_RT) > 0)) {
            print("ERROR!!! Problem with SFT0011! Incorrect column headers in the reference spreadsheet -> The following columns should be detected in the spreadsheet : 'FileName', 'MolcularFormula', 'IonizationPathway', `RetentionTime(min)` - case sensitive")
            checkpoint_parameter <- 0
          }
        } else {
          print("ERROR!!! Problem with SFT0011! The reference spreadsheet not found!")
          checkpoint_parameter <- 0
        }
      }
      ##
      x0012 <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0012'), 2])
      if (is.na(x0012)) {
        print("ERROR!!! Problem with SFT0012!")
        checkpoint_parameter <- 0
      } else {
        if (x0012 < 0) {
          print("ERROR!!! Problem with SFT0012!")
          checkpoint_parameter <- 0
        }
      }
      ##
      x0013 <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0013'), 2])
      if (is.na(x0013)) {
        print("ERROR!!! Problem with SFT0013!")
        checkpoint_parameter <- 0
      } else {
        if (x0013 < 0) {
          print("ERROR!!! Problem with SFT0013!")
          checkpoint_parameter <- 0
        }
      }
      ##
      x0014 <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0014'), 2])
      if (is.na(x0014)) {
        print("ERROR!!! Problem with SFT0014!")
        checkpoint_parameter <- 0
      } else {
        if (x0014 <= 0) {
          print("ERROR!!! Problem with SFT0014!")
          checkpoint_parameter <- 0
        }
      }
      ##
      x0015 <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0015'), 2])
      if (is.na(x0015)) {
        print("ERROR!!! Problem with SFT0015!")
        checkpoint_parameter <- 0
      } else {
        if (x0015 <= 0) {
          print("ERROR!!! Problem with SFT0015!")
          checkpoint_parameter <- 0
        }
      }
      ##
      x0016 <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0016'), 2])
      if (is.na(x0016)) {
        print("ERROR!!! Problem with SFT0016!")
        checkpoint_parameter <- 0
      } else {
        if (x0016 < 0) {
          print("ERROR!!! Problem with SFT0016!")
          checkpoint_parameter <- 0
        }
      }
      ##
      x0017 <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0017'), 2])
      if (is.na(x0017)) {
        print("ERROR!!! Problem with SFT0017!")
        checkpoint_parameter <- 0
      } else {
        if (x0017 <= 0) {
          print("ERROR!!! Problem with SFT0017!")
          checkpoint_parameter <- 0
        }
      }
    }
    ########################### Genetics Algorithm #############################
    ############################################################################
    SFT_function_obj_func_test <- function(PARAM_SFT, checkpoint_parameter) {
      x0018 <- which(PARAM_SFT[, 1] == 'SFT0018')
      SFT0018 <- PARAM_SFT[x0018, 2]
      if (is.na(SFT0018)) {
        print("ERROR!!! Problem with SFT0018!")
        checkpoint_parameter <- 0
      } else {
        SFT0018 <- gsub(" ", "", tolower(SFT0018))
        if (SFT0018 == "toprank" | SFT0018 == "overalrank") {
          PARAM_SFT[x0018, 2] <- SFT0018
        } else {
          print("ERROR!!! Problem with SFT0018!")
          checkpoint_parameter <- 0
        }
      }
      if (!is.na(SFT0018)) {
        if (SFT0018 == "toprank") {
          max_rank <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0019'), 2])
          if (is.na(max_rank)) {
            print("ERROR!!! Problem with SFT0019! This parameter should be a positive integer greater than or equal to 1 !")
            checkpoint_parameter <- 0
          } else {
            if (max_rank <= 0) {
              print("ERROR!!! Problem with SFT0019! This parameter should be a positive integer greater than or equal to 1 !")
              checkpoint_parameter <- 0
            } else {
              if ((max_rank %% 1) != 0) {
                print("ERROR!!! Problem with SFT0019! This parameter should be a positive integer greater than or equal to 1 !")
                checkpoint_parameter <- 0
              }
            }
          }
        }
      }
      return(checkpoint_parameter)
    }
    ##
    if (tolower(x0003) == "yes") {
      checkpoint_parameter <- SFT_function_obj_func_test(PARAM_SFT, checkpoint_parameter)
      ##
      SFT0020 <- tryCatch(eval(parse(text = PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0020'), 2])), error = function(e){NULL})
      if (is.null(SFT0020)) {
        print("ERROR!!! Problem with SFT0020! This parameter should be a vector of five positive numbers!")
        checkpoint_parameter <- 0
      } else {
        if (length(SFT0020) != 5) {
          print("ERROR!!! Problem with SFT0020! This parameter should be a vector of five positive numbers!")
          checkpoint_parameter <- 0
        }
      }
      ##
      SFT0021 <- tryCatch(eval(parse(text = PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0021'), 2])), error = function(e){NULL})
      if (is.null(SFT0021)) {
        print("ERROR!!! Problem with SFT0021! This parameter should be a vector of five positive numbers!")
        checkpoint_parameter <- 0
      } else {
        if (length(SFT0021) != 5) {
          print("ERROR!!! Problem with SFT0021! This parameter should be a vector of five positive numbers!")
          checkpoint_parameter <- 0
        }
      }
      ##
      diff_limits <- SFT0021 - SFT0020
      x_neg <- which(diff_limits < 0)
      if (length(x_neg) > 0) {
        print("ERROR!!! Visit SFT0021 and SFT0020! Upper limits must be greater than lower limits!")
        checkpoint_parameter <- 0
      }
      ##
      x0022 <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0022'), 2])
      if (is.na(x0022)) {
        print("ERROR!!! Problem with SFT0022!")
        checkpoint_parameter <- 0
      } else {
        if (x0022 <= 0) {
          print("ERROR!!! Problem with SFT0022!")
          checkpoint_parameter <- 0
        }
      }
      ##
      x0023 <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0023'), 2])
      if (is.na(x0023)) {
        print("ERROR!!! Problem with SFT0023!")
        checkpoint_parameter <- 0
      } else {
        if (x0023 <= 0) {
          print("ERROR!!! Problem with SFT0023!")
          checkpoint_parameter <- 0
        }
      }
    }
    ############################################################################
    ####################### Genetics Algorithm Evaluation ######################
    if (tolower(x0004) == "yes") {
      x0014 <- as.numeric(PARAM_SFT[which(PARAM_SFT[, 1] == "SFT0014"), 2])
      if (is.na(x0014)) {
        print("ERROR!!! Problem with SFT0014!")
        checkpoint_parameter <- 0
      }
      else {
        if (x0014 <= 0) {
          print("ERROR!!! Problem with SFT0014!")
          checkpoint_parameter <- 0
        }
      }
      ##
      checkpoint_parameter <- SFT_function_obj_func_test(PARAM_SFT, checkpoint_parameter)
    }
    ##
  }
  if (checkpoint_parameter == FALSE) {
    print("Please visit   https://ufa.idsl.me   for instructions!")
    PARAM_SFT <- NULL
  }
  return(PARAM_SFT)
}
