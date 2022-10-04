UFA_xlsxAnalyzer <- function(spreadsheet) {
  ##
  print("Initiated testing the spreadsheet consistency!")
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
    ######################## Global parameters #################################
    x0001 <- PARAM[which(PARAM[, 1] == 'PARAM0001'), 2]
    if (length(x0001) == 0) {
      print("ERROR!!! Problem with PARAM0001!")
      checkpoint_parameter <- FALSE
      x0001 <- 0
    } else {
      if (!(tolower(x0001) == "yes" | tolower(x0001) == "no")) {
        print("ERROR!!! Problem with PARAM0001!")
        checkpoint_parameter <- FALSE
      }
    }
    if (tolower(x0001) == "yes") {
      x0002 <- PARAM[which(PARAM[, 1] == 'PARAM0002'), 2]
      if (length(x0002) == 0) {
        print("ERROR!!! Problem with PARAM0003!")
        checkpoint_parameter <- FALSE
        x0002 <- 0
      } else {
        if (!(tolower(x0002) == "yes" | tolower(x0002) == "no")) {
          print("ERROR!!! Problem with PARAM0003!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0003 <- PARAM[which(PARAM[, 1] == 'PARAM0003'), 2]
      if (length(x0003) == 0) {
        print("ERROR!!! Problem with PARAM0003!")
        checkpoint_parameter <- FALSE
        x0003 <- 0
      } else {
        if (!(tolower(x0003) == "yes" | tolower(x0003) == "no")) {
          print("ERROR!!! Problem with PARAM0003!")
          checkpoint_parameter <- FALSE
        }
      }
      if ((tolower(x0002) == "yes" & tolower(x0003) == "yes") | (tolower(x0002) == "no" & tolower(x0003) == "no")) {
        x0002 <- 0
        x0003 <- 0
        print("ERROR!!! Problem with PARAM0002 & PARAM0003!")
        print("ERROR!!! Both PARAM0002 & PARAM0003 cannot be 'YES' and 'NO'! Only choose one method to generate isotopic profiles database (IPDB)!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0004 <- which(PARAM[, 1] == 'PARAM0004')
    x0005 <- PARAM[which(PARAM[, 1] == 'PARAM0005'), 2]
    if (tolower(x0001) == "no" & (tolower(x0005) == "yes")) {
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
    if (length(x0005) == 0) {
      print("ERROR!!! Problem with PARAM0005!")
      checkpoint_parameter <- FALSE
    } else {
      if (!(tolower(x0005) == "yes" | tolower(x0005) == "no")) {
        print("ERROR!!! Problem with PARAM0005!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0006 <- PARAM[which(PARAM[, 1] == 'PARAM0006'), 2]
    if (length(x0006) == 0) {
      print("ERROR!!! Problem with PARAM0006!")
      checkpoint_parameter <- FALSE
    } else {
      if (!(tolower(x0006) == "yes" | tolower(x0006) == "no")) {
        print("ERROR!!! Problem with PARAM0006!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0007 <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
    if (length(x0007) == 0) {
      print("ERROR!!! Problem with PARAM0007!")
      checkpoint_parameter <- FALSE
    } else {
      if (!(tolower(x0007) == "yes" | tolower(x0007) == "no")) {
        print("ERROR!!! Problem with PARAM0007!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0008 <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
    if (length(x0008) == 0) {
      print("ERROR!!! Problem with PARAM0008!")
      checkpoint_parameter <- FALSE
    } else {
      if (!(tolower(x0008) == "yes" | tolower(x0008) == "no")) {
        print("ERROR!!! Problem with PARAM0008!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0009 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2])
    if (length(x0009) == 0) {
      print("ERROR!!! Problem with PARAM0009! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (x0009 >= 1) {
        if ((x0009 %% 1) != 0) {
          print("ERROR!!! Problem with PARAM0009! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        print("ERROR!!! Problem with PARAM0009! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ################################## Data ####################################
    if (tolower(x0005) == "yes") {
      x0010 <- which(PARAM[, 1] == 'PARAM0010')
      if (length(x0010) == 0) {
        print("ERROR!!! Problem with PARAM0010!")
        checkpoint_parameter <- FALSE
      } else {
        input_path_hrms <- PARAM[x0010, 2]
        input_path_hrms <- gsub("\\", "/", input_path_hrms, fixed = TRUE)
        PARAM[x0010, 2] <- input_path_hrms
        if (!dir.exists(input_path_hrms)) {
          print("ERROR!!! Problem with PARAM0010! Please make sure the full path is provided!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0011 <- which(PARAM[, 1] == 'PARAM0011')
      if (is.null(PARAM[x0011, 2])) {
        print("ERROR!!! Problem with PARAM0011!")
        checkpoint_parameter <- FALSE
      } else {
        if (tolower(PARAM[x0011, 2]) != "all") {
          samples_string <- PARAM[x0011, 2]
          name <- strsplit(samples_string, ";")[[1]]
          ID <- sapply(1:length(name), function(i) {
            ID_name <- paste0(input_path_hrms, "/", name[i])
            as.numeric(file.exists(ID_name))
          })
          x_ID <- which(ID == 0)
          if (length(x_ID) > 0) {
            print("ERROR!!! Problem with PARAM0011! not detected the following file(s) (case sensitive even for file extensions):")
            for (i in 1:length(x_ID)) {
              print(name[x_ID[i]])
            }
            checkpoint_parameter <- FALSE
          }
        }
        ##
        if (tolower(PARAM[x0011, 2]) == "all") {
          x0012 <- PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]
          if (is.na(x0012)) {
            print("ERROR!!! Problem with PARAM0012!")
            checkpoint_parameter <- FALSE
          } else {
            if (!(tolower(x0012) == "mzml" | tolower(x0012) == "mzxml" | tolower(x0012) == "cdf")) {
              print("ERROR!!! Problem with PARAM0012! HRMS data are incompatible!")
              checkpoint_parameter <- FALSE
            }
          }
        }
        ##
        x0013 <- which(PARAM[, 1] == 'PARAM0013')
        if (length(x0013) == 0) {
          print("ERROR!!! Problem with PARAM0013!")
          checkpoint_parameter <- FALSE
        } else {
          input_path_pl <- PARAM[x0013, 2]
          input_path_pl <- gsub("\\", "/", input_path_pl, fixed = TRUE)
          PARAM[x0013, 2] <- input_path_pl
          if (!dir.exists(input_path_pl)) {
            print("ERROR!!! Problem with PARAM0013! Please make sure the full path is provided!")
            checkpoint_parameter <- FALSE
          }
        }
      }
      #### To see if the entire peaklists were generated for all HRMS files ####
      if (dir.exists(input_path_pl)) {
        file_names_peaklist1 <- dir(path = input_path_pl, pattern = ".Rdata")
        file_names_peaklist2 <- dir(path = input_path_pl, pattern = "peaklist_")
        file_names_peaklist <- file_names_peaklist1[file_names_peaklist1 %in% file_names_peaklist2]
        ##
        if (dir.exists(input_path_hrms)) {
          if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0011'), 2]) == "all") {
            file_name_hrms <- dir(path = input_path_hrms)
            file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
          } else {
            samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0011'), 2]
            file_name_hrms <- strsplit(samples_string, ";")[[1]]
          }
          ##
          file_names_peaklist_hrms1 <- gsub(".Rdata", "", file_names_peaklist)
          file_names_peaklist_hrms2 <- gsub("peaklist_", "", file_names_peaklist_hrms1)
          file_names_peaklist_hrms <- file_name_hrms %in% file_names_peaklist_hrms2
          L_PL <- length(which(file_names_peaklist_hrms == TRUE))
          if (length(file_name_hrms) != L_PL) {
            checkpoint_parameter <- FALSE
            print("Error!!! peaklist files are not available for the selected HRMS files!")
          }
        }
      }
    }
    if ((tolower(x0005) == "yes") | (tolower(x0006) == "yes")) {
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
    ################## Molecular formula annotation criteria ###################
    ############################################################################
    if (tolower(x0005) == "yes") {
      ##
      x0015 <- which(PARAM[, 1] == 'PARAM0015')
      PARAM0015 <- PARAM[x0015, 2]
      if (is.na(PARAM0015)) {
        print("ERROR!!! Problem with PARAM0015!")
        checkpoint_parameter <- FALSE
      } else {
        PARAM0015 <- gsub(" ", "", tolower(PARAM0015))
        if (PARAM0015 == "samplemode" | PARAM0015 == "peakmode") {
          PARAM[x0015, 2] <- PARAM0015
        } else {
          print("ERROR!!! Problem with PARAM0015!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0016 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0016'), 2])
      if (is.na(x0016)) {
        print("ERROR!!! Problem with PARAM0016!")
        checkpoint_parameter <- FALSE
      } else {
        if (x0016 < 0) {
          print("ERROR!!! Problem with PARAM0016!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0017 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0017'), 2])
      if (is.na(x0017)) {
        print("ERROR!!! Problem with PARAM0017!")
        checkpoint_parameter <- FALSE
      } else {
        if (x0017 <= 0) {
          print("ERROR!!! Problem with PARAM0017!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0018 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0018'), 2])
      if (is.na(x0018)) {
        print("ERROR!!! Problem with PARAM0018!")
        checkpoint_parameter <- FALSE
      } else {
        if (x0018 <= 0) {
          print("ERROR!!! Problem with PARAM0018!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0019 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0019'), 2])
      if (is.na(x0019)) {
        print("ERROR!!! Problem with PARAM0019!")
        checkpoint_parameter <- FALSE
      } else {
        if (x0019 < 0) {
          print("ERROR!!! Problem with PARAM0019!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0020 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0020'), 2])
      if (is.na(x0020)) {
        print("ERROR!!! Problem with PARAM0020!")
        checkpoint_parameter <- FALSE
      } else {
        if (x0020 <= 0) {
          print("ERROR!!! Problem with PARAM0020!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0021 <- tryCatch(eval(parse(text = PARAM[which(PARAM[, 1] == 'PARAM0021'), 2])), error = function(e){NULL})
      if (is.null(x0021)) {
        print("ERROR!!! Problem with PARAM0021! This parameter should be a vector of five positive numbers!")
        checkpoint_parameter <- FALSE
      } else {
        if (length(x0021) != 5) {
          print("ERROR!!! Problem with PARAM0021! This parameter should be a vector of five positive numbers!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0022 <- PARAM[which(PARAM[, 1] == 'PARAM0022'), 2]
      if (length(x0022) == 0) {
        print("ERROR!!! Problem with PARAM0022!")
        checkpoint_parameter <- FALSE
      } else {
        if (!(tolower(x0022) == "yes" | tolower(x0022) == "no")) {
          print("ERROR!!! Problem with PARAM0022!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      if (tolower(x0022) == "yes") {
        x0023 <- tryCatch(eval(parse(text = paste0("c(", PARAM[which(PARAM[, 1] == 'PARAM0023'), 2], ")"))), error = function(e){NULL})
        if (is.null(x0023)) {
          print("ERROR!!! Problem with PARAM0023!")
          checkpoint_parameter <- FALSE
        }
        ##
        x0024 <- which(PARAM[, 1] == 'PARAM0024')
        if (length(x0024) == 0) {
          print("ERROR!!! Problem with PARAM0024! PubChem library data is not available! You should use the 'molecular_formula_library_generator' module to produce the molecular formula library!")
          checkpoint_parameter <- FALSE
        } else {
          PubChem_library_path <- gsub("\\", "/", PARAM[x0024, 2], fixed = TRUE)
          PARAM[x0024, 2] <- PubChem_library_path
          if (!file.exists(PubChem_library_path)) {
            print("ERROR!!! Problem with PARAM0024! PubChem library data is not available! You should use the 'molecular_formula_library_generator' module to produce the molecular formula library!")
            checkpoint_parameter <- FALSE
          }
        }
      }
    }
    ############### Aligned table molecular formula annotation #################
    if (tolower(x0006) == "yes") {
      x0025 <- which(PARAM[, 1] == 'PARAM0025')
      if (length(x0025) == 0) {
        print("ERROR!!! Problem with PARAM0025! The aligned indixed peak table is not available! This file usually is named 'peak_Xcol.Rdata'")
        checkpoint_parameter <- FALSE
      } else {
        ipa_Xcol_path <- gsub("\\", "/", PARAM[x0025, 2], fixed = TRUE)
        PARAM[x0025, 2] <- ipa_Xcol_path
        if (!file.exists(ipa_Xcol_path)) {
          print("ERROR!!! Problem with PARAM0025! The aligned indixed peak table is not available! This file usually is named 'peak_Xcol.Rdata'")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0026 <- which(PARAM[, 1] == 'PARAM0026')
      if (length(x0026) == 0) {
        print("ERROR!!! Problem with PARAM0026! The aligned peak property table is not available!")
        checkpoint_parameter <- FALSE
      } else {
        ipa_height_path <- gsub("\\", "/", PARAM[x0026, 2], fixed = TRUE)
        PARAM[x0026, 2] <- ipa_height_path
        if (!file.exists(ipa_height_path)) {
          print("ERROR!!! Problem with PARAM0026! The aligned peak property table is not available!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0027 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0027'), 2])
      if (is.na(x0027)) {
        print("ERROR!!! Problem with PARAM0027!")
        checkpoint_parameter <- FALSE
      } else {
        if (x0027 <= 0) {
          print("ERROR!!! Problem with PARAM0027!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0028 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0028'), 2])
      if (is.na(x0028)) {
        print("ERROR!!! Problem with PARAM0028!")
        checkpoint_parameter <- FALSE
      } else {
        if (x0028 <= 0) {
          print("ERROR!!! Problem with PARAM0028!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0029 <- PARAM[which(PARAM[, 1] == 'PARAM0029'), 2]
      if (length(x0029) == 0) {
        print("ERROR!!! Problem with PARAM0029!")
        checkpoint_parameter <- FALSE
      } else {
        if (!(tolower(x0029) == "yes" | tolower(x0029) == "no")) {
          print("ERROR!!! Problem with PARAM0029!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ############################ IPDB production ###############################
    if (tolower(x0001) == "yes") {
      if (tolower(x0002) == "yes") {
        print("Initiated analyzing the `enumerated_chemical_space` tab!")
        ##
        PARAM_MF <- readxl::read_xlsx(spreadsheet, sheet = "enumerated_chemical_space")
        x_address_IPDB <- which(PARAM_MF$Parameter == "IPDB output address")
        x_name_IPDB <- which(PARAM_MF$Parameter == "IPDB file name")
        address_IPDB <- paste0(PARAM_MF$`User input 2`[x_address_IPDB], "/", PARAM_MF$`User input 2`[x_name_IPDB], ".Rdata")
        address_IPDB <- gsub("\\", "/", address_IPDB, fixed = TRUE)
        PARAM[x0004, 2] <- address_IPDB
        ##
        checkpoint_parameter <- UFA_enumerated_chemical_space_xlsxAnalyzer(PARAM_MF)
        print("Completed analyzing the `enumerated_chemical_space` tab!")
      }
      if (tolower(x0003) == "yes") {
        print("Initiated analyzing the `formula_source` tab!")
        PARAM_SF <- readxl::read_xlsx(spreadsheet, sheet = "formula_source")
        PARAM_SF <- cbind(PARAM_SF[, 2], PARAM_SF[, 4])
        ##
        fs0001 <- which(PARAM_SF[, 1] == 'FS0001')
        if (length(fs0001) == 0) {
          print("ERROR!!! Problem with FS0001! Molecular formula source file is not available!")
          checkpoint_parameter <- FALSE
        } else {
          Molecular_formula_source_file <- gsub("\\", "/", PARAM_SF[fs0001, 2], fixed = TRUE)
          PARAM_SF[fs0001, 2] <- Molecular_formula_source_file
          if (!file.exists(Molecular_formula_source_file)) {
            print("ERROR!!! Problem with FS0001! Molecular formula source file is not available!")
            checkpoint_parameter <- FALSE
          }
        }
        ##
        fs0002 <- as.numeric(PARAM_SF[which(PARAM_SF[, 1] == 'FS0002'), 2])
        if (is.na(fs0002)) {
          print("ERROR!!! Problem with FS0002!")
          checkpoint_parameter <- FALSE
        } else {
          if (fs0002 < 0) {
            print("ERROR!!! Problem with FS0002!")
            checkpoint_parameter <- FALSE
          }
        }
        ##
        fs0003 <- PARAM_SF[which(PARAM_SF[, 1] == 'FS0003'), 2]
        if (is.na(fs0003)) {
          print("ERROR!!! Problem with FS0003!")
          checkpoint_parameter <- FALSE
        }
        ##
        fs0004 <- tryCatch(eval(parse(text = paste0("c(", PARAM_SF[which(PARAM_SF[, 1] == 'FS0004'), 2], ")"))), error = function(e){NULL})
        if (length(fs0004) != 3) {
          print("ERROR!!! Problem with FS0004! This parameter should be a vector of three positive numbers!")
          checkpoint_parameter <- FALSE
        }
        ##
        fs0005 <- tryCatch(eval(parse(text = paste0("c(", PARAM_SF[which(PARAM_SF[, 1] == 'FS0005'), 2], ")"))), error = function(e){NULL})
        if (is.null(fs0005)) {
          print("ERROR!!! Problem with FS0005!")
          checkpoint_parameter <- FALSE
        }
        ##
        fs0006 <- as.numeric(PARAM_SF[which(PARAM_SF[, 1] == 'FS0006'), 2])
        if (length(fs0006) == 0) {
          print("ERROR!!! Problem with FS0006! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        } else {
          if (fs0006 >= 1) {
            if ((fs0006 %% 1) != 0) {
              print("ERROR!!! Problem with FS0005! This parameter should be a positive integer!")
              checkpoint_parameter <- FALSE
            }
          } else {
            print("ERROR!!! Problem with FS0005! This parameter should be at least 1 !")
            checkpoint_parameter <- FALSE
          }
        }
        ##
        x_address <- PARAM_SF[which(PARAM_SF[, 1] == "FS0007"), 2]
        x_file <- PARAM_SF[which(PARAM_SF[, 1] == "FS0008"), 2]
        address_IPDB <- paste0(x_address, "/", x_file, ".Rdata")
        address_IPDB <- gsub("\\", "/", address_IPDB, fixed = TRUE)
        PARAM[x0004, 2] <- address_IPDB
        ##
        print("Completed analyzing the `formula_source` tab!")
      }
    }
    if (tolower(x0007) == "yes") {
      PARAM_SFT <- UFA_score_function_optimization_xlsxAnalyzer(spreadsheet)
      if (is.null(PARAM_SFT)) {
        checkpoint_parameter <- FALSE
      }
    }
    if (tolower(x0008) == "yes") {
      PARAM_SA <- UFA_profile_visualizer_xlsxAnalyzer(spreadsheet)
      if (is.null(PARAM_SA)) {
        checkpoint_parameter <- FALSE
      }
    }
  }
  if (!checkpoint_parameter) {
    print("Please visit   https://ufa.idsl.me    for instructions!")
    PARAM <- NULL
  } else {
    print("The spreadsheet is consistent with the IDSL.UFA workflow!")
  }
  return(PARAM)
}
