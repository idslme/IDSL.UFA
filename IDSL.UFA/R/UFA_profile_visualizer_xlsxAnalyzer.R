UFA_profile_visualizer_xlsxAnalyzer <- function (spreadsheet) {
  checkpoint_parameter <- FALSE
  if (length(spreadsheet) >= 4) {
    if (typeof(spreadsheet) == "list") {
      PARAM_SA <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
    } else {
      print("The UFA spreadsheet was not produced properly!")
    }
  } else if (length(spreadsheet) == 1) {
    if (typeof(spreadsheet) == "character") {
      if (file.exists(spreadsheet)) {
        spreadsheet_SA <- readxl::read_xlsx(spreadsheet, sheet = "profile_visualization")
        PARAM_SA <- cbind(spreadsheet_SA[, 2], spreadsheet_SA[, 4])
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
  ##############################################################################
  if (checkpoint_parameter) {
    x0010 <- which(PARAM_SA[, 1] == "PARAM0010")
    if (length(x0010) == 0) {
      print("ERROR!!! Problem with PARAM0010!")
      checkpoint_parameter <- FALSE
    } else {
      input_path_hrms <- PARAM_SA[x0010, 2]
      input_path_hrms <- gsub("\\", "/", input_path_hrms, fixed = TRUE)
      PARAM_SA[x0010, 2] <- input_path_hrms
      if (!dir.exists(input_path_hrms)) {
        print("ERROR!!! Problem with PARAM0010! Please make sure the full path is provided!")
        checkpoint_parameter <- FALSE
      }
    }
    x0011 <- which(PARAM_SA[, 1] == "PARAM0011")
    if (is.null(PARAM_SA[x0011, 2])) {
      print("ERROR!!! Problem with PARAM0011!")
      checkpoint_parameter <- FALSE
    } else {
      if (tolower(PARAM_SA[x0011, 2]) != "all") {
        samples_string <- PARAM_SA[x0011, 2]
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
      if (tolower(PARAM_SA[x0011, 2]) == "all") {
        x0012 <- PARAM_SA[which(PARAM_SA[, 1] == "PARAM0012"), 2]
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
      x0013 <- which(PARAM_SA[, 1] == "PARAM0013")
      if (length(x0013) == 0) {
        print("ERROR!!! Problem with PARAM0013!")
        checkpoint_parameter <- FALSE
      } else {
        input_path_pl <- PARAM_SA[x0013, 2]
        input_path_pl <- gsub("\\", "/", input_path_pl, fixed = TRUE)
        PARAM_SA[x0013, 2] <- input_path_pl
        if (!dir.exists(input_path_pl)) {
          print("ERROR!!! Problem with PARAM0013! Please make sure the full path is provided!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    if (dir.exists(input_path_pl)) {
      file_names_peaklist1 <- dir(path = input_path_pl, pattern = ".Rdata")
      file_names_peaklist2 <- dir(path = input_path_pl, pattern = "peaklist_")
      file_names_peaklist <- file_names_peaklist1[file_names_peaklist1 %in% file_names_peaklist2]
      if (dir.exists(input_path_hrms)) {
        if (tolower(PARAM_SA[which(PARAM_SA[, 1] == "PARAM0011"), 2]) == "all") {
          file_name_hrms <- dir(path = input_path_hrms)
          file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM_SA[which(PARAM_SA[, 1] == "PARAM0012"), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
        } else {
          samples_string <- PARAM_SA[which(PARAM_SA[, 1] == "PARAM0011"), 2]
          file_name_hrms <- strsplit(samples_string, ";")[[1]]
        }
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
    ##
    x0014 <- which(PARAM_SA[, 1] == "PARAM0014")
    if (length(x0014) == 0) {
      print("ERROR!!! Problem with PARAM0014!")
      checkpoint_parameter <- FALSE
    } else {
      output_path <- PARAM_SA[x0014, 2]
      output_path <- gsub("\\", "/", output_path, fixed = TRUE)
      PARAM_SA[x0014, 2] <- output_path
      if (!dir.exists(output_path)) {
        print("ERROR!!! Problem with PARAM0014! Please make sure the full path is provided!")
        checkpoint_parameter <- FALSE
      }
    }
    ############################################################################
    xsa0001 <- tryCatch(eval(parse(text = paste0("c(", PARAM_SA[which(PARAM_SA[, 1] == "SA0001"), 2], ")"))), error = function(e){NULL})
    if (is.null(xsa0001)) {
      print("ERROR!!! Problem with SA0001!")
      checkpoint_parameter <- FALSE
    }
    xsa0002 <- tryCatch(eval(parse(text = paste0("c(", PARAM_SA[which(PARAM_SA[, 1] == "SA0002"), 2], ")"))), error = function(e){NULL})
    if (is.null(xsa0002)) {
      print("ERROR!!! Problem with SA0002!")
      checkpoint_parameter <- FALSE
    }
    if (length(xsa0001) != length(xsa0002)) {
      print("ERROR!!! Problem with SA0001 and  SA0002! These two parameters should be in the same length")
      checkpoint_parameter <- FALSE
    }
    SA0003 <- as.numeric(PARAM_SA[which(PARAM_SA[, 1] == "SA0003"), 2])
    if (length(SA0003) == 0) {
      print("ERROR!!! Problem with SA0003! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (SA0003 < 0) {
        print("ERROR!!! Problem with SA0003! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    xsa0004 <- tryCatch(eval(parse(text = paste0("c(", PARAM_SA[which(PARAM_SA[, 1] == "SA0004"), 2], ")"))), error = function(e){NULL})
    if (is.null(xsa0004)) {
      print("ERROR!!! Problem with SA0004!")
      checkpoint_parameter <- FALSE
    }
    SA0005 <- as.numeric(PARAM_SA[which(PARAM_SA[, 1] == "SA0005"), 2])
    if (length(SA0005) == 0) {
      print("ERROR!!! Problem with SA0005! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (SA0005 < 0) {
        print("ERROR!!! Problem with SA0005! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    SA0006 <- PARAM_SA[which(PARAM_SA[, 1] == "SA0006"), 2]
    if (is.na(SA0006)) {
      print("ERROR!!! Problem with SA0006!")
      checkpoint_parameter <- FALSE
    }
    SA0007 <- tryCatch(eval(parse(text = paste0("c(", PARAM_SA[which(PARAM_SA[, 1] == "SA0007"), 2], ")"))), error = function(e){NULL})
    if (length(SA0007) != 2) {
      print("ERROR!!! Problem with SA0008! This parameter should be a vector of two positive numbers!")
      checkpoint_parameter <- FALSE
    }
    SA0008 <- as.numeric(PARAM_SA[which(PARAM_SA[, 1] == "SA0008"), 2])
    if (length(SA0008) == 0) {
      print("ERROR!!! Problem with SA0008! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (SA0008 < 0) {
        print("ERROR!!! Problem with SA0008! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    SA0009 <- as.numeric(PARAM_SA[which(PARAM_SA[, 1] == "SA0009"), 2])
    if (length(SA0009) == 0) {
      print("ERROR!!! Problem with SA0009! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (SA0009 >= 1) {
        if ((SA0009 %% 1) != 0) {
          print("ERROR!!! Problem with SA0009! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        print("ERROR!!! Problem with SA0009! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    x0010 <- PARAM_SA[which(PARAM_SA[, 1] == "SA0010"), 2]
    if (length(x0010) == 0) {
      print("ERROR!!! Problem with SA0010!")
      checkpoint_parameter <- FALSE
      x0010 <- 0
    } else {
      if (tolower(x0010) == "yes" | tolower(x0010) == "no") {
      } else {
        print("ERROR!!! Problem with SA0010!")
        checkpoint_parameter <- FALSE
      }
    }
    x0011 <- PARAM_SA[which(PARAM_SA[, 1] == "SA0011"), 2]
    if (length(x0011) == 0) {
      print("ERROR!!! Problem with SA0011!")
      checkpoint_parameter <- FALSE
      x0011 <- 0
    } else {
      if (tolower(x0011) == "yes" | tolower(x0011) == "no") {
      } else {
        print("ERROR!!! Problem with SA0011!")
        checkpoint_parameter <- FALSE
      }
    }
    if (tolower(x0010) == "no" & tolower(x0011) == "no") {
      print("ERROR!!! Problem with SA0010 and SA0011! Both can not be 'NO'!")
      checkpoint_parameter <- FALSE
    }
  }
  if (checkpoint_parameter == FALSE) {
    print("Please visit   https://ufa.idsl.me    for instructions!")
    PARAM_SA <- c()
  }
  return(PARAM_SA)
}
