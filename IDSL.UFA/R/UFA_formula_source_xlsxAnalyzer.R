UFA_formula_source_xlsxAnalyzer <- function(spreadsheet) {
  ##
  IPA_message("Initiated analyzing the `formula_source` tab!", failedMessage= FALSE)
  ##
  checkpoint_parameter <- FALSE
  ##
  if (typeof(spreadsheet) == "list") {
    if (ncol(spreadsheet) >= 4) {
      PARAM_FormSource <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
      ##
    } else if (ncol(spreadsheet) == 2) {
      PARAM_FormSource <- spreadsheet
      checkpoint_parameter <- TRUE
      ##
    } else {
      IPA_message("The `formula_source` spreadsheet tab was not produced properly!")
    }
  } else if (typeof(spreadsheet) == "character") {
    if (length(spreadsheet) == 1) {
      if (file.exists(spreadsheet)) {
        PARAM_FormSource <- readxl::read_xlsx(spreadsheet, sheet = "formula_source")
        PARAM_FormSource <- cbind(PARAM_FormSource[, 2], PARAM_FormSource[, 4])
        checkpoint_parameter <- TRUE
      } else {
        IPA_message("The `formula_source` spreadsheet tab not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      IPA_message("The `formula_source` spreadsheet tab was not produced properly!")
    }
  } else {
    IPA_message("The `formula_source` spreadsheet tab was not produced properly!")
  }
  ##############################################################################
  if (checkpoint_parameter) {
    ##
    x0001 <- which(PARAM_FormSource[, 1] == 'FS0001')
    if (length(x0001) == 0) {
      IPA_message("ERROR!!! Problem with FS0001! Molecular formula source file is not detected!")
      checkpoint_parameter <- FALSE
    } else {
      Molecular_formula_source_file <- gsub("\\", "/", PARAM_FormSource[x0001, 2], fixed = TRUE)
      PARAM_FormSource[x0001, 2] <- Molecular_formula_source_file
      if (!file.exists(Molecular_formula_source_file)) {
        IPA_message("ERROR!!! Problem with FS0001! Molecular formula source file is not detected!")
        checkpoint_parameter <- FALSE
      } else {
        if (!grepl(".xlsx$|.csv$|.txt$", Molecular_formula_source_file, ignore.case = TRUE)) {
          IPA_message("ERROR!!! Problem with FS0001! Inconsistent format for the formula source!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    x_address <- which(PARAM_FormSource[, 1] == "FS0002")
    formula_source_address <- gsub("\\", "/", PARAM_FormSource[x_address, 2], fixed = TRUE)
    PARAM_FormSource[x_address, 2] <- formula_source_address
    if (!dir.exists(formula_source_address)) {
      dir.create(formula_source_address, recursive = TRUE)
      if (!dir.exists(formula_source_address)) {
        IPA_message(paste0("ERROR!!! Problem with FS0002! Can't create `", formula_source_address, "` folder!"))
        checkpoint_parameter <- FALSE
      }
    }
    ##
    number_processing_threads <- as.numeric(PARAM_FormSource[which(PARAM_FormSource[, 1] == 'FS0004'), 2])
    if (length(number_processing_threads) == 0) {
      IPA_message("ERROR!!! Problem with FS0004! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (number_processing_threads >= 1) {
        if ((number_processing_threads %% 1) != 0) {
          IPA_message("ERROR!!! Problem with FS0004! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        IPA_message("ERROR!!! Problem with FS0004! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0005 <- which(PARAM_FormSource[, 1] == 'FS0005')
    allowedMustRunCalculation <- tolower(gsub(" ", "", PARAM_FormSource[x0005, 2]))
    if (allowedMustRunCalculation == "1" | allowedMustRunCalculation == "t" | allowedMustRunCalculation == "true") {
      allowedMustRunCalculation <- TRUE
    } else {
      allowedMustRunCalculation <- FALSE
    }
    PARAM_FormSource[x0005, 2] <- allowedMustRunCalculation
    ##
    IonPathways <- tryCatch(eval(parse(text = paste0("c(", PARAM_FormSource[which(PARAM_FormSource[, 1] == 'FS0006'), 2], ")"))), error = function(e){NULL})
    if (is.null(IonPathways)) {
      IPA_message("ERROR!!! Problem with FS0006!")
      checkpoint_parameter <- FALSE
    }
    ##
    intensity_cutoff_str <- PARAM_FormSource[which(PARAM_FormSource[, 1] == 'FS0007'), 2]
    if (is.na(intensity_cutoff_str)) {
      IPA_message("ERROR!!! Problem with FS0007!")
      checkpoint_parameter <- FALSE
    } else {
      ##
      c <- 5
      b <- 5
      br <- 5
      cl <- 5
      k <- 5
      s <- 5
      se <- 5
      si <- 5
      ##
      checkStrInt <- FALSE
      tryCatch(eval(parse(text = intensity_cutoff_str)), error = function(e) {checkStrInt <- TRUE})
      if (checkStrInt) {
        IPA_message("ERROR!!! Problem with FS0007!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    peak_spacing <- as.numeric(PARAM_FormSource[which(PARAM_FormSource[, 1] == 'FS0008'), 2])
    if (is.na(peak_spacing)) {
      IPA_message("ERROR!!! Problem with FS0008!")
      checkpoint_parameter <- FALSE
    } else {
      if (peak_spacing < 0) {
        IPA_message("ERROR!!! Problem with FS0008!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    UFA_IP_memeory_variables <- tryCatch(eval(parse(text = paste0("c(", PARAM_FormSource[which(PARAM_FormSource[, 1] == 'FS0009'), 2], ")"))), error = function(e){NULL})
    if (length(UFA_IP_memeory_variables) != 3) {
      IPA_message("ERROR!!! Problem with FS0009! This parameter should be a vector of three positive numbers!")
      checkpoint_parameter <- FALSE
    }
  }
  ##############################################################################
  if (!checkpoint_parameter) {
    PARAM_FormSource <- NULL
  } else {
    IPA_message("Completed analyzing the `formula_source` tab!", failedMessage= FALSE)
  }
  ##
  return(PARAM_FormSource)
}
