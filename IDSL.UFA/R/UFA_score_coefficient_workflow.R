UFA_score_coefficient_workflow <- function(spreadsheet) {
  print("Initiated testing the `score coefficient spreadsheet` tab consistency!")
  PARAM_SFT <- UFA_score_function_optimization_xlsxAnalyzer(spreadsheet)
  if (!is.null(PARAM_SFT)) {
    print("The `score coefficient spreadsheet` tab is consistent with the score optimization workflow!")
    ##
    ##########################################################################
    ## To create log record for IDSL.UFA
    initiation_time <- Sys.time()
    timeZone <- tryCatch(Sys.timezone(), warning = function(w) {"UTC"}, error = function(e) {"UTC"})
    input_path <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0006'), 2]
    output_path <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0010'), 2]
    .GlobalEnv$logIPA <- paste0(output_path, "/logUFA_score_function_optimization.txt")
    IPA_logRecorder(paste0(rep("", 100), collapse = "="))
    IPA_logRecorder(paste0("mzML/mzXML/netCDF:  ", input_path))
    IPA_logRecorder(paste0("OUTPUT:  ", output_path))
    IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
    IPA_logRecorder("Initiated score function optimization workflow!")
    IPA_logRecorder(paste0(as.character(initiation_time), " ", timeZone))
    IPA_logRecorder("", printMessage = FALSE)
    IPA_logRecorder("", printMessage = FALSE)
    IPA_logRecorder(paste0(PARAM_SFT[, 1], "\t", PARAM_SFT[, 2]),  printMessage = FALSE)
    IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
    ##
    ##########################################################################
    ##
    x0001 <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0001'), 2]
    if (tolower(x0001) == "yes") {
      zero_score_function(PARAM_SFT)
    }
    ##
    x0003 <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0003'), 2]
    if (tolower(x0003) == "yes") {
      score_coefficients_optimization(PARAM_SFT)
    }
    ##
    x0004 <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0004'), 2]
    if (tolower(x0004) == "yes") {
      score_coefficient_evaluation(PARAM_SFT)
    }
    ##
    ##########################################################################
    ##
    completion_time <- Sys.time()
    IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
    required_time <- completion_time - initiation_time
    print(required_time)
    IPA_logRecorder(paste0(as.character(completion_time), " ", timeZone), printMessage = FALSE)
    IPA_logRecorder("", printMessage = FALSE)
    IPA_logRecorder("", printMessage = FALSE)
    IPA_logRecorder("Completed score function optimization workflow!")
    IPA_logRecorder(paste0(rep("", 100), collapse = "="), printMessage = FALSE)
    ##
    ##########################################################################
    ##
  }
}
