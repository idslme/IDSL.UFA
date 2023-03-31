UFA_score_function_optimization <- function(PARAM_ScoreFunc) {
  ##
  ##############################################################################
  ## To create log record for IDSL.UFA
  initiation_time <- Sys.time()
  timeZone <- tryCatch(Sys.timezone(), warning = function(w) {"UTC"}, error = function(e) {"UTC"})
  input_path <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0005'), 2]
  excelfile_address <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0006'), 2]
  output_path <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0008'), 2]
  .logIPA <- NULL
  .logIPA <<- paste0(output_path, "/logUFA_score_function_optimization.txt")
  IPA_logRecorder(paste0(rep("", 100), collapse = "="))
  IPA_logRecorder("Type <<< citation('IDSL.UFA') >>> for citing this R package in publications.")
  IPA_logRecorder(paste0("mzML/mzXML/netCDF:  ", input_path))
  IPA_logRecorder(paste0("Reference spreadsheet (.xlsx)", excelfile_address))
  IPA_logRecorder(paste0("OUTPUT:  ", output_path))
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Initiated score function optimization workflow!")
  IPA_logRecorder(paste0(as.character(initiation_time), " ", timeZone))
  IPA_logRecorder("", allowedPrinting = FALSE)
  IPA_logRecorder("", allowedPrinting = FALSE)
  IPA_logRecorder(paste0(PARAM_ScoreFunc[, 1], "\t", PARAM_ScoreFunc[, 2]),  allowedPrinting = FALSE)
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ##############################################################################
  ##
  x0001 <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0001'), 2]
  if (tolower(x0001) == "yes") {
    scoreCoefficientsReplicate(PARAM_ScoreFunc)
  }
  ##
  x0003 <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0003'), 2]
  if (tolower(x0003) == "yes") {
    scoreCoefficientsOptimization(PARAM_ScoreFunc)
  }
  ##
  x0004 <- PARAM_ScoreFunc[which(PARAM_ScoreFunc[, 1] == 'SFT0004'), 2]
  if (tolower(x0004) == "yes") {
    scoreCoefficientsEvaluation(PARAM_ScoreFunc)
  }
  ##
  ##############################################################################
  ##
  completion_time <- Sys.time()
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  required_time <- completion_time - initiation_time
  IPA_logRecorder(paste0("The required processing time was `", required_time, " ", attributes(required_time)$units, "`"))
  IPA_logRecorder(paste0(as.character(completion_time), " ", timeZone), allowedPrinting = FALSE)
  IPA_logRecorder("", allowedPrinting = FALSE)
  IPA_logRecorder("", allowedPrinting = FALSE)
  IPA_logRecorder("Completed score coefficients optimization workflow!")
  IPA_logRecorder(paste0(rep("", 100), collapse = "="), allowedPrinting = FALSE)
  ##
  ##############################################################################
  ##
  return()
}
