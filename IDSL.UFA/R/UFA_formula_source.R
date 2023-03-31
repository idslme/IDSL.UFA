UFA_formula_source <- function(PARAM_FormSource) {
  ##
  xCSVfile <- which(PARAM_FormSource[, 1] == "FS0001")
  formula_source_file <- gsub("\\", "/", PARAM_FormSource[xCSVfile, 2], fixed = TRUE)
  xAddressIPDB <- which(PARAM_FormSource[, 1] == "FS0002")
  output_path <- PARAM_FormSource[xAddressIPDB, 2]
  output_path <- gsub("\\", "/", output_path, fixed = TRUE)
  xNameIPDB <- which(PARAM_FormSource[, 1] == "FS0003")
  IPDB_file_name <- PARAM_FormSource[xNameIPDB, 2]
  ##
  ##############################################################################
  ## To create log record for IDSL.UFA
  initiation_time <- Sys.time()
  timeZone <- tryCatch(Sys.timezone(), warning = function(w) {"UTC"}, error = function(e) {"UTC"})
  .logIPA <- NULL
  .logIPA <<- paste0(output_path, "/logIPDB_", IPDB_file_name, ".txt")
  IPA_logRecorder(paste0(rep("", 100), collapse = "="))
  IPA_logRecorder("Type <<< citation('IDSL.UFA') >>> for citing this R package in publications.")
  IPA_logRecorder(paste0("csv/txt:  ", formula_source_file))
  IPA_logRecorder(paste0("OUTPUT:  ", output_path))
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Initiated isotopic profile database (IPDB) production from a source of known molecular formulas!")
  IPA_logRecorder(paste0(as.character(initiation_time), " ", timeZone))
  IPA_logRecorder("", allowedPrinting = FALSE)
  IPA_logRecorder("", allowedPrinting = FALSE)
  IPA_logRecorder(paste0(PARAM_FormSource[, 1], "\t", PARAM_FormSource[, 2]),  allowedPrinting = FALSE)
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ##############################################################################
  ##
  strMoleFormFileLocation <- strsplit(formula_source_file, "[.]")[[1]]
  moleFormFileFormat <- tolower(strMoleFormFileLocation[length(strMoleFormFileLocation)])
  ##
  if (moleFormFileFormat == "csv") {
    formulaSourceFile <- data.frame(read.csv(formula_source_file, header = FALSE))
    molecularFormulaDatabase <- as.vector(formulaSourceFile[, 1])
  } else if (moleFormFileFormat == "xlsx") {
    formulaSourceFile <- data.frame(readxl::read_xlsx(formula_source_file, col_names = FALSE))
    molecularFormulaDatabase <- as.vector(formulaSourceFile[, 1])
  } else if (moleFormFileFormat == "txt") {
    molecularFormulaDatabase <- readLines(formula_source_file, warn = FALSE)
  } else {
    stop(IPA_logRecorder("Inconsistent format for 'FS0001'!"))
  }
  ##
  retentionTime <- NULL
  if ((moleFormFileFormat == "csv") | (moleFormFileFormat == "xlsx")) {
    if (dim(formulaSourceFile)[2] >= 2) {
      retentionTime <- as.vector(formulaSourceFile[, 2])
      ##
      nMolecularFormulaDatabase <- length(molecularFormulaDatabase)
      if ((nMolecularFormulaDatabase != length(retentionTime)) | (nMolecularFormulaDatabase == 0)) {
        stop(IPA_logRecorder(paste0("The first and second columns in `", formula_source_file,"` are not in the same size!")))
      }
    }
  }
  formulaSourceFile <- NULL
  ##
  molecularFormulaDatabase <-  gsub(" ", "", molecularFormulaDatabase)
  molecularFormulaDatabase <-  gsub("[+]", "", molecularFormulaDatabase)
  molecularFormulaDatabase <-  gsub("-", "", molecularFormulaDatabase)
  ##
  number_processing_threads <- as.numeric(PARAM_FormSource[which(PARAM_FormSource[, 1] == "FS0004"), 2])
  allowedMustRunCalculation <- eval(parse(text = PARAM_FormSource[which(PARAM_FormSource[, 1] == "FS0005"), 2]))
  IonPathways <- eval(parse(text = paste0("c(", PARAM_FormSource[which(PARAM_FormSource[, 1] == "FS0006"), 2], ")")))
  intensity_cutoff_str <- PARAM_FormSource[which(PARAM_FormSource[, 1] == "FS0007"), 2]
  peak_spacing <- as.numeric(PARAM_FormSource[which(PARAM_FormSource[, 1] == "FS0008"), 2])
  UFA_IP_memeory_variables <- eval(parse(text = paste0("c(", PARAM_FormSource[which(PARAM_FormSource[, 1] == "FS0009"), 2], ")")))
  ##
  IPDB <- molecularFormula2IPdb(molecularFormulaDatabase, retentionTime, peak_spacing, intensity_cutoff_str, IonPathways,
  				number_processing_threads, UFA_IP_memeory_variables, allowedMustRunCalculation, allowedVerbose = TRUE)
  ##
  strSourceFile <- strsplit(formula_source_file, "/")[[1]]
  source_file <- strSourceFile[length(strSourceFile)]
  ##
  PARAM_FormSource[xCSVfile, 2] <- source_file
  PARAM_FormSource[xAddressIPDB, 2] <- NA
  primaryNamesIPDB <- names(IPDB)
  IPDB <- c(list(PARAM_FormSource), IPDB)
  names(IPDB) <- c("logIPDB", primaryNamesIPDB)
  ##
  addressIPDB <- paste0(output_path, "/", IPDB_file_name, ".Rdata")
  IPA_logRecorder("Initiated saving the isotopic profile database (IPDB)!")
  save(IPDB, file = addressIPDB)
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
  IPA_logRecorder("Successfully stored isotopic profile database (IPDB) from known molecular formulas!")
  IPA_logRecorder(paste0(rep("", 100), collapse = "="), allowedPrinting = FALSE)
  ##
  ##############################################################################
  ##
  gc()
  closeAllConnections()
  ##
  return()
}
