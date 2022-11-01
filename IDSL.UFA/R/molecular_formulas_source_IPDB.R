molecular_formulas_source_IPDB <- function(PARAM_SF) {
  ##
  x_csv_file <- which(PARAM_SF[, 1] == "FS0001")
  Molecular_formula_source_file <- gsub("\\", "/", PARAM_SF[x_csv_file, 2], fixed = TRUE)
  x_address_IPDB <- which(PARAM_SF[, 1] == "FS0007")
  output_path <- PARAM_SF[x_address_IPDB, 2]
  output_path <- gsub("\\", "/", output_path, fixed = TRUE)
  x_name_IPDB <- which(PARAM_SF[, 1] == "FS0008")
  IPDB_file_name <- PARAM_SF[x_name_IPDB, 2]
  ##
  ##############################################################################
  ## To create log record for IDSL.UFA
  initiation_time <- Sys.time()
  timeZone <- tryCatch(Sys.timezone(), warning = function(w) {"UTC"}, error = function(e) {"UTC"})
  .GlobalEnv$logIPA <- paste0(output_path, "/logIPDB_", IPDB_file_name, ".txt")
  IPA_logRecorder(paste0(rep("", 100), collapse = "="))
  IPA_logRecorder(paste0("csv/txt:  ", Molecular_formula_source_file))
  IPA_logRecorder(paste0("OUTPUT:  ", output_path))
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Initiated isotopic profile database (IPDB) production from a source of known molecular formulas!")
  IPA_logRecorder(paste0(as.character(initiation_time), " ", timeZone))
  IPA_logRecorder("", printMessage = FALSE)
  IPA_logRecorder("", printMessage = FALSE)
  IPA_logRecorder(paste0(PARAM_SF[, 1], "\t", PARAM_SF[, 2]),  printMessage = FALSE)
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ##############################################################################
  ##
  strMoleFormFileLocation <- strsplit(Molecular_formula_source_file, "[.]")[[1]]
  moleFormFileFormat <- tolower(strMoleFormFileLocation[length(strMoleFormFileLocation)])
  ##
  if (moleFormFileFormat == "csv") {
    molecular_formula <- data.frame(V1 = as.vector(read.csv(Molecular_formula_source_file, header = FALSE)))
    molecular_formula <- molecular_formula[, 1]
  } else {
    molecular_formula <- readLines(Molecular_formula_source_file, warn = FALSE)
  }
  ##
  molecular_formula <-  gsub(" ", "", molecular_formula)
  molecular_formula <-  gsub("[+]", "", molecular_formula)
  molecular_formula <-  gsub("-", "", molecular_formula)
  ##
  peak_spacing <- as.numeric(PARAM_SF[which(PARAM_SF[, 1] == "FS0002"), 2])
  intensity_cutoff_str <- PARAM_SF[which(PARAM_SF[, 1] == "FS0003"), 2]
  UFA_IP_memeory_variables <- eval(parse(text = paste0("c(", PARAM_SF[which(PARAM_SF[, 1] == "FS0004"), 2], ")")))
  IonPathways <- eval(parse(text = paste0("c(", PARAM_SF[which(PARAM_SF[, 1] == "FS0005"), 2], ")")))
  number_processing_threads <- as.numeric(PARAM_SF[which(PARAM_SF[, 1] == "FS0006"), 2])
  address_IPDB <- paste0(output_path, "/", IPDB_file_name, ".Rdata")
  ##
  IPDB <- isotopic_profile_molecular_formula_feeder(molecular_formula, peak_spacing, intensity_cutoff_str, UFA_IP_memeory_variables, IonPathways, number_processing_threads)
  PARAM_SF[x_csv_file, 2] <- NA
  PARAM_SF[x_address_IPDB, 2] <- NA
  IPDB <- c(list(PARAM_SF), IPDB)
  names(IPDB) <- c("logIPDB", "AggregatedList", "MassMAIso", "MolecularFormulaDB", "IsotopicProfile", "R13C", "IndexMAIso", "IPsize")
  ##
  IPA_logRecorder("Initiated saving the isotopic profile database (IPDB)!")
  save(IPDB, file = address_IPDB)
  ##
  ##############################################################################
  ##
  completion_time <- Sys.time()
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  required_time <- completion_time - initiation_time
  print(required_time)
  IPA_logRecorder(paste0(as.character(completion_time), " ", timeZone), printMessage = FALSE)
  IPA_logRecorder("", printMessage = FALSE)
  IPA_logRecorder("", printMessage = FALSE)
  IPA_logRecorder("Stored isotopic profile database (IPDB) from the source of known molecular formulas!")
  IPA_logRecorder(paste0(rep("", 100), collapse = "="), printMessage = FALSE)
  ##
  ##############################################################################
  ##
  gc()
  closeAllConnections()
  #
  return()
}
