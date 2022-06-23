molecular_formulas_source_IPDB <- function(PARAM_SF) {
  print("Initiated isotopic profile database (IPDB) production from a source of known molecular formulas!")
  ##
  x_csv_file <- which(PARAM_SF[, 1] == "FS0001")
  Molecular_formula_source_file <- gsub("\\", "/", PARAM_SF[x_csv_file, 2], fixed = TRUE)
  ##
  molecular_formula <- data.frame(V1 = as.vector(read.csv(Molecular_formula_source_file, header = FALSE)))
  molecular_formula <-  gsub(" ", "", molecular_formula[, 1])
  molecular_formula <-  gsub("[+]", "", molecular_formula)
  molecular_formula <-  gsub("-", "", molecular_formula)
  ##
  peak_spacing <- as.numeric(PARAM_SF[which(PARAM_SF[, 1] == "FS0002"), 2])
  intensity_cutoff_str <- PARAM_SF[which(PARAM_SF[, 1] == "FS0003"), 2]
  UFA_IP_memeory_variables <- eval(parse(text = paste0("c(", PARAM_SF[which(PARAM_SF[, 1] == "FS0004"), 2], ")")))
  IonPathways <- eval(parse(text = paste0("c(", PARAM_SF[which(PARAM_SF[, 1] == "FS0005"), 2], ")")))
  number_processing_threads <- as.numeric(PARAM_SF[which(PARAM_SF[, 1] == "FS0006"), 2])
  x_address_IPDB <- which(PARAM_SF[, 1] == "FS0007")
  x_name_IPDB <- which(PARAM_SF[, 1] == "FS0008")
  address_IPDB <- paste0(PARAM_SF[x_address_IPDB, 2], "/", PARAM_SF[x_name_IPDB, 2], ".Rdata")
  address_IPDB <- gsub("\\", "/", address_IPDB, fixed = TRUE)
  ##
  IPDB <- isotopic_profile_molecular_formula_feeder(molecular_formula, peak_spacing, intensity_cutoff_str, UFA_IP_memeory_variables, IonPathways, number_processing_threads)
  PARAM_SF[x_csv_file, 2] <- NA
  PARAM_SF[x_address_IPDB, 2] <- NA
  IPDB <- c(IPDB, list(PARAM_SF))
  names(IPDB) <- c("MassMAIso", "MolecularFormulaDB", "IsotopicProfile", "R13C", "IndexMAIso", "IPsize", "logIPDB")
  ##
  print("Initiated saving the isotopic profile database")
  save(IPDB, file = address_IPDB)
  #
  gc()
  closeAllConnections()
  #
  print("Saved isotopic profile database (IPDB) from the source of known molecular formulas!")
}
