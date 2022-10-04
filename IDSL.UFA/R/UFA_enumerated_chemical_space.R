UFA_enumerated_chemical_space <- function(PARAM_MF) {
  ##
  x_name_IPDB <- which(PARAM_MF$Parameter == "IPDB file name")
  IPDB_file_name <- PARAM_MF$`User input 2`[x_name_IPDB]
  x_address_IPDB <- which(PARAM_MF$Parameter == "IPDB output address")
  output_path <- gsub("\\", "/", PARAM_MF$`User input 2`[x_address_IPDB], fixed = TRUE)
  ##
  ##############################################################################
  ## To create log record for IDSL.UFA
  initiation_time <- Sys.time()
  timeZone <- tryCatch(Sys.timezone(), warning = function(w) {"UTC"}, error = function(e) {"UTC"})
  .GlobalEnv$logUFA <- paste0(output_path, "/logIPDB_", IPDB_file_name, ".txt")
  UFA_logRecorder(paste0(rep("", 100), collapse = "="))
  UFA_logRecorder(paste0("OUTPUT:  ", output_path))
  UFA_logRecorder(paste0(rep("", 100), collapse = "-"))
  UFA_logRecorder("Initiated isotopic profile database (IPDB) production through the enumerating chemical space approach!")
  UFA_logRecorder(paste0(as.character(initiation_time), " ", timeZone))
  UFA_logRecorder("", printMessage = FALSE)
  UFA_logRecorder("", printMessage = FALSE)
  UFA_logRecorder(paste0(PARAM_MF$`Parameter`, "\t", PARAM_MF$`User input 1`, "\t", PARAM_MF$`User input 2`),  printMessage = FALSE)
  UFA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ##############################################################################
  ##
  x_peak_spacing <- which(PARAM_MF$Parameter == "Peak spacing (Da)")
  peak_spacing <- as.numeric(PARAM_MF$`User input 2`[x_peak_spacing])
  ##
  x_IP_mem_usage <- which(PARAM_MF$Parameter == "Isotopic profile calculations memory usage")
  IP_mem_usage <- PARAM_MF$`User input 2`[x_IP_mem_usage]
  UFA_IP_memeory_variables <- eval(parse(text = IP_mem_usage))
  ##
  x_npc <- which(PARAM_MF$Parameter == "Number of parallel threads")
  number_processing_threads <- as.numeric(PARAM_MF$`User input 2`[x_npc])
  ##
  mr_x <- which(PARAM_MF$Parameter == "Mass range (Da)")
  LOWEST_mass <- as.numeric(PARAM_MF$`User input 1`[mr_x])
  HIGHEST_mass <- as.numeric(PARAM_MF$`User input 2`[mr_x])
  #
  c_x <- which(PARAM_MF$Parameter == "Carbon")
  c_xyz1 <- as.numeric(PARAM_MF$`User input 1`[c_x])
  c_xyz2 <- as.numeric(PARAM_MF$`User input 2`[c_x])
  ##
  Ess_MolVecMat_call <- "Ess_MolVecMat_call <- function(c) { # Carbon range
    do.call(rbind, lapply((b_xyz1):(b_xyz2), function(b) { # Boron range
      do.call(rbind, lapply((br_xyz1):(br_xyz2), function(br) { # Bromine range
        do.call(rbind, lapply((cl_1):(cl_2), function(cl) { # Chlorine range
          SUMbrcl <- br + cl
          if ((SUMbrcl >= (sum_br_cl_xyz_1)) & (SUMbrcl <= (sum_br_cl_xyz_2))) { # Optional condition to remove low probable compounds
            do.call(rbind, lapply((k_1):(k_2), function(k) { # Potassium range
              do.call(rbind, lapply((s_1):(s_2), function(s) { # Sulfur range
                do.call(rbind, lapply((se_1):(se_2), function(se) { # Selenium range
                  do.call(rbind, lapply((si_1):(si_2), function(si) { # Silicon range
                    if ((COND1)) {
                      c(c, b, br, cl, k, s, se, si)
                    }
                  }))
                }))
              }))
            }))
          }
        }))
      }))
    }))
  }"
  ##
  b_x <- which(PARAM_MF$Parameter == "Boron")
  b_xyz1 <- as.numeric(PARAM_MF$`User input 1`[b_x])
  b_xyz2 <- as.numeric(PARAM_MF$`User input 2`[b_x])
  Ess_MolVecMat_call <- gsub("b_xyz1", b_xyz1, Ess_MolVecMat_call)
  Ess_MolVecMat_call <- gsub("b_xyz2", b_xyz2, Ess_MolVecMat_call)
  #
  br_x <- which(PARAM_MF$Parameter == "Bromine")
  br_xyz1 <- as.numeric(PARAM_MF$`User input 1`[br_x])
  br_xyz2 <- as.numeric(PARAM_MF$`User input 2`[br_x])
  Ess_MolVecMat_call <- gsub("br_xyz1", br_xyz1, Ess_MolVecMat_call)
  Ess_MolVecMat_call <- gsub("br_xyz2", br_xyz2, Ess_MolVecMat_call)
  #
  cl_x <- which(PARAM_MF$Parameter == "Chlorine")
  cl_1 <- as.numeric(PARAM_MF$`User input 1`[cl_x])
  cl_2 <- as.numeric(PARAM_MF$`User input 2`[cl_x])
  Ess_MolVecMat_call <- gsub("cl_1", cl_1, Ess_MolVecMat_call)
  Ess_MolVecMat_call <- gsub("cl_2", cl_2, Ess_MolVecMat_call)
  #
  sum_Rule3_x <- grep("Rule 3", PARAM_MF$Parameter, ignore.case = TRUE)
  sum_br_cl_xyz_x <- sum_Rule3_x[1]
  sum_br_cl_xyz_1 <- as.numeric(PARAM_MF$`User input 1`[sum_br_cl_xyz_x])
  sum_br_cl_xyz_2 <- as.numeric(PARAM_MF$`User input 2`[sum_br_cl_xyz_x])
  Ess_MolVecMat_call <- gsub("sum_br_cl_xyz_1", sum_br_cl_xyz_1, Ess_MolVecMat_call)
  Ess_MolVecMat_call <- gsub("sum_br_cl_xyz_2", sum_br_cl_xyz_2, Ess_MolVecMat_call)
  #
  k_x <- which(PARAM_MF$Parameter == "Potassium")
  k_1 <- PARAM_MF$`User input 1`[k_x]
  k_2 <- PARAM_MF$`User input 2`[k_x]
  Ess_MolVecMat_call <- gsub("k_1", k_1, Ess_MolVecMat_call)
  Ess_MolVecMat_call <- gsub("k_2", k_2, Ess_MolVecMat_call)
  #
  s_x <- which(PARAM_MF$Parameter == "Sulfur")
  s_1 <- PARAM_MF$`User input 1`[s_x]
  s_2 <- PARAM_MF$`User input 2`[s_x]
  Ess_MolVecMat_call <- gsub("s_1", s_1, Ess_MolVecMat_call)
  Ess_MolVecMat_call <- gsub("s_2", s_2, Ess_MolVecMat_call)
  #
  se_y <- which(PARAM_MF$Parameter == "Selenium")
  se_1 <- PARAM_MF$`User input 1`[se_y]
  se_2 <- PARAM_MF$`User input 2`[se_y]
  Ess_MolVecMat_call <- gsub("se_1", se_1, Ess_MolVecMat_call)
  Ess_MolVecMat_call <- gsub("se_2", se_2, Ess_MolVecMat_call)
  #
  si_x <- which(PARAM_MF$Parameter == "Silicon")
  si_1 <- PARAM_MF$`User input 1`[si_x]
  si_2 <- PARAM_MF$`User input 2`[si_x]
  Ess_MolVecMat_call <- gsub("si_1", si_1, Ess_MolVecMat_call)
  Ess_MolVecMat_call <- gsub("si_2", si_2, Ess_MolVecMat_call)
  #
  COND1_x <- which(PARAM_MF$Parameter == "Condition1")
  COND1 <- PARAM_MF$`User input 2`[COND1_x]
  if (COND1 == "" | COND1 == " " | COND1 == "1" | tolower(COND1) == "t" | tolower(COND1) == "true") {
    COND1 <- "TRUE"
  }
  Ess_MolVecMat_call <- gsub("COND1", COND1, Ess_MolVecMat_call)
  ##############################################################################
  eval(parse(text = Ess_MolVecMat_call))
  ##############################################################################
  Ess_IP_call <- "Ess_IP_call <- function(counter) {
    c <- Ess_MolVecMat[counter, 1]
    b <- Ess_MolVecMat[counter, 2]
    br <- Ess_MolVecMat[counter, 3]
    cl <- Ess_MolVecMat[counter, 4]
    k <- Ess_MolVecMat[counter, 5]
    s <- Ess_MolVecMat[counter, 6]
    se <- Ess_MolVecMat[counter, 7]
    si <- Ess_MolVecMat[counter, 8]
    Ess_MolVec <- c(c, b, br, cl, k, s, se, si)
    intensity_cutoff <- int_cutoff_str
    isotopic_profile_calculator(Ess_MolVec, Essential_Elements_mass_abundance, peak_spacing, intensity_cutoff, UFA_IP_memeory_variables)
  }"
  ##
  x_int_cutoff_str <- which(PARAM_MF$Parameter == "Theoretical isotopic profile intensity cutoff (%)")
  int_cutoff_str <- PARAM_MF$`User input 2`[x_int_cutoff_str]
  Ess_IP_call <- gsub("int_cutoff_str", int_cutoff_str, Ess_IP_call)
  ##############################################################################
  eval(parse(text = Ess_IP_call))
  ##
  Essential_Elements <- c("C", "B", "Br", "Cl", "K", "S", "Se", "Si")
  L_Essential_Elements <- length(Essential_Elements)
  EL_ESSE <- element_sorter(ElementList = Essential_Elements, ElementOrder = "same")
  Essential_Elements_mass_abundance <- EL_ESSE[["massAbundanceList"]]
  ##############################################################################
  Ess_IPDB_mat_call <- function(counter) {
    IPP <- Ess_IP[[counter]]
    x_100 <- which.max(IPP[, 2])
    L_IPP <- length(IPP[, 2])
    ##
    IP_R13C <- 0
    if (L_IPP > x_100) {
      M13C <- abs(IPP[, 1] - IPP[x_100, 1] - 1.00335484)
      M13C <- M13C[(x_100 + 1):L_IPP]
      x_101 <- which.min(M13C)[1]
      if (M13C[x_101] <= 0.015) {
        x_101 <- x_101 + x_100
        IP_R13C <- IPP[x_101, 2]/IPP[x_100, 2]*100
      }
    }
    c(IPP[L_IPP, 1], IPP[x_100, 1], IP_R13C, x_100, L_IPP)
  }
  ##############################################################################
  MolVecMat_call <- "MolVecMat_call <- function(counter) {
    c <- Ess_MolVecMat[counter, 1]
    b <- Ess_MolVecMat[counter, 2]
    br <- Ess_MolVecMat[counter, 3]
    cl <- Ess_MolVecMat[counter, 4]
    k <- Ess_MolVecMat[counter, 5]
    s <- Ess_MolVecMat[counter, 6]
    se <- Ess_MolVecMat[counter, 7]
    si <- Ess_MolVecMat[counter, 8]
    ##
    do.call(rbind, lapply((n_1):(n_2), function(n) { # Nitrogen range
      do.call(rbind, lapply((h_1):(h_2), function(h) { # Hydrogen range
        do.call(rbind, lapply((as_xyz1):(as_xyz2), function(as) { # Arsenic range
          do.call(rbind, lapply((f_xyz1):(f_xyz2), function(f) { # Fluorine range
            do.call(rbind, lapply((i_xyz1):(i_xyz2), function(i) { # Iodine range
              ##
              SUMbrclfi <- br + cl + f + i
              if ((SUMbrclfi >= (sum_br_cl_f_i_1)) & (SUMbrclfi <= (sum_br_cl_f_i_2))) { # Optional condition to remove low probable compounds
                ##
                if (rule1Check) {
                  SUMhbrclfi <- SUMbrclfi + h
                  rule1 <- ifelse((SUMhbrclfi >= (c/2 - n - 1) & SUMhbrclfi <= (2*c + 3*n + 6)), TRUE, FALSE)
                } else {
                  rule1 <- TRUE
                }
                if (rule1) {
                  ##
                  do.call(rbind, lapply((na_1):(na_2), function(na) { # Sodium range
                    do.call(rbind, lapply((o_1):(o_2), function(o) { # Oxygen range
                      do.call(rbind, lapply((p_1):(p_2), function(p) { # Phosphorus range
                        if ((COND2)) {
                          ##
                          NUM14elements <- length(which(c(c, h, as, b, br, cl, f, i, n, o, p, s, se, si) > 0)) ## Na and K were mot included in this equation
                          if (NUM14elements <= maxNUM14elements) {
                            ##
                            mol_vec <- c(c, h, as, b, br, cl, f, i, k, n, na, o, p, s, se, si)
                            if (extended_SENIOR_rule_str) {
                              rule2 <- extended_SENIOR_rule_check(mol_vec, valence_vec, ionization_correction = ipw_n)
                            } else {
                              rule2 <- TRUE
                            }
                            if (rule2) {
                              c(mol_vec, counter)
                            }
                          }
                        }
                      }))
                    }))
                  }))
                }
              }
            }))
          }))
        }))
      }))
    }))
  }"
  ##
  n_x <- which(PARAM_MF$Parameter == "Nitrogen")
  n_1 <- PARAM_MF$`User input 1`[n_x]
  n_2 <- PARAM_MF$`User input 2`[n_x]
  MolVecMat_call <- gsub("n_1", n_1, MolVecMat_call)
  MolVecMat_call <- gsub("n_2", n_2, MolVecMat_call)
  #
  h_x <- which(PARAM_MF$Parameter == "Hydrogen")
  h_1 <- PARAM_MF$`User input 1`[h_x]
  h_2 <- PARAM_MF$`User input 2`[h_x]
  MolVecMat_call <- gsub("h_1", h_1, MolVecMat_call)
  MolVecMat_call <- gsub("h_2", h_2, MolVecMat_call)
  #
  as_xyz <- which(PARAM_MF$Parameter == "Arsenic")
  as_xyz1 <- PARAM_MF$`User input 1`[as_xyz]
  as_xyz2 <- PARAM_MF$`User input 2`[as_xyz]
  MolVecMat_call <- gsub("as_xyz1", as_xyz1, MolVecMat_call)
  MolVecMat_call <- gsub("as_xyz2", as_xyz2, MolVecMat_call)
  #
  f_xyz <- which(PARAM_MF$Parameter == "Fluorine")
  f_xyz1 <- PARAM_MF$`User input 1`[f_xyz]
  f_xyz2 <- PARAM_MF$`User input 2`[f_xyz]
  MolVecMat_call <- gsub("f_xyz1", f_xyz1, MolVecMat_call)
  MolVecMat_call <- gsub("f_xyz2", f_xyz2, MolVecMat_call)
  #
  i_xyz <- which(PARAM_MF$Parameter == "Iodine")
  i_xyz1 <- PARAM_MF$`User input 1`[i_xyz]
  i_xyz2 <- PARAM_MF$`User input 2`[i_xyz]
  MolVecMat_call <- gsub("i_xyz1", i_xyz1, MolVecMat_call)
  MolVecMat_call <- gsub("i_xyz2", i_xyz2, MolVecMat_call)
  #
  sum_br_cl_f_i_x <- sum_Rule3_x[2]
  sum_br_cl_f_i_1 <- PARAM_MF$`User input 1`[sum_br_cl_f_i_x]
  sum_br_cl_f_i_2 <- PARAM_MF$`User input 2`[sum_br_cl_f_i_x]
  MolVecMat_call <- gsub("sum_br_cl_f_i_1", sum_br_cl_f_i_1, MolVecMat_call)
  MolVecMat_call <- gsub("sum_br_cl_f_i_2", sum_br_cl_f_i_2, MolVecMat_call)
  #
  na_x <- which(PARAM_MF$Parameter == "Sodium")
  na_1 <- PARAM_MF$`User input 1`[na_x]
  na_2 <- PARAM_MF$`User input 2`[na_x]
  MolVecMat_call <- gsub("na_1", na_1, MolVecMat_call)
  MolVecMat_call <- gsub("na_2", na_2, MolVecMat_call)
  #
  o_x <- which(PARAM_MF$Parameter == "Oxygen")
  o_1 <- PARAM_MF$`User input 1`[o_x]
  o_2 <- PARAM_MF$`User input 2`[o_x]
  MolVecMat_call <- gsub("o_1", o_1, MolVecMat_call)
  MolVecMat_call <- gsub("o_2", o_2, MolVecMat_call)
  #
  p_x <- which(PARAM_MF$Parameter == "Phosphorus")
  p_1 <- PARAM_MF$`User input 1`[p_x]
  p_2 <- PARAM_MF$`User input 2`[p_x]
  MolVecMat_call <- gsub("p_1", p_1, MolVecMat_call)
  MolVecMat_call <- gsub("p_2", p_2, MolVecMat_call)
  #
  x_rule1 <- grep("Rule 1", PARAM_MF$Parameter, ignore.case = TRUE)
  str_rule1 <- PARAM_MF$`User input 2`[x_rule1]
  if (str_rule1 == "1" | tolower(str_rule1) == "t" | tolower(str_rule1) == "true") {
    rule1Check <- "TRUE"
  } else {
    rule1Check <- "FALSE"
  }
  MolVecMat_call <- gsub("rule1Check", rule1Check, MolVecMat_call)
  #
  ext_sen_x <- grep("Extended SENIOR rule", PARAM_MF$Parameter, ignore.case = TRUE)
  ext_sen <- PARAM_MF$`User input 2`[ext_sen_x]
  if (ext_sen == "1" | tolower(ext_sen) == "t" | tolower(ext_sen) == "true") {
    extended_SENIOR_rule_str <- "TRUE"
  } else {
    extended_SENIOR_rule_str <- "FALSE"
  }
  MolVecMat_call <- gsub("extended_SENIOR_rule_str", extended_SENIOR_rule_str, MolVecMat_call)
  #
  if (extended_SENIOR_rule_str == "TRUE") {
    ipw_x <- which(PARAM_MF$Parameter == "Ionization pathway")
    ipw <- PARAM_MF$`User input 2`[ipw_x]
    if (ipw == "[M+H/K/Na]" | ipw == "[M+H]" | ipw == "[M+K]" | ipw == "[M+Na]") {
      ipw_n <- "-1"
    } else if (ipw == "[M-H]") {
      ipw_n <- "+1"
    } else {
      ipw_n <- "0"
    }
    MolVecMat_call <- gsub("ipw_n", ipw_n, MolVecMat_call)
  }
  #
  maxNUM14elements <- as.numeric(PARAM_MF$`User input 2`[grep("Rule 4", PARAM_MF$Parameter, ignore.case = TRUE)])
  #
  COND2_x <- which(PARAM_MF$Parameter == "Condition2")
  COND2 <- PARAM_MF$`User input 2`[COND2_x]
  if (COND2 == "" | COND2 == " " | COND2 == "1" | tolower(COND2) == "t" | tolower(COND2) == "true") {
    COND2 <- "TRUE"
  }
  MolVecMat_call <- gsub("COND2", COND2, MolVecMat_call)
  ##############################################################################
  eval(parse(text = MolVecMat_call))
  ##
  ElementsAlphabetical <- c("C", "H", "As", "B", "Br", "Cl", "F", "I", "K", "N", "Na", "O", "P", "S", "Se", "Si")
  L_ElementsAlphabetical <- length(ElementsAlphabetical)
  L_ElementsAlphabetical_1 <- L_ElementsAlphabetical + 1
  if (extended_SENIOR_rule_str == "TRUE") {
    EL_Alphabetical <- element_sorter(ElementList = ElementsAlphabetical, ElementOrder = "same")
    valence_vec <- EL_Alphabetical[["Valence"]]
  }
  ##############################################################################
  NoNEss_mass_call <- function(counter) {
    NoNEss_Molecule <- c(h[counter], as[counter], f[counter], i[counter], n[counter], na[counter], o[counter], p[counter])
    monoisotopic_mass_calculator(NoNEss_Molecule, NonEssential_Elements_mass_abundance)
  }
  ##
  NonEssential_Elements <- c("H", "As", "F", "I", "N", "Na", "O", "P")
  EL_NoNESSE <- element_sorter(ElementList = NonEssential_Elements, ElementOrder = "same")
  NonEssential_Elements_mass_abundance <- EL_NoNESSE[["massAbundanceList"]]
  ##############################################################################
  UFA_logRecorder("Initiated calculating isotopic profiles for essential elements!")
  ##
  if (number_processing_threads == 1) {
    ##
    Ess_MolVecMat <- do.call(rbind, lapply(c_xyz1:c_xyz2, function(c) {
      Ess_MolVecMat_call(c)
    }))
    Ess_MolVecMat <- matrix(Ess_MolVecMat, ncol = L_Essential_Elements)
    L_Ess_MolVecMat <- dim(Ess_MolVecMat)[1]
    ##
    Ess_IP <- lapply(1:L_Ess_MolVecMat, function(counter) {
      Ess_IP_call(counter)
    })
    ##
    Ess_IPDB_mat <- do.call(rbind, lapply(1:L_Ess_MolVecMat, function(counter) {
      Ess_IPDB_mat_call(counter)
    }))
    #
    x_Ess_IP <- which(Ess_IPDB_mat[, 1] <= HIGHEST_mass)
    Ess_MolVecMat <- matrix(Ess_MolVecMat[x_Ess_IP, ], ncol = L_Essential_Elements)
    Ess_IPDB_mat <- matrix(Ess_IPDB_mat[x_Ess_IP, ], ncol = 5)
    Ess_IPDB_mat <- matrix(Ess_IPDB_mat[, -1], ncol = 4)
    Ess_IP <- Ess_IP[x_Ess_IP]
    L_Ess_MolVecMat <- length(x_Ess_IP)
    UFA_logRecorder(paste0("Completed calculating isotopic profiles for essential elements with ", L_Ess_MolVecMat, " combinations!"))
    ##
    UFA_logRecorder("Initiated enumerating molecular formulas!")
    MolVecMat <- do.call(rbind, lapply(1:L_Ess_MolVecMat, function(counter) {
      MolVecMat_call(counter)
    }))
    MolVecMat <- matrix(MolVecMat, ncol = L_ElementsAlphabetical_1)
    L_MolVecMat <- dim(MolVecMat)[1]
    ##
    UFA_logRecorder(paste0("Initiated calculating masses for ", L_MolVecMat, " enumerated molecular formulas!"))
    h <- MolVecMat[, 2]
    as <- MolVecMat[, 3]
    f <- MolVecMat[, 7]
    i <- MolVecMat[, 8]
    n <- MolVecMat[, 10]
    na <- MolVecMat[, 11]
    o <- MolVecMat[, 12]
    p <- MolVecMat[, 13]
    #
    NoNEss_mass <- do.call(rbind, lapply(1:L_MolVecMat, function(counter) {
      NoNEss_mass_call(counter)
    }))
    #
    h <- 0
    as <- 0
    f <- 0
    i <- 0
    n <- 0
    na <- 0
    o <- 0
    p <- 0
    #
    x_mass <- which(NoNEss_mass <= (HIGHEST_mass - 12*c_xyz1)) # which greater than one carbon atom mass
    L_NoNEss_mass <- length(x_mass)
    if (L_NoNEss_mass < L_MolVecMat) {
      NoNEss_mass <- NoNEss_mass[x_mass]
      MolVecMat <- matrix(MolVecMat[x_mass, ], ncol = L_ElementsAlphabetical_1)
      L_MolVecMat <- dim(MolVecMat)[1]
    }
    UFA_logRecorder(paste0("Completed calculating masses for ", L_NoNEss_mass, " molecular formula combinations!"))
    ##
    UFA_logRecorder("Initiated creating the isotopic profile database!")
    x_IP <- c(0, which(abs(diff(MolVecMat[, L_ElementsAlphabetical_1])) > 0), L_MolVecMat)
    #
    L_IP_combination <- length(x_IP) - 1
    ID_IP <- do.call(rbind, lapply(2:(L_IP_combination + 1), function(counter) {
      row_number_first <- (x_IP[counter - 1] + 1)
      c(row_number_first, x_IP[counter], MolVecMat[row_number_first, L_ElementsAlphabetical_1])
    }))
    ##
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      Ess_MolVecMat <- do.call(rbind, mclapply(c_xyz1:c_xyz2, function(c) {
        Ess_MolVecMat_call(c)
      }, mc.cores = number_processing_threads))
      Ess_MolVecMat <- matrix(Ess_MolVecMat, ncol = L_Essential_Elements)
      L_Ess_MolVecMat <- dim(Ess_MolVecMat)[1]
      ##
      Ess_IP <- mclapply(1:L_Ess_MolVecMat, function(counter) {
        Ess_IP_call(counter)
      }, mc.cores = number_processing_threads)
      ##
      Ess_IPDB_mat <- do.call(rbind, mclapply(1:L_Ess_MolVecMat, function(counter) {
        Ess_IPDB_mat_call(counter)
      }, mc.cores = number_processing_threads))
      #
      x_Ess_IP <- which(Ess_IPDB_mat[, 1] <= HIGHEST_mass)
      Ess_MolVecMat <- matrix(Ess_MolVecMat[x_Ess_IP, ], ncol = L_Essential_Elements)
      Ess_IPDB_mat <- matrix(Ess_IPDB_mat[x_Ess_IP, ], ncol = 5)
      Ess_IPDB_mat <- matrix(Ess_IPDB_mat[, -1], ncol = 4)
      Ess_IP <- Ess_IP[x_Ess_IP]
      L_Ess_MolVecMat <- length(x_Ess_IP)
      UFA_logRecorder(paste0("Completed calculating isotopic profiles for essential elements with ", L_Ess_MolVecMat, " combinations!"))
      ##
      UFA_logRecorder("Initiated enumerating molecular formulas!")
      MolVecMat <- do.call(rbind, mclapply(1:L_Ess_MolVecMat, function(counter) {
        MolVecMat_call(counter)
      }, mc.cores = number_processing_threads))
      MolVecMat <- matrix(MolVecMat, ncol = L_ElementsAlphabetical_1)
      L_MolVecMat <- dim(MolVecMat)[1]
      ##
      UFA_logRecorder(paste0("Initiated calculating masses for ", L_MolVecMat, " enumerated molecular formulas!"))
      h <- MolVecMat[, 2]
      as <- MolVecMat[, 3]
      f <- MolVecMat[, 7]
      i <- MolVecMat[, 8]
      n <- MolVecMat[, 10]
      na <- MolVecMat[, 11]
      o <- MolVecMat[, 12]
      p <- MolVecMat[, 13]
      #
      NoNEss_mass <- do.call(rbind, mclapply(1:L_MolVecMat, function(counter) {
        NoNEss_mass_call(counter)
      }, mc.cores = number_processing_threads))
      #
      h <- 0
      as <- 0
      f <- 0
      i <- 0
      n <- 0
      na <- 0
      o <- 0
      p <- 0
      #
      x_mass <- which(NoNEss_mass <= (HIGHEST_mass - 12*c_xyz1)) # which greater than one carbon atom mass
      L_NoNEss_mass <- length(x_mass)
      if (L_NoNEss_mass < L_MolVecMat) {
        NoNEss_mass <- NoNEss_mass[x_mass]
        MolVecMat <- matrix(MolVecMat[x_mass, ], ncol = L_ElementsAlphabetical_1)
        L_MolVecMat <- dim(MolVecMat)[1]
      }
      UFA_logRecorder(paste0("Completed calculating masses for ", L_NoNEss_mass, " molecular formula combinations!"))
      ##
      UFA_logRecorder("Initiated creating the isotopic profile database!")
      x_IP <- c(0, which(abs(diff(MolVecMat[, L_ElementsAlphabetical_1])) > 0), L_MolVecMat)
      #
      L_IP_combination <- length(x_IP) - 1
      ID_IP <- do.call(rbind, mclapply(2:(L_IP_combination + 1), function(counter) {
        row_number_first <- (x_IP[counter - 1] + 1)
        c(row_number_first, x_IP[counter], MolVecMat[row_number_first, L_ElementsAlphabetical_1])
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      Ess_MolVecMat <- foreach(c = c_xyz1:c_xyz2, .combine = 'rbind', .verbose = FALSE) %dopar% {
        Ess_MolVecMat_call(c)
      }
      Ess_MolVecMat <- matrix(Ess_MolVecMat, ncol = L_Essential_Elements)
      L_Ess_MolVecMat <- dim(Ess_MolVecMat)[1]
      ##
      Ess_IP <- foreach(counter = 1:L_Ess_MolVecMat, .verbose = FALSE) %dopar% {
        Ess_IP_call(counter)
      }
      ##
      Ess_IPDB_mat <- foreach(counter = 1:L_Ess_MolVecMat, .combine = 'rbind', .verbose = FALSE) %dopar% {
        Ess_IPDB_mat_call(counter)
      }
      #
      x_Ess_IP <- which(Ess_IPDB_mat[, 1] <= HIGHEST_mass)
      Ess_MolVecMat <- matrix(Ess_MolVecMat[x_Ess_IP, ], ncol = L_Essential_Elements)
      Ess_IPDB_mat <- matrix(Ess_IPDB_mat[x_Ess_IP, ], ncol = 5)
      Ess_IPDB_mat <- matrix(Ess_IPDB_mat[, -1], ncol = 4)
      Ess_IP <- Ess_IP[x_Ess_IP]
      L_Ess_MolVecMat <- length(x_Ess_IP)
      UFA_logRecorder(paste0("Completed calculating isotopic profiles for essential elements with ", L_Ess_MolVecMat, " combinations!"))
      ##
      UFA_logRecorder("Initiated enumerating molecular formulas!")
      MolVecMat <- foreach(counter = 1:L_Ess_MolVecMat, .combine = 'rbind', .verbose = FALSE) %dopar% {
        MolVecMat_call(counter)
      }
      MolVecMat <- matrix(MolVecMat, ncol = L_ElementsAlphabetical_1)
      L_MolVecMat <- dim(MolVecMat)[1]
      ##
      UFA_logRecorder(paste0("Initiated calculating masses for ", L_MolVecMat, " enumerated molecular formulas!"))
      h <- MolVecMat[, 2]
      as <- MolVecMat[, 3]
      f <- MolVecMat[, 7]
      i <- MolVecMat[, 8]
      n <- MolVecMat[, 10]
      na <- MolVecMat[, 11]
      o <- MolVecMat[, 12]
      p <- MolVecMat[, 13]
      #
      NoNEss_mass <- foreach(counter = 1:L_MolVecMat, .combine = 'rbind', .verbose = FALSE) %dopar% {
        NoNEss_mass_call(counter)
      }
      #
      h <- 0
      as <- 0
      f <- 0
      i <- 0
      n <- 0
      na <- 0
      o <- 0
      p <- 0
      #
      x_mass <- which(NoNEss_mass <= (HIGHEST_mass - 12*c_xyz1)) # which greater than one carbon atom mass
      L_NoNEss_mass <- length(x_mass)
      if (L_NoNEss_mass < L_MolVecMat) {
        NoNEss_mass <- NoNEss_mass[x_mass]
        MolVecMat <- matrix(MolVecMat[x_mass, ], ncol = L_ElementsAlphabetical_1)
        L_MolVecMat <- dim(MolVecMat)[1]
      }
      UFA_logRecorder(paste0("Completed calculating masses for ", L_NoNEss_mass, " molecular formula combinations!"))
      ##
      UFA_logRecorder("Initiated creating the isotopic profile database!")
      x_IP <- c(0, which(diff(MolVecMat[, L_ElementsAlphabetical_1]) > 0), L_MolVecMat)
      #
      L_IP_combination <- length(x_IP) - 1
      ID_IP <- foreach(counter = 2:(L_IP_combination + 1), .combine = 'rbind', .verbose = FALSE) %dopar% {
        row_number_first <- (x_IP[counter - 1] + 1)
        c(row_number_first, x_IP[counter], MolVecMat[row_number_first, L_ElementsAlphabetical_1])
      }
      ##
      stopCluster(clust)
      ##
    }
  }
  ##
  ##############################################################################
  ##
  MolVecMat <- matrix(MolVecMat[, -L_ElementsAlphabetical_1], ncol = L_ElementsAlphabetical)
  ##
  progressBARboundaries <- txtProgressBar(min = 0, max = L_IP_combination, initial = 0, style = 3)
  ip_export <- lapply(1:L_IP_combination, function(i) {
    setTxtProgressBar(progressBARboundaries, i)
    counter <- ID_IP[i, 3]
    Counter_IP <- 0
    counter_ip_export <- list()
    IPP <- Ess_IP[[counter]]
    IP_profile <- round(IPP[, 2], 3)
    IP_R13C <- round(Ess_IPDB_mat[counter, 2], 3)
    x_100 <- Ess_IPDB_mat[counter, 3]
    L_IPP <- Ess_IPDB_mat[counter, 4]
    for (j in ID_IP[i, 1]:ID_IP[i, 2]) {
      IP_mass <- IPP[, 1] + NoNEss_mass[j]
      if (IP_mass[L_IPP] <= HIGHEST_mass & IP_mass[1] >= LOWEST_mass) {
        Counter_IP <- Counter_IP + 1
        IP_mass <- round(IP_mass, 6)
        IP <- cbind(IP_mass, IP_profile)
        counter_ip_export[[Counter_IP]] <- list(ipl = j, ipm = IP_mass[x_100], ipp = IP, r13c = IP_R13C, i_100 = x_100, lip = L_IPP)
      }
    }
    counter_ip_export
  })
  close(progressBARboundaries)
  Ess_IP <- NULL
  NoNEss_mass <- 0
  ip_export <- unlist(ip_export, recursive = FALSE)
  L_ip_export <- length(ip_export)
  UFA_logRecorder(paste0("There are ", L_ip_export, " molecular formulas in this isotopic profile database (IPDB) after applying the entire criteria!"))
  ##############################################################################
  gc()
  ##
  xMolVecMat <- do.call(c, lapply(1:L_ip_export, function(x) {
    ip_export[[x]]$ipl
  }))
  MolVecMat <- matrix(MolVecMat[xMolVecMat, ], ncol = L_ElementsAlphabetical)
  #
  x_element_non0 <- do.call(c, lapply(1:L_ElementsAlphabetical, function(x) {
    x_non0 <- which(MolVecMat[, x] > 0)
    if (length(x_non0) > 0) {
      x
    }
  }))
  #
  IP_library <- list(Elements = ElementsAlphabetical[x_element_non0], MolVecMat[, x_element_non0])
  MolVecMat <- NULL
  names(IP_library) <- c("Elements", "MolecularFormulaMatrix")
  ##
  IP_Mass <- do.call(c, lapply(1:L_ip_export, function(x) {
    ip_export[[x]]$ipm
  }))
  ##
  IsotopicProfile <- lapply(1:L_ip_export, function(x) {
    ip_export[[x]]$ipp
  })
  ##
  IP_R13C <- do.call(c, lapply(1:L_ip_export, function(x) {
    ip_export[[x]]$r13c
  }))
  ##
  Index_MAIso <- do.call(c, lapply(1:L_ip_export, function(x) {
    ip_export[[x]]$i_100
  }))
  ##
  IP_size <- do.call(c, lapply(1:L_ip_export, function(x) {
    ip_export[[x]]$lip
  }))
  ##
  ip_export <- NULL
  ##
  IDroundMass <- cbind(round(IP_Mass, digits = 2), seq(1, L_ip_export, 1))
  IDroundMass <- IDroundMass[order(IDroundMass[, 1], decreasing = FALSE), ]
  xDiff <- c(0, which(abs(diff(IDroundMass[, 1])) > 0), L_ip_export)
  LDiff <- length(xDiff) - 1
  AggregatedList <- lapply(1:LDiff, function(i) {
    IDroundMass[(xDiff[i] + 1):xDiff[i + 1], 2]
  })
  names(AggregatedList) <- IDroundMass[(xDiff[1:LDiff] + 1), 1]
  ##
  UFA_logRecorder("Completed generating the isotopic profile database (IPDB)!")
  ##
  PARAM_MF$`User input 2`[x_address_IPDB] <- NA
  IPDB <- list(PARAM_MF, AggregatedList, IP_Mass, IP_library, IsotopicProfile, IP_R13C, Index_MAIso, IP_size)
  names(IPDB) <- c("logIPDB", "AggregatedList", "MassMAIso", "MolecularFormulaDB", "IsotopicProfile", "R13C", "IndexMAIso", "IPsize")
  ##
  UFA_logRecorder("Initiated saving the isotopic profile database!")
  address_IPDB <- paste0(output_path, "/", IPDB_file_name, ".Rdata")
  save(IPDB, file = address_IPDB)
  ##
  ##############################################################################
  ##
  completion_time <- Sys.time()
  UFA_logRecorder(paste0(rep("", 100), collapse = "-"))
  required_time <- completion_time - initiation_time
  print(required_time)
  UFA_logRecorder(paste0(as.character(completion_time), " ", timeZone), printMessage = FALSE)
  UFA_logRecorder("", printMessage = FALSE)
  UFA_logRecorder("", printMessage = FALSE)
  UFA_logRecorder("Stored isotopic profile database (IPDB) from the enumerating chemical space approach!")
  UFA_logRecorder(paste0(rep("", 100), collapse = "="), printMessage = FALSE)
  ##
  ##############################################################################
  ##
  gc()
  closeAllConnections()
  #
  return()
}
