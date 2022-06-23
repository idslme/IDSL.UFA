UFA_enumerated_chemical_space_xlsxAnalyzer <- function(PARAM_MF) {
  tr1_sub_pattern <- Sys.time()
  ##
  checkpoint_parameter <- 1
  ##
  x_address_IPDB <- which(PARAM_MF$Parameter == "IPDB output address")
  address_IPDB <- PARAM_MF$`User input 2`[x_address_IPDB]
  address_IPDB <- gsub("\\", "/", address_IPDB, fixed = TRUE)
  if (!dir.exists(address_IPDB)) {
    checkpoint_parameter <- 0
    print("ERROR!!! Please make sure the full path is provided for the 'IPDB output address' in the 'enumerated_chemical_space' tab!")
  }
  ##
  x_peak_spacing <- which(PARAM_MF$Parameter == "Peak spacing (Da)")
  peak_spacing <- as.numeric(PARAM_MF$`User input 2`[x_peak_spacing])
  if (is.na(peak_spacing)) {
    checkpoint_parameter <- 0
    print("Error!!! Peak Spacing in the 'enumerated_chemical_space' tab should be a positive number!")
  }
  ##
  x_IP_mem_usage <- which(PARAM_MF$Parameter == "Isotopic profile calculations memory usage")
  IP_mem_usage <- PARAM_MF$`User input 2`[x_IP_mem_usage]
  UFA_IP_memeory_variables <- eval(parse(text = IP_mem_usage))
  if (length(UFA_IP_memeory_variables) != 2) {
    checkpoint_parameter <- 0
    print("Error!!! Isotopic profile calculations memory usage in the 'enumerated_chemical_space' tab should be a vector of two positive numbers!")
  }
  ##
  x_npc <- which(PARAM_MF$Parameter == "Number of parallel threads")
  number_processing_threads <- as.numeric(PARAM_MF$`User input 2`[x_npc])
  if (is.na(number_processing_threads)) {
    checkpoint_parameter <- 0
    print("Error!!! Number of parallel threads in the 'enumerated_chemical_space' tab should be an integer!")
  }
  ##
  mr_x <- which(PARAM_MF$Parameter == "Mass range (Da)")
  LOWEST_mass <- as.numeric(PARAM_MF$`User input 1`[mr_x])
  HIGHEST_mass <- as.numeric(PARAM_MF$`User input 2`[mr_x])
  if (is.na(LOWEST_mass) | is.na(HIGHEST_mass)) {
    checkpoint_parameter <- 0
    print("Error!!! Check 'Mass range' in the 'enumerated_chemical_space' tab!")
  } else {
    if (!is.na(LOWEST_mass) & !is.na(HIGHEST_mass)) {
      if (LOWEST_mass > HIGHEST_mass) {
        checkpoint_parameter <- 0
        print("Error!!! Check 'Mass range' in the 'enumerated_chemical_space' tab!")
      }
    }
  }
  #
  c_x <- which(PARAM_MF$Parameter == "Carbon")
  c_xyz1 <- as.numeric(PARAM_MF$`User input 1`[c_x])
  c_xyz2 <- as.numeric(PARAM_MF$`User input 2`[c_x])
  if (is.na(c_xyz1) | is.na(c_xyz2)) {
    checkpoint_parameter <- 0
    print("Error!!! Check the 'Carbon' range in the 'enumerated_chemical_space' tab!")
  } else {
    if (!is.na(c_xyz1) & !is.na(c_xyz2)) {
      if (c_xyz1 > c_xyz2) {
        checkpoint_parameter <- 0
        print("Error!!! Check the 'Carbon' range in the 'enumerated_chemical_space' tab!")
      }
    }
  }
  ##
  Ess_MolVecMat_call <- "Ess_MolVecMat_call <- function(c) { # Carbon range
    do.call(rbind, lapply((b_xyz1):(b_xyz2), function(b) { # Boron range
      do.call(rbind, lapply((br_xyz1):(br_xyz2), function(br) { # Bromine range
        do.call(rbind, lapply((cl_1):(cl_2), function(cl) { # Chlorine range
          if (((br + cl) >= (sum_br_cl_xyz_1)) & ((br + cl) <= (sum_br_cl_xyz_2))) { # Optional condition to remove low probable compounds
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
  sum_br_cl_xyz_x <- which(PARAM_MF$Parameter == "Sum Br + Cl")
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
  Essential_Elements <- c("C", "B", "Br", "Cl", "K", "S", "Se", "Si")
  L_Essential_Elements <- length(Essential_Elements)
  ##
  eval(parse(text = Ess_MolVecMat_call))
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
    do.call(rbind, lapply((n_1):(n_2), function(n) { # Nitrogen range
      do.call(rbind, lapply((h_1):(h_2), function(h) { # Hydrogen range
        do.call(rbind, lapply((as_xyz1):(as_xyz2), function(as) { # Arsenic range
          do.call(rbind, lapply((f_xyz1):(f_xyz2), function(f) { # Fluorine range
            do.call(rbind, lapply((i_xyz1):(i_xyz2), function(i) { # Iodine range
              if (((br + cl + f + i) >= (sum_br_cl_f_i_1)) & ((br + cl + f + i) <= (sum_br_cl_f_i_2))) { # Optional condition to remove low probable compounds
                do.call(rbind, lapply((na_1):(na_2), function(na) { # Sodium range
                  do.call(rbind, lapply((o_1):(o_2), function(o) { # Oxygen range
                    do.call(rbind, lapply((p_1):(p_2), function(p) { # Phosphorus range
                      if ((COND2)) {
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
                    }))
                  }))
                }))
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
  sum_br_cl_f_i_x <- which(PARAM_MF$Parameter == "Sum Br + Cl + F + I")
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
  COND2_x <- which(PARAM_MF$Parameter == "Condition2")
  COND2 <- PARAM_MF$`User input 2`[COND2_x]
  if (COND2 == "" | COND2 == " " | COND2 == "1" | tolower(COND2) == "t" | tolower(COND2) == "true") {
    COND2 <- "TRUE"
  }
  MolVecMat_call <- gsub("COND2", COND2, MolVecMat_call)
  #
  ext_sen_x <- which(PARAM_MF$Parameter == "Extended SENIOR rule")
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
    if (ipw == "[M+H/K/Na]" | ipw == "[M+H]") {
      ipw_n <- "-1"
    } else if (ipw == "[M-H]") {
      ipw_n <- "+1"
    } else if (ipw == "[M]") {
      ipw_n <- "0"
    } else {
      checkpoint_parameter <- 0
      print("Error!!! Check 'Ionization pathway' in the 'enumerated_chemical_space' tab!")
      ipw_n <- "0"
    }
    MolVecMat_call <- gsub("ipw_n", ipw_n, MolVecMat_call)
  }
  ##############################################################################
  eval(parse(text = MolVecMat_call))
  ##
  ElementsAlphabetical <- c("C", "H", "As", "B", "Br", "Cl", "F", "I", "K", "N", "Na", "O", "P", "S", "Se", "Si")
  L_ElementsAlphabetical <- length(ElementsAlphabetical)
  if (extended_SENIOR_rule_str == "TRUE") {
    EL_Alphabetical <- element_sorter(ElementList = ElementsAlphabetical, ElementOrder = "same")
    valence_vec <- EL_Alphabetical[[3]]
  }
  ##############################################################################
  osType <- Sys.info()[['sysname']]
  if (osType == "Windows") {
    clust <- makeCluster(number_processing_threads)
    registerDoParallel(clust)
    ##
    print("Initiated counting essential elements combinations!")
    Ess_MolVecMat <- foreach(c = c_xyz1:c_xyz2, .combine = 'rbind', .verbose = FALSE) %dopar% {
      Ess_MolVecMat_call(c)
    }
    Ess_MolVecMat <- matrix(Ess_MolVecMat, ncol = L_Essential_Elements)
    L_Ess_MolVecMat <- dim(Ess_MolVecMat)[1]
    print(paste0("There are ", L_Ess_MolVecMat, " essential elements combinations!"))
    ##
    print("Initiated enumerating molecular formulas!")
    MolVecMat <- foreach(counter = 1:L_Ess_MolVecMat, .verbose = FALSE) %dopar% {
      MolVecMat_call(counter)
    }
    MolVecMat <- matrix(MolVecMat, ncol = (L_ElementsAlphabetical + 1))
    L_MolVecMat <- dim(MolVecMat)[1]
    ##
    stopCluster(clust)
    ##
  } else if (osType == "Linux") {
    ##
    print("Initiated counting essential elements combinations!")
    Ess_MolVecMat <- do.call(rbind, mclapply(c_xyz1:c_xyz2, function(c) {
      Ess_MolVecMat_call(c)
    }, mc.cores = number_processing_threads))
    Ess_MolVecMat <- matrix(Ess_MolVecMat, ncol = L_Essential_Elements)
    L_Ess_MolVecMat <- dim(Ess_MolVecMat)[1]
    print(paste0("There are ", L_Ess_MolVecMat, " essential elements combinations!"))
    ##
    print("Initiated enumerating molecular formulas!")
    MolVecMat <- do.call(rbind, mclapply(1:L_Ess_MolVecMat, function(counter) {
      MolVecMat_call(counter)
    }, mc.cores = number_processing_threads))
    MolVecMat <- matrix(MolVecMat, ncol = (L_ElementsAlphabetical + 1))
    L_MolVecMat <- dim(MolVecMat)[1]
    ##
    closeAllConnections()
  }
  ##
  if (checkpoint_parameter == 1) {
    print(paste0("The approximate maximum number of the candidate molecular formulas is ", L_MolVecMat, " !"))
    ##
    print("The required time only to iterate loops through this chemical space is :")
    tr2_sub_pattern <- Sys.time()
    print(tr2_sub_pattern - tr1_sub_pattern)
  }
  return(checkpoint_parameter)
}
