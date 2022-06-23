UFAx_workflow <- function(spreadsheet) {
  ##
  exhaustive_chemical_enumeration_annotated_table <- c()
  gc()
  initiation_time <- Sys.time()
  message("Initiated testing the spreadsheet consistency!")
  ##
  checkpoint_parameter <- FALSE
  #
  if (length(spreadsheet) >= 4) {
    if (typeof(spreadsheet) == "list") {
      PARAM_ECE <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
    } else {
      message("The UFAx spreadsheet was not produced properly!")
    }
  } else if (length(spreadsheet) == 1) {
    if (typeof(spreadsheet) == "character") {
      if (file.exists(spreadsheet)) {
        spreadsheet <- readxl::read_xlsx(spreadsheet, sheet = "exhaustive_chemical_enumeration")
        PARAM_ECE <- cbind(spreadsheet[, 2], spreadsheet[, 4])
        checkpoint_parameter <- TRUE
      } else {
        message("The UFAx spreadsheet not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      message("The UFAx spreadsheet was not produced properly!")
    }
  } else {
    message("The UFAx spreadsheet was not produced properly!")
  }
  ##
  if (checkpoint_parameter == TRUE) {
    x0001 <- which(PARAM_ECE[, 1] == "ECE0001")
    address_pl <- gsub("\\", "/", PARAM_ECE[x0001, 2], fixed = TRUE)
    PARAM_ECE[x0001, 2] <- address_pl
    if (dir.exists(address_pl)) {
      ECE0002 <- PARAM_ECE[which(PARAM_ECE[, 1] == "ECE0002"), 2]
      pl_file <- paste0(address_pl, "/", ECE0002)
      if (!file.exists(pl_file)) {
        checkpoint_parameter <- FALSE
        message("ERROR!!! Problem with ECE0002! peaklist is not available!")
      }
    } else {
      checkpoint_parameter <- FALSE
      message("ERROR!!! Problem with ECE0001! Folder of peaklist is not available!")
    }
    peaklist <- IDSL.IPA::loadRdata(pl_file)
    peaklist <- matrix(peaklist, ncol = 24)
    n_peaks <- dim(peaklist)[1]
    ##
    ECE0003 <- PARAM_ECE[which(PARAM_ECE[, 1] == "ECE0003"), 2]
    if (is.na(ECE0003)) {
      checkpoint_parameter <- FALSE
      message("ERROR!!! Problem with ECE0003! This parameter should be 'All' or a vector of indices!")
    } else {
      if (gsub(" ", "", tolower(ECE0003)) == "all") {
        selected_IPA_peaks <- 1:n_peaks
        message("The enitre 12C m/z values in the peaklist were placed in the processing row! Annotated molecular formulas for peak IDs are kept in the 'log_ECE_annotation_' folder!")
      } else {
        selected_IPA_peaks <- tryCatch(eval(parse(text = paste0("c(", ECE0003, ")"))), error = function(e){NULL})
        if (is.null(selected_IPA_peaks) | (max(selected_IPA_peaks) > n_peaks)) {
          checkpoint_parameter <- FALSE
          message("ERROR!!! Problem with ECE0003! The range of indices are out of the peaklist dimension!")
        } else {
          message("The following peak IDs were selected for processing: ")
          for (id in 1:length(selected_IPA_peaks)) {
            message(paste0(selected_IPA_peaks[id], " - ", peaklist[selected_IPA_peaks[id], 3],  " - ", peaklist[selected_IPA_peaks[id], 8]))
          }
        }
      }
    }
    ##
    x0004 <- which(PARAM_ECE[, 1] == "ECE0004")
    address_mass_spec_file <- gsub("\\", "/", PARAM_ECE[x0004, 2], fixed = TRUE)
    PARAM_ECE[x0004, 2] <- address_mass_spec_file
    if (dir.exists(address_mass_spec_file)) {
      ECE0005 <- PARAM_ECE[which(PARAM_ECE[, 1] == "ECE0005"), 2]
      mass_spec_file <- paste0(address_mass_spec_file, "/", ECE0005)
      if (!file.exists(mass_spec_file)) {
        checkpoint_parameter <- FALSE
        message("ERROR!!! Problem with ECE0005! HRMS is not available!")
      }
    } else {
      checkpoint_parameter <- FALSE
      message("ERROR!!! Problem with ECE0004! Folder of HRMS file is not available!")
    }
    ##
    x0006 <- which(PARAM_ECE[, 1] == "ECE0006")
    output_path <- gsub("\\", "/", PARAM_ECE[x0006, 2], fixed = TRUE)
    PARAM_ECE[x0006, 2] <- output_path
    if (!dir.exists(output_path)) {
      tryCatch(dir.create(output_path), warning = function(w){message("WARNING!!! Problem with ECE0006! R can only create one folder!")})
      if (!dir.exists(output_path)) {
        checkpoint_parameter <- FALSE
      }
    }
    ##
    number_processing_threads <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0007'), 2])
    if (length(number_processing_threads) == 0) {
      message("ERROR!!! Problem with ECE0007! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (number_processing_threads >= 1) {
        if ((number_processing_threads %% 1) != 0) {
          message("ERROR!!! Problem with ECE0007! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        message("ERROR!!! Problem with ECE0007! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    r13c_threshold <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0008'), 2])
    if (length(r13c_threshold) == 0) {
      message("ERROR!!! Problem with ECE0008! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (r13c_threshold <= 0) {
        message("ERROR!!! Problem with ECE0008! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ######################### Chemical space ###################################
    B_MAX <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0009'), 2])
    if (length(B_MAX) == 0) {
      message("ERROR!!! Problem with ECE0009! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (B_MAX < 0) {
        message("ERROR!!! Problem with ECE0009! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    Br_MAX <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0010'), 2])
    if (length(Br_MAX) == 0) {
      message("ERROR!!! Problem with ECE0010! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (Br_MAX < 0) {
        message("ERROR!!! Problem with ECE0010! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    Cl_MAX <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0011'), 2])
    if (length(Cl_MAX) == 0) {
      message("ERROR!!! Problem with ECE0011! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (Cl_MAX < 0) {
        message("ERROR!!! Problem with ECE0011! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    S_MAX <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0012'), 2])
    if (length(S_MAX) == 0) {
      message("ERROR!!! Problem with ECE0012! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (S_MAX < 0) {
        message("ERROR!!! Problem with ECE0012! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    Si_MAX <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0013'), 2])
    if (length(Si_MAX) == 0) {
      message("ERROR!!! Problem with ECE0013! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (Si_MAX < 0) {
        message("ERROR!!! Problem with ECE0013! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    N_MAX <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0014'), 2])
    if (length(N_MAX) == 0) {
      message("ERROR!!! Problem with ECE0014! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (N_MAX < 0) {
        message("ERROR!!! Problem with ECE0014! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    As_MAX <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0015'), 2])
    if (length(As_MAX) == 0) {
      message("ERROR!!! Problem with ECE0015! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (As_MAX < 0) {
        message("ERROR!!! Problem with ECE0015! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    I_MAX <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0016'), 2])
    if (length(I_MAX) == 0) {
      message("ERROR!!! Problem with ECE0016! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (I_MAX < 0) {
        message("ERROR!!! Problem with ECE0016! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    O_MAX <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0017'), 2])
    if (length(O_MAX) == 0) {
      message("ERROR!!! Problem with ECE0017! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (O_MAX < 0) {
        message("ERROR!!! Problem with ECE0017! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    P_MAX <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0018'), 2])
    if (length(P_MAX) == 0) {
      message("ERROR!!! Problem with ECE0018! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (P_MAX < 0) {
        message("ERROR!!! Problem with ECE0018! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    Na_MAX <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0019'), 2])
    if (length(Na_MAX) == 0) {
      message("ERROR!!! Problem with ECE0019! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (Na_MAX < 0) {
        message("ERROR!!! Problem with ECE0019! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    K_MAX <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0020'), 2])
    if (length(K_MAX) == 0) {
      message("ERROR!!! Problem with ECE0020! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (K_MAX < 0) {
        message("ERROR!!! Problem with ECE0020! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ############################################################################
    ipw <- PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0021'), 2]
    if (ipw == "[M+H/K/Na]" | ipw == "[M+H]") {
      ipw_n <- -1
    } else if (ipw == "[M-H]") {
      ipw_n <- +1
    } else if (ipw == "[M]") {
      ipw_n <- 0
    } else {
      message("ERROR!!! Problem with ECE0021! This parameter should be any of '[M+H/K/Na]', '[M-H]', '[M]' ")
      checkpoint_parameter <- FALSE
    }
    ##
    peak_spacing <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0022'), 2])
    if (length(peak_spacing) == 0) {
      message("ERROR!!! Problem with ECE0022! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (peak_spacing < 0) {
        message("ERROR!!! Problem with ECE0022! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    int_cutoff_str <- PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0023'), 2]
    if (length(int_cutoff_str) == 0) {
      message("ERROR!!! Problem with ECE0023!")
      checkpoint_parameter <- FALSE
    }
    ##
    UFA_IP_memeory_variables <- tryCatch(eval(parse(text = paste0("c(", PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0024'), 2], ")"))), error = function(e){NULL})
    if (length(UFA_IP_memeory_variables) != 2) {
      message("ERROR!!! Problem with ECE0024! This parameter should be a vector of two positive numbers!")
      checkpoint_parameter <- FALSE
    }
    ##
    halogen_rule <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0025'), 2])
    if (length(halogen_rule) == 0) {
      message("ERROR!!! Problem with ECE0025! This parameter should be a number!")
      checkpoint_parameter <- FALSE
    } else {
      if (halogen_rule < 0) {
        message("ERROR!!! Problem with ECE0025! This parameter should be a number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    Na_K_x <- PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0026'), 2]
    if (Na_K_x == "1" | tolower(Na_K_x) == "t" | tolower(Na_K_x) == "true") {
      Na_K_rule <- TRUE
    } else if (Na_K_x == "" | Na_K_x == " " | Na_K_x == "0" | tolower(Na_K_x) == "f" | tolower(Na_K_x) == "false") {
      Na_K_rule <- FALSE
    } else {
      message("ERROR!!! Problem with ECE0026!")
      checkpoint_parameter <- FALSE
    }
    ##
    MaxR13C <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0027'), 2]) + r13c_threshold
    if (length(MaxR13C) == 0) {
      message("ERROR!!! Problem with ECE0027! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (MaxR13C <= 0) {
        message("ERROR!!! Problem with ECE0027! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ############### Molecular formula annotation criteria ######################
    mass_accuracy <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0028'), 2])
    if (length(mass_accuracy) == 0) {
      message("ERROR!!! Problem with ECE0028! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (mass_accuracy <= 0) {
        message("ERROR!!! Problem with ECE0028! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    maxNEME <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0029'), 2])
    if (length(maxNEME) == 0) {
      message("ERROR!!! Problem with ECE0029! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (maxNEME <= 0) {
        message("ERROR!!! Problem with ECE0029! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    minPCS <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0030'), 2])
    if (length(minPCS) == 0) {
      message("ERROR!!! Problem with ECE0030! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (minPCS <= 0) {
        message("ERROR!!! Problem with ECE0030! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    minNDCS <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0031'), 2])
    if (length(minNDCS) == 0) {
      message("ERROR!!! Problem with ECE0031! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (minNDCS <= 0) {
        message("ERROR!!! Problem with ECE0031! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    minRCS <- as.numeric(PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0032'), 2])
    if (length(minRCS) == 0) {
      message("ERROR!!! Problem with ECE0032! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (minRCS <= 0) {
        message("ERROR!!! Problem with ECE0032! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    Score_coeff <- tryCatch(eval(parse(text = PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0033'), 2])), error = function(e){NULL})
    if (is.null(Score_coeff)) {
      message("ERROR!!! Problem with ECE0033! This parameter should be a vector of five positive numbers!")
      checkpoint_parameter <- FALSE
    } else {
      if (length(Score_coeff) != 5) {
        message("ERROR!!! Problem with ECE0033! This parameter should be a vector of five positive numbers!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    ECE0034 <- PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0034'), 2]
    if (length(ECE0034) == 0) {
      message("ERROR!!! Problem with ECE0034!")
      checkpoint_parameter <- FALSE
    } else {
      if (!(tolower(ECE0034) == "yes" | tolower(ECE0034) == "no")) {
        message("ERROR!!! Problem with ECE0034!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (tolower(ECE0034) == "yes") {
      IonPathways <- tryCatch(eval(parse(text = paste0("c(", PARAM_ECE[which(PARAM_ECE[, 1] == 'ECE0035'), 2], ")"))), error = function(e){NULL})
      if (is.null(IonPathways)) {
        message("ERROR!!! Problem with ECE0035!")
        checkpoint_parameter <- FALSE
      }
      ##
      ECE0036 <- which(PARAM_ECE[, 1] == 'ECE0036')
      if (length(ECE0036) == 0) {
        message("ERROR!!! Problem with ECE0036! PubChem library data is not available! You should use the 'molecular_formula_library_generator' module to produce the molecular formula library!")
        checkpoint_parameter <- FALSE
      } else {
        PubChem_library_path <- gsub("\\", "/", PARAM_ECE[ECE0036, 2], fixed = TRUE)
        PARAM_ECE[ECE0036, 2] <- PubChem_library_path
        if (!file.exists(PubChem_library_path)) {
          message("ERROR!!! Problem with ECE0036! PubChem library data is not available! You should use the 'molecular_formula_library_generator' module to produce the molecular formula library!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ############################################################################
    if (checkpoint_parameter == TRUE) {
      message("The spreadsheet is consistent with the IDSL.UFAx workflow!")
      ##
      Elements <- c("As", "Br", "Cl", "Na", "Si", "B",  "C", "F",  "H",  "I",  "K",  "N",  "O",  "P",  "S") # DO NOT change this order!
      x_c_el <- which(Elements == "C") # index number of Carbon
      L_Elements <- length(Elements)
      EL <- element_sorter(ElementList = Elements, ElementOrder = "same")
      Elements_mass_abundance <- EL[[2]]
      valence_vec <- EL[[3]]
      Elements_mass <- sapply(1:L_Elements, function(el) {
        Elements_mass_abundance[[el]][[1]][1]
      })
      ##########################################################################
      get_formula_mz <- function(query_mz, qvec, carbon_masses, atom_vec, mass_accuracy) {
        ##
        set_1 <- comboGeneral(qvec,1,FALSE,constraintFun = "sum",
                              comparisonFun = c(">=","<="),
                              limitConstraints = c((query_mz - mass_accuracy), (query_mz + mass_accuracy)),
                              keepResults = TRUE)
        
        set_2 <- comboGeneral(qvec,2,FALSE,constraintFun = "sum",
                              comparisonFun = c(">=","<="),
                              limitConstraints = c((query_mz - mass_accuracy), (query_mz + mass_accuracy)),
                              keepResults = TRUE)
        set_2_a <- set_2
        set_2_a[!set_2_a %in% carbon_masses] <- 0
        set_2_a[set_2_a %in% carbon_masses] <- 1
        set_2 <- set_2[which(rowSums(set_2_a) > 0),]
        
        set_3 <- comboGeneral(qvec,3,FALSE,constraintFun = "sum",
                              comparisonFun = c(">=","<="),
                              limitConstraints = c((query_mz - mass_accuracy), (query_mz + mass_accuracy)),
                              keepResults = TRUE)
        set_3_a <- set_3
        set_3_a[!set_3_a %in% carbon_masses] <- 0
        set_3_a[set_3_a %in% carbon_masses] <- 1
        set_3 <- set_3[which(rowSums(set_3_a) > 0),]
        
        set_4 <- comboGeneral(qvec,4,FALSE,constraintFun = "sum",
                              comparisonFun = c(">=","<="),
                              limitConstraints = c((query_mz - mass_accuracy), (query_mz + mass_accuracy)),
                              keepResults = TRUE)
        set_4_a <- set_4
        set_4_a[!set_4_a %in% carbon_masses] <- 0
        set_4_a[set_4_a %in% carbon_masses] <- 1
        set_4 <- set_4[which(rowSums(set_4_a) > 0),]
        
        set_5 <- comboGeneral(qvec,5,FALSE,constraintFun = "sum",
                              comparisonFun = c(">=","<="),
                              limitConstraints = c((query_mz - mass_accuracy), (query_mz + mass_accuracy)),
                              keepResults = TRUE)
        set_5_a <- set_5
        set_5_a[!set_5_a %in% carbon_masses] <- 0
        set_5_a[set_5_a %in% carbon_masses] <- 1
        set_5 <- set_5[which(rowSums(set_5_a) > 0),]
        
        set_6 <- comboGeneral(qvec,6,FALSE,constraintFun = "sum",
                              comparisonFun = c(">=","<="),
                              limitConstraints = c((query_mz - mass_accuracy), (query_mz + mass_accuracy)),
                              keepResults = TRUE)
        set_6_a <- set_6
        set_6_a[!set_6_a %in% carbon_masses] <- 0
        set_6_a[set_6_a %in% carbon_masses] <- 1
        set_6 <- set_6[which(rowSums(set_6_a) > 0), ]
        ##
        get_formula_text <- function(sety) {
          if (is.matrix(sety) == F) {
            sety <- t(as.matrix(sety))
          }
          if (nrow(sety) > 0) {
            form_mat <- sapply(1:(ncol(sety) - 1), function(z) {
              as.character(atom_vec[as.character(sety[,z])])
            })
            if (!is.matrix(form_mat)) {
              form_mat <-  t(as.matrix(form_mat ))
            }
            apply(form_mat , 1, paste0, collapse = "")
          }
        }
        c(get_formula_text(set_1),get_formula_text(set_2),get_formula_text(set_3),get_formula_text(set_4),get_formula_text(set_5),get_formula_text(set_6))
      }
      ##########################################################################
      exhaustive_chemical_enumeration_call <- "exhaustive_chemical_enumeration_call <- function (i_mz) {
        ex_chem_enum_spa_annontated_table <- c()
        ##
        query_mz <- as.numeric(peaklist[i_mz, 8])
        R13C_PL <- as.numeric(peaklist[i_mz, 11])
        ##
        c_min <- max(c(1, floor((R13C_PL - r13c_threshold)/1.0816)))
        c_max <- floor((R13C_PL + r13c_threshold)/1.0816) + 1
        b_max <- min(c(B_MAX, max(c(0, floor((query_mz-c_min*Elements_mass[x_c_el])/Elements_mass[6])))))
        br_max <- min(c(Br_MAX, max(c(0, floor((query_mz-c_min*Elements_mass[x_c_el])/Elements_mass[2])))))
        cl_max <- min(c(Cl_MAX, max(c(0, floor((query_mz-c_min*Elements_mass[x_c_el])/Elements_mass[3])))))
        k_max <- K_MAX
        s_max <- min(c(S_MAX, max(c(0, floor((query_mz-c_min*Elements_mass[x_c_el])/Elements_mass[15])))))
        si_max <- min(c(Si_MAX, max(c(0, floor((query_mz-c_min*Elements_mass[x_c_el])/Elements_mass[5])))))
        n_max <- min(c(N_MAX, max(c(0, floor((query_mz-c_min*Elements_mass[x_c_el])/Elements_mass[12])))))
        as_max <- min(c(As_MAX, floor((query_mz-c_min*Elements_mass[x_c_el])/Elements_mass[1])))
        f_max <- min(c(2*c_max+3*n_max+6, max(c(0, floor((query_mz-c_min*Elements_mass[x_c_el])/Elements_mass[8])))))
        i_max <- min(c(I_MAX, floor((query_mz-c_min*Elements_mass[x_c_el])/Elements_mass[10])))
        h_max <- min(c(2*c_max+3*n_max+6, max(c(0, floor((query_mz-c_min*Elements_mass[x_c_el])/Elements_mass[9])))))
        na_max <- Na_MAX
        o_max <- min(c(O_MAX, max(c(0, floor((query_mz-c_min*Elements_mass[x_c_el])/Elements_mass[13])))))
        p_max <- min(c(P_MAX, max(c(0, floor((query_mz-c_min*Elements_mass[x_c_el])/Elements_mass[14])))))
        ##
        minfreq <- rep(0, L_Elements)
        minfreq[x_c_el] <- c_min        
        maxfreq <- c(as_max, br_max, cl_max, na_max, si_max, b_max,  c_max, f_max,  h_max,  i_max,  k_max,  n_max,  o_max,  p_max,  s_max)
        #
        qvec <- unlist(lapply(1:L_Elements, function(x) {
          xvec <- round(sapply(c(minfreq[x]:maxfreq[x]), function(yy) { yy*Elements_mass[x] }), digits = 3)
          names(xvec) <- sapply(c(minfreq[x]:maxfreq[x]), function(yy) { paste0(Elements[x],yy) })
          xvec
        }))
        qvec <- qvec[which(qvec != 0)]
        atom_vec <- names(qvec)
        names(atom_vec) <- as.character(qvec)
        carbon_masses <- Elements_mass[x_c_el]*(minfreq[x_c_el]:maxfreq[x_c_el])
        molecular_formula <- get_formula_mz(query_mz, qvec, carbon_masses, atom_vec, mass_accuracy)
        ########################################################################
        if (length(molecular_formula) > 0) {
          MoleFormVecMat <- do.call(rbind, lapply(1:length(molecular_formula), function (molf) {
            molvec <- formula_vector_generator(molecular_formula[molf], Elements, L_Elements)
            if (molvec[x_c_el] > 0) {
              molvec
            }
          }))
          ######################################################################
          if (length(MoleFormVecMat) > 0) {
            MoleFormVecMat <- matrix(MoleFormVecMat, ncol = L_Elements)
            br <- MoleFormVecMat[, 2]
            cl <- MoleFormVecMat[, 3]
            c <- MoleFormVecMat[, x_c_el] # x_c_el = 7
            f <- MoleFormVecMat[, 8]
            h <- MoleFormVecMat[, 9]
            i <- MoleFormVecMat[, 10]
            n <- MoleFormVecMat[, 12]
            #
            x_Condition <- which(((h + cl + br + f + i) >= (c/2 - n - 1) & (h + cl + br + f + i) <= (2*c + 3*n + 6)) == TRUE)
            ##
            if (Na_K_rule) {
              na <- MoleFormVecMat[, 4]
              k <- MoleFormVecMat[, 11]
              x_Condition2 <- which(((na + k) <= 1) == TRUE)
              x_Condition <- x_Condition[which((x_Condition%in%x_Condition2) == TRUE)]
            }
            ##
            if (halogen_rule > 0) {
              x_Condition3 <- which(((br + cl + f + i) <= halogen_rule) == TRUE)
              x_Condition <- x_Condition[which((x_Condition%in%x_Condition3) == TRUE)]
            }
            ##
            L_MoleFormVecMat <- length(x_Condition)
            if (L_MoleFormVecMat > 0) {
              MoleFormVecMat <- matrix(MoleFormVecMat[x_Condition, ], ncol = L_Elements)
              x_sen_rule <- sapply(1:L_MoleFormVecMat, function(el) {# Extended SENOIR rule
                extended_SENIOR_rule_check(MoleFormVecMat[el, ], valence_vec, ionization_correction = ipw_n)
              })
              if (length(which(x_sen_rule == TRUE)) > 0) {
                MoleFormVecMat <- matrix(MoleFormVecMat[x_sen_rule, ], ncol = L_Elements)
                ##
                MoleFormVecMat <- matrix(unique(as.matrix(MoleFormVecMat)), ncol = L_Elements) # To remove redundant rows
                ##
                L_MoleFormVecMat <- dim(MoleFormVecMat)[1]
                ##
                IP <- lapply(1:L_MoleFormVecMat, function (i_mat) {
                  ##
                  br <- MoleFormVecMat[i_mat, 2]
                  cl <- MoleFormVecMat[i_mat, 3]
                  si <- MoleFormVecMat[i_mat, 5]
                  b <- MoleFormVecMat[i_mat, 6]
                  c <- MoleFormVecMat[i_mat, x_c_el] # x_c_el = 7
                  k <- MoleFormVecMat[i_mat, 11]
                  s <- MoleFormVecMat[i_mat, 15]
                  ##
                  intensity_cutoff <- int_cutoff_str
                  isotopic_profile_calculator(MoleFormVecMat[i_mat, ], Elements_mass_abundance, peak_spacing, intensity_cutoff, UFA_IP_memeory_variables)
                })
                ##
                ip_db_mat <- do.call(rbind, lapply(1:L_MoleFormVecMat, function(i_mat) {
                  IPP <- IP[[i_mat]]
                  L_IPP <- length(IPP[, 2])
                  ##
                  r13c_ip <- 0
                  if (L_IPP > 1) {
                    M13C <- abs(IPP[, 1] - IPP[1, 1] - 1.00335484)
                    M13C <- M13C[2:L_IPP]
                    x_11 <- which.min(M13C)[1]
                    if (M13C[x_11] <= 0.015) {
                      x_11 <- x_11 + 1
                      r13c_ip <- IPP[x_11, 2]/IPP[1, 2]*100
                    }
                  }
                  c(IPP[1, 1], r13c_ip, L_IPP)
                }))
                ################################################################
                x_in_range <- which((abs(query_mz - ip_db_mat[, 1]) <= mass_accuracy) & (ip_db_mat[, 2] <= MaxR13C))
                L_x_in_range <- length(x_in_range)
                if (L_x_in_range > 0) {
                  IsotopicProfile_DataBase <- IP[x_in_range]
                  R13C_DataBase <- ip_db_mat[x_in_range, 2]
                  SizeIP_IsotopicProfile_DataBase <- ip_db_mat[x_in_range, 3]
                  ##
                  RangeScan <- peaklist[i_mz, 1]:peaklist[i_mz, 2]
                  NumberScans <- length(RangeScan)
                  mzList2 <- do.call(rbind, lapply(1:L_x_in_range, function(j) {
                    A <- c()
                    IsotopicProfile <- IsotopicProfile_DataBase[[j]]
                    size_IP <- SizeIP_IsotopicProfile_DataBase[j]
                    MW_exp <- matrix(rep(0, size_IP*NumberScans), ncol = NumberScans)
                    INT_exp <- MW_exp
                    for (sc in 1:NumberScans) {
                      PEAKS <- spectraList[[RangeScan[sc]]]
                      for (Iso in 1:size_IP) {
                        x_Iso <- which(abs(PEAKS[, 1] - IsotopicProfile[Iso, 1]) <= mass_accuracy)
                        if (length(x_Iso) > 0) {
                          if (length(x_Iso) > 1) {
                            x_Is0 <- which.min(abs(PEAKS[x_Iso, 1] - IsotopicProfile[Iso, 1]))
                            x_Iso <- x_Iso[x_Is0[1]]
                          }
                          MW_exp[Iso, sc] <- PEAKS[x_Iso, 1]
                          INT_exp[Iso, sc] <- PEAKS[x_Iso, 2]
                        }
                      }
                    }
                    sum_INT_exp <- rowSums(INT_exp)
                    x_INT_0 <- which(sum_INT_exp == 0)
                    if (length(x_INT_0) == 0) {
                      Ave_MW_exp <- rowSums(MW_exp*INT_exp)/sum_INT_exp
                      NEME <- sqrt(sum((Ave_MW_exp - IsotopicProfile[, 1])^2)/size_IP)*1000 # in mDa
                      if (NEME <= maxNEME) {
                        PCS <- sum(sum_INT_exp*IsotopicProfile[, 2])/sqrt(sum(sum_INT_exp^2)*sum(IsotopicProfile[, 2]^2))*1000 # in per-mille
                        if (PCS >= minPCS) {
                          MW_exp[which(MW_exp > 0)] <- 1
                          nd <- colSums(MW_exp)
                          Int_100 <- INT_exp[1, ]
                          max_Int <- max(Int_100)
                          x_80 <- which(Int_100/max_Int > 0.2)
                          NDCS <- length(which(nd[x_80] == size_IP))
                          if (NDCS >= minNDCS) {
                            L_80 <- x_80[length(x_80)] - x_80[1] + 1
                            RCS <- NDCS/L_80*100
                            if (RCS >= minRCS) {
                              R13C_IP <- R13C_DataBase[j]
                              IdentificationScore <- (size_IP^Score_coeff[1])*((PCS/1000)^Score_coeff[2])*((RCS/100)^Score_coeff[3])/((NEME/maxNEME)^Score_coeff[4])/(exp(abs(log(R13C_PL/R13C_IP)))^Score_coeff[5])
                              A <- c(i_mz, j,
                                     IsotopicProfile[1, 1],
                                     peaklist[i_mz, 8],
                                     peaklist[i_mz, 3],
                                     max_Int,
                                     NEME,
                                     PCS,
                                     R13C_PL,
                                     R13C_IP,
                                     NDCS,
                                     RCS,
                                     IdentificationScore)
                            }
                          }
                        }
                      }
                    }
                    A
                  }))
                  if (!is.null(mzList2)) {
                    mzList2 <- matrix(mzList2, ncol = 13)
                    mzList <- mzList2[order(mzList2[, 13], decreasing = TRUE), ]
                    mzList <- matrix(mzList, ncol = 13)
                    mzList <- matrix(cbind(mzList, 1:nrow(mzList)), ncol = 14)
                    mzList <- mzList[, -13]
                    mzList <- matrix(mzList, ncol = 13)
                    ##
                    mzList[, 3] <- round(mzList[, 3], 5)
                    mzList[, 4] <- round(mzList[, 4], 5)
                    mzList[, 5] <- round(mzList[, 5], 3)
                    mzList[, 6] <- round(mzList[, 6], 0)
                    mzList[, 7] <- round(mzList[, 7], 2)
                    mzList[, 8] <- round(mzList[, 8], 0)
                    mzList[, 9] <- round(mzList[, 9], 2)
                    mzList[, 10] <- round(mzList[, 10], 2)
                    mzList[, 12] <- round(mzList[, 12], 2)
                    ##
                    MolFlist <- hill_molecular_formula_printer(Elements, MoleFormVecMat[x_in_range[mzList[, 2]], ])
                    mzList[, 2] <- SizeIP_IsotopicProfile_DataBase[mzList[, 2]]
                    ##
                    ex_chem_enum_spa_annontated_table <- cbind(matrix(mzList[, 1:2], ncol = 2), MolFlist, matrix(mzList[, 3:13], ncol = 11))
                  }
                }
              }
            }
          }
        }
        save(ex_chem_enum_spa_annontated_table, file = paste0(output_path_log, i_mz, '_ECE_annotation.Rdata'))
        ##
        return(ex_chem_enum_spa_annontated_table)
      }"
      ##########################################################################
      exhaustive_chemical_enumeration_call <- gsub("int_cutoff_str", int_cutoff_str, exhaustive_chemical_enumeration_call)
      eval(parse(text = exhaustive_chemical_enumeration_call))
      ##########################################################################
      outputer003 <- IDSL.IPA::IPA_MSdeconvoluter(HRMS_path = mass_spec_file, MSfile = "")
      spectraList <- outputer003[[1]]
      ## Creating the log folder
      output_path_log <- paste0(output_path, "/log_ECE_annotation_", ECE0005, "/")
      if (!dir.exists(output_path_log)) {
        tryCatch(dir.create(output_path_log), warning = function(w){stop("Can't create the log folder!!!")})
      }
      ##
      message("Initiated the exhaustive chemical enumeration analysis!!!")
      ##
      if (number_processing_threads == 1) {
        ##
        exhaustive_chemical_enumeration_annotated_table <- do.call(rbind, lapply(selected_IPA_peaks, function(i_mz) {
          exhaustive_chemical_enumeration_call(i_mz)
        }))
        ##
      } else {
        ##
        osType <- Sys.info()[['sysname']]
        if (osType == "Windows") {
          clust <- makeCluster(number_processing_threads)
          registerDoParallel(clust)
          ##
          exhaustive_chemical_enumeration_annotated_table <- foreach(i_mz = selected_IPA_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
            exhaustive_chemical_enumeration_call(i_mz)
          }
          ##
          stopCluster(clust)
          ##
        } else if (osType == "Linux") {
          ##
          exhaustive_chemical_enumeration_annotated_table <- do.call(rbind, mclapply(selected_IPA_peaks, function(i_mz) {
            exhaustive_chemical_enumeration_call(i_mz)
          }, mc.cores = number_processing_threads))
          closeAllConnections()
        }
      }
      ##
      exhaustive_chemical_enumeration_annotated_table <- data.frame(exhaustive_chemical_enumeration_annotated_table)
      colnames(exhaustive_chemical_enumeration_annotated_table) <- c("PeakID", "sizeIP", "FormulaIon", "m/z Isotopic Profile", "m/z peaklist", "RT", "PeakHeight", "NEME (mDa)", "PCS", "R13C peakList", "R13C Isotopic Profile", "NDCS", "RCS (%)", "Rank")
      rownames(exhaustive_chemical_enumeration_annotated_table) <- c()
      if (tolower(ECE0034) == "yes") {
        message("Initiated searching in the library of molecular formula!!!")
        MF_library <- IDSL.IPA::loadRdata(PubChem_library_path)
        MF_ex_chem_enum <- exhaustive_chemical_enumeration_annotated_table[, 3]
        MF_library_matched <- UFAx_molecular_formula_library_search(MF_ex_chem_enum, IonPathways, Elements, MF_library, number_processing_threads)
        exhaustive_chemical_enumeration_annotated_table <- cbind(exhaustive_chemical_enumeration_annotated_table, MF_library_matched)
        message("Completed searching in the library of molecular formula!!!")
      }
      save(exhaustive_chemical_enumeration_annotated_table, file = paste0(output_path, "/exhaustive_chemical_enumeration_annotated_table_", ECE0005, ".Rdata"))
      write.csv(exhaustive_chemical_enumeration_annotated_table, file = paste0(output_path, "/exhaustive_chemical_enumeration_annotated_table_", ECE0005, ".csv"))
      message("Completed the exhaustive chemical enumeration analysis!!!")
      required_time <- Sys.time() - initiation_time
      print(required_time)
    } else {
      message("Please visit   https://ufa.idsl.me    for instructions!!!")
    }
  }
  return(exhaustive_chemical_enumeration_annotated_table)
}