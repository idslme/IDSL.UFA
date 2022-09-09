detect_formula_sets  <- function(molecular_formulas, ratio_delta_HBrClFI_C, mixed.HBrClFI.allowed, min_molecular_formula_class, max_number_formula_class, number_processing_threads = 1) {
  ##
  classes_formula_matrix <- NULL
  ##
  molecular_formulas <- unique(molecular_formulas)
  ## To detect halogenated-saturated compounds such as PFOS or mixed halogenated compounds with hydrogen
  if (mixed.HBrClFI.allowed == TRUE) {
    N_mixed.HBrClFI.allowed <- 0
  } else if (mixed.HBrClFI.allowed == FALSE) {
    N_mixed.HBrClFI.allowed <- 1
  }
  ##
  MMFC <- max(2, min_molecular_formula_class - 1)
  MNC <- max_number_formula_class + 1
  ##############################################################################
  Elements <- c("As", "Br", "Cl", "Na", "Se", "Si", "B", "C", "F", "H", "I", "K", "N", "O", "P", "S")
  L_Elements <- length(Elements)
  ##
  molecular_formulasMat_call <- function(k) {
    mol_vec <- formula_vector_generator(molecular_formulas[k], Elements, L_Elements)
    x_neg <- which(mol_vec < 0)
    if (length(x_neg) > 0) {
      mol_vec <- NULL
    }
    mol_vec
  }
  ##
  if (number_processing_threads == 1) {
    molecular_formulasMat <- do.call(rbind, lapply(1:length(molecular_formulas), function(k) {
      molecular_formulasMat_call(k)
    }))
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      molecular_formulasMat <- foreach(k = 1:length(molecular_formulas), .combine = 'rbind', .verbose = FALSE) %dopar% {
        molecular_formulasMat_call(k)
      }
      ##
      stopCluster(clust)
      ##
    } else if (osType == "Linux") {
      ##
      molecular_formulasMat <- do.call(rbind, mclapply(1:length(molecular_formulas), function(k) {
        molecular_formulasMat_call(k)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
    }
  }
  ##
  molecular_formulasMat <- unique(as.matrix(molecular_formulasMat)) # To remove redundant rows
  #
  x_c_el <- which(Elements == "C")
  x_h_el <- which(Elements == "H")
  x_br_el <- which(Elements == "Br")
  x_cl_el <- which(Elements == "Cl")
  x_f_el <- which(Elements == "F")
  x_i_el <- which(Elements == "I")
  #
  molecular_formulasMat <- molecular_formulasMat[order(molecular_formulasMat[, x_c_el]), ]
  colnames(molecular_formulasMat) <- Elements
  ##
  Backbone_Elements <- c("C", "H", "Br", "Cl", "F", "I")
  MolFormMat1 <- do.call(cbind, lapply(1:L_Elements, function(k) {
    i_el <- molecular_formulasMat[, k]
    x_el <- which(Backbone_Elements == Elements[k])
    if (length(x_el) > 0) {
      if (x_el == x_c_el) {
        i_el[i_el > 0] <- 1
      } else {
        i_el[i_el > 0] <- N_mixed.HBrClFI.allowed # to to merge h+br+cl+f+i and 1 to separate them
      }
    }
    i_el
  }))
  #
  colnames(MolFormMat1) <- Elements
  molecular_formulas1 <- hill_molecular_formula_printer(Elements, MolFormMat1, number_processing_threads)
  molecular_formulas1[is.na(molecular_formulas1)] <- "NoOrganicGroup"
  unique_molecular_formulas1 <- unique(molecular_formulas1)
  ##
  if (ratio_delta_HBrClFI_C > 0) {
    I_call <- function(k) {
      J <- list()
      class_index <- which(molecular_formulas1 == unique_molecular_formulas1[k])
      L_class_index <- length(class_index)
      if (L_class_index > MMFC) {
        molecular_formulasMat_subset <- molecular_formulasMat[class_index, ]
        hbrclfi <- molecular_formulasMat_subset[, x_h_el] + molecular_formulasMat_subset[, x_br_el] + molecular_formulasMat_subset[, x_cl_el] + molecular_formulasMat_subset[, x_f_el] + molecular_formulasMat_subset[, x_i_el]
        ccc <- molecular_formulasMat_subset[, x_c_el]
        Counter_C <- 0
        for (j in 1:(L_class_index - 1)) {
          if (ccc[j] != 0) {
            a1 <- hbrclfi[j]
            K <- NULL
            for (k in (j + 1):L_class_index) {
              if ((ccc[k] != ccc[j]) & (ccc[k] > 0) & (ccc[j] > 0)) {
                an <- hbrclfi[k]
                Cn_chain <- (an - a1)/(ccc[k] - ccc[j])
                if((Cn_chain == ratio_delta_HBrClFI_C)) { # The golden ratio of difference hbrclfi/c
                  K <- c(K, class_index[k])
                  ccc[k] <- 0
                  hbrclfi[k] <- 0
                }
              }
            }
            ccc[j] <- 0
            hbrclfi[j] <- 0
            if (length(K) >= MMFC) { # minimum number of congeners
              Counter_C <- Counter_C + 1
              J[[Counter_C]] <- c(class_index[j], K)
            }
          }
        }
      }
      J
    }
  }
  ##
  if (ratio_delta_HBrClFI_C == 0) {
    I_call <- function(k) {
      J <- list()
      class_index <- which(molecular_formulas1 == unique_molecular_formulas1[k])
      L_class_index <- length(class_index)
      if (L_class_index > MMFC) {
        molecular_formulasMat_subset <- molecular_formulasMat[class_index, ]
        hbrclfi <- molecular_formulasMat_subset[, x_h_el] + molecular_formulasMat_subset[, x_br_el] + molecular_formulasMat_subset[, x_cl_el] + molecular_formulasMat_subset[, x_f_el] + molecular_formulasMat_subset[, x_i_el]
        ccc <- molecular_formulasMat_subset[, x_c_el]
        Counter_C <- 0
        for (j in 1:(L_class_index - 1)) {
          if (ccc[j] != 0) {
            a1 <- hbrclfi[j]
            K <- NULL
            for (k in (j + 1):L_class_index) {
              if ((ccc[k] == ccc[j]) & (hbrclfi[k] == hbrclfi[j])  & (ccc[k] > 0) & (ccc[j] > 0)) {
                K <- c(K, class_index[k])
                ccc[k] <- 0
                hbrclfi[k] <- 0
              }
            }
            ccc[j] <- 0
            hbrclfi[j] <- 0
            if (length(K) >= MMFC) { # minimum number of congeners
              Counter_C <- Counter_C + 1
              J[[Counter_C]] <- c(class_index[j], K)
            }
          }
        }
      }
      J
    }
  }
  ##
  if (ratio_delta_HBrClFI_C > 0) {
    classes_call <- function(k) {
      index_molf <- I[[k]]
      if (length(index_molf) > 0) {
        MolVec <- rbind(MolFormMat1[index_molf[1], ], molecular_formulasMat[index_molf, ])
        molf1 <- hill_molecular_formula_printer(Elements, MolVec)
        molf1[is.na(molf1)] <- "NoOrganicGroup"
        molf1
      }
    }
  } else if (ratio_delta_HBrClFI_C == 0) {
    classes_call <- function(k) {
      index_molf <- I[[k]]
      if (length(index_molf) > 0) {
        NumC <- molecular_formulasMat[index_molf[1], x_c_el]
        sumX <- molecular_formulasMat[index_molf[1], x_h_el] + molecular_formulasMat[index_molf[1], x_br_el] + molecular_formulasMat[index_molf[1], x_cl_el] + molecular_formulasMat[index_molf[1], x_f_el] + molecular_formulasMat[index_molf[1], x_i_el]
        vec_non0 <- molecular_formulasMat[index_molf[1], ]
        vec_non0[c(x_c_el, x_h_el, x_br_el, x_cl_el, x_f_el, x_i_el)] <- 0
        molf0 <- hill_molecular_formula_printer(Elements, vec_non0)
        if (is.na(molf0)) {
          molf0 <- paste0("C", NumC, "X", sumX)
        } else {
          molf0 <- paste0("C", NumC, "X", sumX, molf0)
        }
        MolVec <- molecular_formulasMat[index_molf, ]
        molf1 <- c(molf0, hill_molecular_formula_printer(Elements, MolVec))
        molf1[is.na(molf1)] <- "NoOrganicGroup"
        molf1
      }
    }
  }
  ##
  classes_formula_matrix_call <- function(k) {
    A <- NULL
    L_C <- length(classes[[k]])
    if (L_C < MNC) {
      A <- matrix(c(classes[[k]], rep("", (MNC - L_C))), ncol = MNC)
    } else {
      A <- matrix(classes[[k]][1:MNC], ncol = MNC)
    }
    A
  }
  ##
  if (number_processing_threads == 1) {
    ##
    I <-  lapply(1:length(unique_molecular_formulas1), function(k) {
      I_call(k)
    })
    I <- unlist(I, recursive = FALSE)
    ##
    if (length(I) > 0) {
      classes <-  lapply(1:length(I), function(k) {
        classes_call(k)
      })
      ##
      classes_formula_matrix <- do.call(rbind, lapply(1:length(classes), function(k) {
        classes_formula_matrix_call(k)
      }))
    }
    ##
  } else {
    if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      I <- foreach(k = 1:length(unique_molecular_formulas1), .verbose = FALSE) %dopar% {
        I_call(k)
      }
      I <- unlist(I, recursive = FALSE)
      ##
      if (length(I) > 0) {
        classes <- foreach(k = 1:length(I), .verbose = FALSE) %dopar% {
          classes_call(k)
        }
        ##
        classes_formula_matrix <- foreach(k = 1:length(classes), .combine = 'rbind', .verbose = FALSE) %dopar% {
          classes_formula_matrix_call(k)
        }
      }
      ##
      stopCluster(clust)
      ##
    } else if (osType == "Linux") {
      ##
      I <-  mclapply(1:length(unique_molecular_formulas1), function(k) {
        I_call(k)
      }, mc.cores = number_processing_threads)
      I <- unlist(I, recursive = FALSE)
      ##
      if (length(I) > 0) {
        classes <-  mclapply(1:length(I), function(k) {
          classes_call(k)
        }, mc.cores = number_processing_threads)
        ##
        classes_formula_matrix <- do.call(rbind, mclapply(1:length(classes), function(k) {
          classes_formula_matrix_call(k)
        }, mc.cores = number_processing_threads))
      }
      ##
      closeAllConnections()
    }
  }
  ##
  Xcolnames <- do.call(c, lapply(1:max_number_formula_class, function(i) {
    paste0("Molecularformula_", i)
  }))
  ##
  classes_formula_matrix <- data.frame(classes_formula_matrix)
  colnames(classes_formula_matrix) <- c("Class_structure", Xcolnames)
  rownames(classes_formula_matrix) <- NULL
  ##
  return(classes_formula_matrix)
}
