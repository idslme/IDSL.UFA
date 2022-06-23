isotopic_profile_molecular_formula_feeder <- function(molecular_formula, peak_spacing = 0, intensity_cutoff_str = 1, UFA_IP_memeory_variables = c(1e30, 1e-12), IonPathways = "[M]+", number_processing_threads = 1) {
  EL <- element_sorter()
  Elements <- EL[[1]]
  Elements_mass_abundance <- EL[[2]]
  L_Elements <- length(Elements)
  ##
  x_el_c <- which(Elements == "C")
  x_el_b <- which(Elements == "B")
  x_el_br <- which(Elements == "Br")
  x_el_cl <- which(Elements == "Cl")
  x_el_k <- which(Elements == "K")
  x_el_s <- which(Elements == "S")
  x_el_se <- which(Elements == "Se")
  x_el_si <- which(Elements == "Si")
  ## IonPathway = [Coeff*M+CO2-H2O+Na-KO2+HCl-...] # Coeff should be an integer between 1-9
  IonPW_DC <- ionization_pathway_deconvoluter(IonPathways, Elements)
  L_PW <- length(IonPathways)
  ##
  molf_deconvoluter <- function (i_molf) {
    molf_deconvoluter_ipw <- c()
    FormulaVector <- formula_vector_generator(molecular_formula[i_molf], Elements, L_Elements)
    if (FormulaVector[x_el_c] > 0) {
      molf_deconvoluter_ipw <- do.call(rbind, lapply(1:L_PW, function (pathway) {
        molv_ipw <- c()
        IonPW <- IonPW_DC[[pathway]]
        Ion_coeff <- IonPW[[1]]
        Ion_adduct <- IonPW[[2]]
        MoleFormVec <- Ion_coeff*FormulaVector + Ion_adduct
        x_neg <- which(MoleFormVec < 0)
        if (length(x_neg) == 0) {
          molv_ipw <- MoleFormVec
        }
        molv_ipw
      }))
    }
    molf_deconvoluter_ipw
  }
  ##
  IP_calculator <- "IP_calculator <- function(i_mat) {
    ##
    c <- MoleFormVecMat[i_mat, x_el_c]
    b <- MoleFormVecMat[i_mat, x_el_b]
    br <- MoleFormVecMat[i_mat, x_el_br]
    cl <- MoleFormVecMat[i_mat, x_el_cl]
    k <- MoleFormVecMat[i_mat, x_el_k]
    s <- MoleFormVecMat[i_mat, x_el_s]
    se <- MoleFormVecMat[i_mat, x_el_se]
    si <- MoleFormVecMat[i_mat, x_el_si]
    ##
    intensity_cutoff <- intensity_cutoff_str
    IPP <- isotopic_profile_calculator(MoleFormVecMat[i_mat, ], Elements_mass_abundance, peak_spacing, intensity_cutoff, UFA_IP_memeory_variables)
    IPP[, 1] <- round(IPP[, 1], 6)
    IPP[, 2] <- round(IPP[, 2], 3)
    return(IPP)
  }"
  IP_calculator <- gsub("intensity_cutoff_str", intensity_cutoff_str, IP_calculator)
  eval(parse(text = IP_calculator))
  ##
  ip_db_function <- function(i) {
    IPP <- IsotopicProfileList[[i]]
    x_100 <- which.max(IPP[, 2])
    L_IPP <- length(IPP[, 2])
    ##
    r13c_ip <- 0
    if (L_IPP > x_100) {
      M13C <- abs(IPP[, 1] - IPP[x_100, 1] - 1.00335484)
      M13C <- M13C[(x_100 + 1):L_IPP]
      x_101 <- which.min(M13C)[1]
      if (M13C[x_101] <= 0.015) {
        x_101 <- x_101 + x_100
        r13c_ip <- IPP[x_101, 2]/IPP[x_100, 2]*100
      }
    }
    c(IPP[x_100, 1], r13c_ip, x_100, L_IPP)
  }
  ##
  L_MolF <- length(molecular_formula)
  print("Initiated deconvoluting molecular formulas!")
  ##
  if (number_processing_threads == 1) {
    ##
    MoleFormVecMat <- do.call(rbind, lapply(1:L_MolF, function(counter) {
      molf_deconvoluter(counter)
    }))
    if (is.null(MoleFormVecMat)) {
      stop("Molecular formulas are not consistent with the ionization pathways!")
    }
    MoleFormVecMat <- unique(as.matrix(MoleFormVecMat)) # To remove redundant rows
    print("Completed deconvoluting molecular formulas!")
    ##
    L_MoleFormVecMat <- dim(MoleFormVecMat)[1]
    print(paste0("There are ", L_MoleFormVecMat, " molecular formula ions for isotopic profile calculations!"))
    print("Initiated calculating isotopic profiles!")
    IsotopicProfileList <- lapply(1:L_MoleFormVecMat, function(counter) {
      IP_calculator(counter)
    })
    print("Completed calculating isotopic profiles!")
    ##
    print("Initiated calculating the database parameters!")
    ip_db_mat <- do.call(rbind, lapply(1:L_MoleFormVecMat, function(counter) {
      ip_db_function(counter)
    }))
    ##
  } else {
    osType <- Sys.info()[['sysname']]
    if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      MoleFormVecMat <- foreach(counter = 1:L_MolF, .combine = 'rbind', .verbose = FALSE) %dopar% {
        molf_deconvoluter(counter)
      }
      if (is.null(MoleFormVecMat)) {
        stop("Molecular formulas are not consistent with the ionization pathways!")
      }
      MoleFormVecMat <- unique(as.matrix(MoleFormVecMat)) # To remove redundant rows
      print("Completed deconvoluting molecular formulas!")
      ##
      L_MoleFormVecMat <- dim(MoleFormVecMat)[1]
      print(paste0("There are ", L_MoleFormVecMat, " molecular formula ions for isotopic profile calculations!"))
      print("Initiated calculating isotopic profiles!")
      IsotopicProfileList <- foreach(counter = 1:L_MoleFormVecMat, .verbose = FALSE) %dopar% {
        IP_calculator(counter)
      }
      print("Completed calculating isotopic profiles!")
      ##
      print("Initiated calculating the database parameters!")
      ip_db_mat <- foreach(counter = 1:L_MoleFormVecMat, .combine = 'rbind', .verbose = FALSE) %dopar% {
        ip_db_function(counter)
      }
      ##
      stopCluster(clust)
      ##
    } else if (osType == "Linux") {
      ##
      MoleFormVecMat <- do.call(rbind, mclapply(1:L_MolF, function(counter) {
        molf_deconvoluter(counter)
      }, mc.cores = number_processing_threads))
      if (is.null(MoleFormVecMat)) {
        stop("Molecular formulas are not consistent with the ionization pathways!")
      }
      MoleFormVecMat <- unique(as.matrix(MoleFormVecMat)) # To remove redundant rows
      print("Completed deconvoluting molecular formulas!")
      ##
      L_MoleFormVecMat <- dim(MoleFormVecMat)[1]
      print(paste0("There are ", L_MoleFormVecMat, " molecular formula ions for isotopic profile calculations!"))
      print("Initiated calculating isotopic profiles!")
      IsotopicProfileList <- mclapply(1:L_MoleFormVecMat, function(counter) {
        IP_calculator(counter)
      }, mc.cores = number_processing_threads)
      print("Completed calculating isotopic profiles!")
      ##
      print("Initiated calculating the database parameters!")
      ip_db_mat <- do.call(rbind, mclapply(1:L_MoleFormVecMat, function(counter) {
        ip_db_function(counter)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
    }
  }
  print("Completed calculating the database parameters!")
  ##############################################################################
  x_element_non0 <- do.call(c, lapply(1:L_Elements, function(counter) {
    x_non0 <- which(MoleFormVecMat[, counter] > 0)
    if (length(x_non0) > 0) {
      counter
    }
  }))
  #
  IP_library <- list(Elements = Elements[x_element_non0], MoleFormVecMat[, x_element_non0])
  MoleFormVecMat <- 0
  names(IP_library) <- c("Elements", "MolecularFormulaMatrix")
  ##
  ip_db_mat <- matrix(ip_db_mat, ncol = 4)
  IP_Mass <- ip_db_mat[, 1]
  R13C_IP <- ip_db_mat[, 2]
  Index_MAIso <- ip_db_mat[, 3]
  IP_size <- ip_db_mat[, 4]
  ##
  IPDB <- list(IP_Mass, IP_library, IsotopicProfileList, R13C_IP, Index_MAIso, IP_size)
  ##
  return(IPDB)
}
