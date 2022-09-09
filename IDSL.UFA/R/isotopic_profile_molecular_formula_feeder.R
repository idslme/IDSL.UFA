isotopic_profile_molecular_formula_feeder <- function(molecular_formula, peak_spacing = 0, intensity_cutoff_str = 1, UFA_IP_memeory_variables = c(1e30, 1e-12, 10), IonPathways = "[M]+", number_processing_threads = 1) {
  ##
  IsoProfInf100 <- matrix(c(Inf, 100), ncol = 2)
  ##
  EL <- element_sorter()
  Elements <- EL[["Elements"]]
  Elements_mass_abundance <- EL[["massAbundanceList"]]
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
  ##
  molf_deconvoluter <- function(i_molf) {
    FormulaVector <- formula_vector_generator(i_molf, Elements, L_Elements)
    if (FormulaVector[x_el_c] > 0) {
      do.call(rbind, lapply(IonPW_DC, function(pathway) {
        Ion_coeff <- pathway[[1]]
        Ion_adduct <- pathway[[2]]
        MoleFormVec <- Ion_coeff*FormulaVector + Ion_adduct
        x_neg <- which(MoleFormVec < 0)
        if (length(x_neg) == 0) {
          MoleFormVec
        }
      }))
    }
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
    ##
    IPP <- tryCatch(isotopic_profile_calculator(MoleFormVecMat[i_mat, ], Elements_mass_abundance, peak_spacing, intensity_cutoff, UFA_IP_memeory_variables),
                    error = function(e) {IsoProfInf100},
                    warning = function(w) {IsoProfInf100})
    ##
    IPP[, 1] <- round(IPP[, 1], 6)
    IPP[, 2] <- round(IPP[, 2], 3)
    ##
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
    c(IPP[x_100, 1], IP_R13C, x_100, L_IPP)
  }
  ##
  print("Initiated deconvoluting molecular formulas!")
  ##
  if (number_processing_threads == 1) {
    ##
    MoleFormVecMat <- do.call(rbind, lapply(molecular_formula, function(i) {
      molf_deconvoluter(i)
    }))
    if (is.null(MoleFormVecMat)) {
      stop("Molecular formulas are not consistent with the ionization pathways!")
    }
    MoleFormVecMat <- unique(as.matrix(MoleFormVecMat)) # To remove redundant rows
    print("Completed deconvoluting molecular formulas!")
    ##
    L_MoleFormVecMat <- dim(MoleFormVecMat)[1]
    print(paste0("There are ", L_MoleFormVecMat, " molecular formula ions for isotopic profile calculations!"))
    ##
    print("Initiated calculating isotopic profiles!")
    progressBARboundaries <- txtProgressBar(min = 0, max = L_MoleFormVecMat, initial = 0, style = 3)
    #
    IsotopicProfileList <- lapply(1:L_MoleFormVecMat, function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      IP_calculator(i)
    })
    close(progressBARboundaries)
    print("Completed calculating isotopic profiles!")
    ##
    print("Initiated calculating the database parameters!")
    ip_db_mat <- do.call(rbind, lapply(1:L_MoleFormVecMat, function(i) {
      ip_db_function(i)
    }))
    ##
  } else {
    osType <- Sys.info()[['sysname']]
    if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      MoleFormVecMat <- foreach(i = molecular_formula, .combine = 'rbind', .verbose = FALSE) %dopar% {
        molf_deconvoluter(i)
      }
      if (is.null(MoleFormVecMat)) {
        stop("Molecular formulas are not consistent with the ionization pathways!")
      }
      MoleFormVecMat <- unique(as.matrix(MoleFormVecMat)) # To remove redundant rows
      print("Completed deconvoluting molecular formulas!")
      ##
      L_MoleFormVecMat <- dim(MoleFormVecMat)[1]
      print(paste0("There are ", L_MoleFormVecMat, " molecular formula ions for isotopic profile calculations!"))
      ##
      print("Initiated calculating isotopic profiles!")
      IsotopicProfileList <- foreach(i = 1:L_MoleFormVecMat, .verbose = FALSE) %dopar% {
        IP_calculator(i)
      }
      print("Completed calculating isotopic profiles!")
      ##
      print("Initiated calculating the database parameters!")
      ip_db_mat <- foreach(i = 1:L_MoleFormVecMat, .combine = 'rbind', .verbose = FALSE) %dopar% {
        ip_db_function(i)
      }
      ##
      stopCluster(clust)
      ##
    } else if (osType == "Linux") {
      ##
      MoleFormVecMat <- do.call(rbind, mclapply(molecular_formula, function(i) {
        molf_deconvoluter(i)
      }, mc.cores = number_processing_threads))
      if (is.null(MoleFormVecMat)) {
        stop("Molecular formulas are not consistent with the ionization pathways!")
      }
      MoleFormVecMat <- unique(as.matrix(MoleFormVecMat)) # To remove redundant rows
      print("Completed deconvoluting molecular formulas!")
      ##
      L_MoleFormVecMat <- dim(MoleFormVecMat)[1]
      print(paste0("There are ", L_MoleFormVecMat, " molecular formula ions for isotopic profile calculations!"))
      ##
      print("Initiated calculating isotopic profiles!")
      IsotopicProfileList <- mclapply(1:L_MoleFormVecMat, function(i) {
        IP_calculator(i)
      }, mc.cores = number_processing_threads)
      print("Completed calculating isotopic profiles!")
      ##
      print("Initiated calculating the database parameters!")
      ip_db_mat <- do.call(rbind, mclapply(1:L_MoleFormVecMat, function(i) {
        ip_db_function(i)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
    }
  }
  print("Completed calculating the database parameters!")
  ##############################################################################
  x_element_non0 <- do.call(c, lapply(1:L_Elements, function(i) {
    x_non0 <- which(MoleFormVecMat[, i] > 0)
    if (length(x_non0) > 0) {
      i
    }
  }))
  #
  IP_library <- list(Elements = Elements[x_element_non0], MoleFormVecMat[, x_element_non0])
  MoleFormVecMat <- NULL
  names(IP_library) <- c("Elements", "MolecularFormulaMatrix")
  ##
  ip_db_mat <- matrix(ip_db_mat, ncol = 4)
  IP_Mass <- ip_db_mat[, 1]
  IP_R13C <- ip_db_mat[, 2]
  Index_MAIso <- ip_db_mat[, 3]
  IP_size <- ip_db_mat[, 4]
  ##
  IDroundMass <- cbind(round(IP_Mass, digits = 2), seq(1, L_MoleFormVecMat, 1))
  IDroundMass <- IDroundMass[order(IDroundMass[, 1], decreasing = FALSE), ]
  xDiff <- c(0, which(abs(diff(IDroundMass[, 1])) > 0), L_MoleFormVecMat)
  LDiff <- length(xDiff) - 1
  AggregatedList <- lapply(1:LDiff, function(i) {
    IDroundMass[(xDiff[i] + 1):xDiff[i + 1], 2]
  })
  names(AggregatedList) <- IDroundMass[(xDiff[1:LDiff] + 1), 1]
  ##
  IPDB <- list(AggregatedList, IP_Mass, IP_library, IsotopicProfileList, IP_R13C, Index_MAIso, IP_size)
  names(IPDB) <- c("AggregatedList", "MassMAIso", "MolecularFormulaDB", "IsotopicProfile", "R13C", "IndexMAIso", "IPsize")
  ##
  return(IPDB)
}
