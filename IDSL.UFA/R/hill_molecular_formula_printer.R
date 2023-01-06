hill_molecular_formula_printer <- function(Elements, MolVecMat, number_processing_threads = 1) {
  ##
  LElements <- length(Elements)
  HillElements <- sort(Elements, decreasing = FALSE)
  ##
  ##############################################################################
  ##
  x_c <- which(HillElements == "C")
  x_12c <- which(HillElements == "[12]C")
  x_13c <- which(HillElements == "[13]C")
  labeledCarbonCheck <- ((length(x_12c) == 1) | (length(x_13c) == 1))
  carbonCheck <- ((length(x_c) == 1) | labeledCarbonCheck)
  ##
  x_h <- which(HillElements == "H")
  x_2h <- which(HillElements == "[2]H")
  x_d <- which(HillElements == "D")
  labeledHydrogenCheck <- ((length(x_2h) == 1) | (length(x_d) == 1))
  hydrogenCheck <- ((length(x_h) == 1) | labeledHydrogenCheck)
  ##
  ##############################################################################
  ##
  if (carbonCheck & hydrogenCheck) {
    if (labeledCarbonCheck | labeledHydrogenCheck) {
      x_CH <- c(x_c, x_12c, x_13c, x_h, x_d, x_2h)
      carbonHydrogenVector <- do.call(c, lapply(x_CH, function(i){
        HillElements[i]
      }))
      ##
      HillElements <- HillElements[-x_CH]
      HillElements <- c(carbonHydrogenVector, HillElements)
    } else {
      x_CH <- c(x_c, x_h)
      HillElements <- HillElements[-x_CH]
      HillElements <- c("C", "H", HillElements)
    }
  } else if (carbonCheck) {
    if (labeledCarbonCheck) {
      x_CH <- c(x_c, x_12c, x_13c)
      carbonVector <- do.call(c, lapply(x_CH, function(i){
        HillElements[i]
      }))
      ##
      HillElements <- HillElements[-x_CH]
      HillElements <- c(carbonVector, HillElements)
    } else {
      x_CH <- c(x_c)
      HillElements <- HillElements[-x_CH]
      HillElements <- c("C", HillElements)
    }
  }
  ##
  ##############################################################################
  ##
  x_35cl <- which(HillElements == "[35]Cl")
  x_37cl <- which(HillElements == "[37]Cl")
  ##
  if ((length(x_35cl) == 1) | (length(x_37cl) == 1)) {
    x_cl <- which(HillElements == "Cl")
    x_LCl <- c(x_35cl, x_37cl)
    chlorineVector <- do.call(c, lapply(x_LCl, function(i){
      HillElements[i]
    }))
    ##
    HillElements <- HillElements[-x_LCl]
    HillElements <- c(HillElements[1:(x_cl - 1)], chlorineVector, HillElements[x_cl:(LElements - length(x_LCl))])
  }
  ##
  ##############################################################################
  ##
  x_LN <- which(HillElements == "[15]N")
  ##
  if (length(x_LN) == 1) {
    x_n <- which(HillElements == "N")
    HillElements <- HillElements[-x_LN]
    HillElements <- c(HillElements[1:(x_n - 1)], "[15]N", HillElements[x_n:(LElements - 1)])
  }
  ##
  ##############################################################################
  ##
  x_16o <- which(HillElements == "[16]O")
  x_17o <- which(HillElements == "[17]O")
  x_18o <- which(HillElements == "[18]O")
  ##
  if ((length(x_16o) == 1) | (length(x_17o) == 1) | (length(x_18o) == 1)) {
    x_o <- which(HillElements == "O")
    x_LO <- c(x_16o, x_17o, x_18o)
    oxygenVector <- do.call(c, lapply(x_LO, function(i){
      HillElements[i]
    }))
    ##
    HillElements <- HillElements[-x_LO]
    HillElements <- c(HillElements[1:(x_o - 1)], oxygenVector, HillElements[x_o:(LElements - length(x_LO))])
  }
  ##
  ##############################################################################
  ##
  x_32s <- which(HillElements == "[32]S")
  x_33s <- which(HillElements == "[33]S")
  x_34s <- which(HillElements == "[34]S")
  x_36s <- which(HillElements == "[36]S")
  ##
  if ((length(x_32s) == 1) | (length(x_33s) == 1) | (length(x_34s) == 1) | (length(x_36s) == 1)) {
    x_s <- which(HillElements == "S")
    x_LS <- c(x_32s, x_33s, x_34s, x_36s)
    sulfurVector <- do.call(c, lapply(x_LS, function(i){
      HillElements[i]
    }))
    ##
    HillElements <- HillElements[-x_LS]
    HillElements <- c(HillElements[1:(x_s - 1)], sulfurVector, HillElements[x_s:(LElements - length(x_LS))])
  }
  ##
  ##############################################################################
  ##
  orderHillElements <- do.call('c', lapply(HillElements, function(i) {
    which(Elements == i)
  }))
  ##
  ##############################################################################
  ##
  molecular_formula_printer <- function(orderHillElements, Elements, MolVec) {
    ##
    strMolecularFormula <- ""
    ##
    for (i in orderHillElements) {
      if (MolVec[i] > 0) {
        if (MolVec[i] == 1) {
          strMolecularFormula <- paste0(strMolecularFormula, Elements[i])
        } else {
          strMolecularFormula <- paste0(strMolecularFormula, Elements[i], MolVec[i])
        }
      }
    }
    ##
    if (strMolecularFormula == "") {
      strMolecularFormula <- NA
    }
    return(strMolecularFormula)
  }
  ##
  ##############################################################################
  ##
  MolVecMat <- matrix(MolVecMat, ncol = LElements)
  ##
  if (number_processing_threads == 1) {
    MolFormList <- do.call('c', lapply(1:dim(MolVecMat)[1], function(i) {
      molecular_formula_printer(orderHillElements, Elements, MolVecMat[i, ])
    }))
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    ############################################################################
    ##
    if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      MolFormList <-  foreach(i = 1:dim(MolVecMat)[1], .combine = 'c', .verbose = FALSE) %dopar% {
        molecular_formula_printer(orderHillElements, Elements, MolVecMat[i, ])
      }
      ##
      stopCluster(clust)
      ##
      ##########################################################################
      ##
    } else if (osType == "Linux") {
      ##
      MolFormList <- do.call('c', mclapply(1:dim(MolVecMat)[1], function(i) {
        molecular_formula_printer(orderHillElements, Elements, MolVecMat[i, ])
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
    }
  }
  ##
  ##############################################################################
  ##
  return(MolFormList)
}
