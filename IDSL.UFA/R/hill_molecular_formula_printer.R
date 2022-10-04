hill_molecular_formula_printer <- function(Elements, MolVecMat, number_processing_threads = 1) {
  HillElements <- sort(Elements, decreasing = FALSE)
  x_c <- which(HillElements == "C")
  x_h <- which(HillElements == "H")
  if ((length(x_c) == 1) & (length(x_h) == 1)) {
    x_c_h <- c(x_c, x_h)
    HillElements <- HillElements[-x_c_h]
    HillElements <- c("C", "H", HillElements)
  } else if ((length(x_c) == 1) & (length(x_h) == 0)) {
    HillElements <- HillElements[-x_c]
    HillElements <- c("C", HillElements)
  }
  ##
  L_Hill_Elements <- length(HillElements)
  x_alpha_hill <- sapply(1:L_Hill_Elements, function(i) {
    which(Elements == HillElements[i])
  })
  ##
  MolVecMat <- matrix(MolVecMat, ncol = L_Hill_Elements)
  ##
  MolFormList_call <- function(k) {
    MolFor <- NULL
    MolVec <- MolVecMat[k, ]
    for (i in 1:L_Hill_Elements) {
      x_vec <- x_alpha_hill[i]
      if (MolVec[x_vec] > 0) {
        if (MolVec[x_vec] == 1) {
          MolFor <- paste0(MolFor, HillElements[i])
        } else {
          MolFor <- paste0(MolFor, HillElements[i], MolVec[x_vec])
        }
      }
    }
    if (is.null(MolFor)) {
      MolFor <- NA
    }
    return(MolFor)
  }
  ##
  if (number_processing_threads == 1) {
    MolFormList <- do.call(c, lapply(1:dim(MolVecMat)[1], function(k) {
      MolFormList_call(k)
    }))
  } else {
    osType <- Sys.info()[['sysname']]
    if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      MolFormList <-  foreach(k = 1:dim(MolVecMat)[1], .combine = 'c', .verbose = FALSE) %dopar% {
        MolFormList_call(k)
      }
      ##
      stopCluster(clust)
      ##
    } else if (osType == "Linux") {
      ##
      MolFormList <- do.call(c, mclapply(1:dim(MolVecMat)[1], function(k) {
        MolFormList_call(k)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
    }
  }
  return(MolFormList)
}
