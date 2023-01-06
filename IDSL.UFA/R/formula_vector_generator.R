formula_vector_generator <- function(molecular_formula, Elements, LElements = length(Elements), allowedRedundantElements = FALSE) {
  ##
  elementsMatchCheck <- FALSE
  i <- 1
  MolecularFormulaVector <- rep(0, LElements)
  while (i <= LElements) {
    Z <- UFA_locate_regex(molecular_formula, Elements[i], fixed = TRUE)
    if (!is.null(Z)) {
      nZ <- dim(Z)[1]
      if (!allowedRedundantElements) {
        if (nZ > 1) {
          break
        }
      }
      for (iz in 1:nZ) {
        nextCharMF <- substr(molecular_formula, (Z[iz, 2] + 1), (Z[iz, 2] + 1))
        nextCharMFcheck <- grepl("[[:digit:]]", nextCharMF)
        if (nextCharMFcheck) {
          Number <- ""
          while (nextCharMFcheck) {
            Number <- paste0(Number, nextCharMF)
            Z[iz, 2] <- Z[iz, 2] + 1
            nextCharMF <- substr(molecular_formula, (Z[iz, 2] + 1), (Z[iz, 2] + 1))
            nextCharMFcheck <- grepl("[[:digit:]]", nextCharMF)
          }
          Number <- as.numeric(Number)
        } else {
          Number <- 1
        }
        MolecularFormulaVector[i] <- MolecularFormulaVector[i] + Number
      }
      ##
      Z <- do.call(c, lapply(1:nZ, function(j) {
        seq(Z[j, 1], Z[j, 2], 1)
      }))
      zPrime <- setdiff(seq(1, nchar(molecular_formula), 1), Z)
      if (length(zPrime) > 0) {
        molecular_formula <- do.call(paste0, lapply(zPrime, function(j) {
          substr(molecular_formula, j, j)
        }))
      } else {
        elementsMatchCheck <- TRUE
        i <- LElements
      }
    }
    i <- i + 1
  }
  if (!elementsMatchCheck) {
    MolecularFormulaVector <- rep(-Inf, LElements)
  }
  return(MolecularFormulaVector)
}
