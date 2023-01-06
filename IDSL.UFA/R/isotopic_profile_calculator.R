isotopic_profile_calculator <- function(MoleFormVec, massAbundanceList, peak_spacing, intensity_cutoff,
                                        UFA_IP_memeory_variables = c(1e30, 1e-12, 100)) {
  ##############################################################################
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))
  ##
  setTimeLimit(cpu = Inf, elapsed = UFA_IP_memeory_variables[3], transient = TRUE)
  ##############################################################################
  combination_formula <- function(n, k) {
    exp(lfactorial(n + k - 1) - lfactorial(k))/factorial(n - 1)
  }
  ##
  non0Elements <- which(MoleFormVec > 0)
  max_memeory_variables <- 1
  for (i in non0Elements) {
    nstableIsotopes <- massAbundanceList[[i]][[3]]
    if (nstableIsotopes > 2) {
      Rep_nCr <- combination_formula(nstableIsotopes, MoleFormVec[i])
      max_memeory_variables <- max_memeory_variables * Rep_nCr
      if (max_memeory_variables > UFA_IP_memeory_variables[1]) {
        IsotopicProfile <- matrix(c(Inf, 100), ncol = 2)
        return(IsotopicProfile)
      }
    }
  }
  ##
  ##############################################################################
  ##
  nmultichoosek <- function(n, k) {
    if (n == 1) {
      A <- matrix(rep(k, n), ncol = k)
    } else if (n == 2) {
      A <- matrix(rep(1, k*(k + 1)), ncol = k)
      for (j in 1:k) {
        for (i in (j + 1):(k + 1)) {
          A[i, j] <- 2
        }
      }
    } else {
      v <- 1:n
      ##########################################################################
      ## CREDIT: This  function was derived from the `combinations` function of the `gtools` package
      combinations_repetition <- function(n, k, v) {
        if (k == 1) {
          cr <- matrix(v, n, 1)
        } else if (n == 1) {
          cr <- matrix(v, 1, k)
        } else {
          cr <- rbind(cbind(v[1], Recall(n, k - 1, v)), Recall(n - 1, k, v[-1]))
        }
        return(cr)
      }
      ##########################################################################
      A <- combinations_repetition(n, k, v)
    }
    return(A)
  }
  ##
  ##############################################################################
  ##
  SUB_ISO_PRO <- function(iMassAbundanceList, numberAtoms) {
    ##
    isotopicCombinations <- nmultichoosek(iMassAbundanceList[[3]], numberAtoms)
    nAk <- dim(isotopicCombinations)[1]
    Ak <- rep(0, nAk)
    Mk <- Ak
    for (k in 1:nAk) {
      lfactValue <- 0
      R <- 1
      ##
      stableIsotopesTable <- table(isotopicCombinations[k, ])
      satbleIsotopeNames <- names(stableIsotopesTable)
      ##
      for (i in satbleIsotopeNames) {
        X <- stableIsotopesTable[[i]]
        ##
        Mk[k] <- Mk[k] + X*iMassAbundanceList[[1]][[i]]
        ##
        lfactValue <- lfactValue + lfactorial(X) ## `lfactorial`is a base R function to calculate logarithmic factorial values for big numbers.
        R <- R*(iMassAbundanceList[[2]][[i]]^X)
      }
      Ak[k] <- R*exp(lfactorial(numberAtoms) - lfactValue)
    }
    ##
    listMkAk <- list(Mk, Ak)
    ##
    return(listMkAk)
  }
  ##
  ##############################################################################
  ## Subfunction for element and isotopic combinations
  SUB_COMB <- function(iMassAbundanceList, numberAtoms, UFA_IP_memeory_variables2) {
    ##
    if (iMassAbundanceList[[3]] == 1) {
      element_combination <- iMassAbundanceList[[1]]*numberAtoms
      RAe_element <- iMassAbundanceList[[2]]
    } else {
      listMkAk <- SUB_ISO_PRO(iMassAbundanceList, numberAtoms)
      element_combination <- listMkAk[[1]]
      RAe_element <- listMkAk[[2]]
      xRAe <- which(RAe_element > UFA_IP_memeory_variables2)
      element_combination <- element_combination[xRAe]
      RAe_element <- RAe_element[xRAe]
    }
    ##
    Combo <- list(element_combination, RAe_element)
    ##
    return(Combo)
  }
  ##
  ##############################################################################
  ##
  nElements <- length(non0Elements)
  massCombinations <- vector(mode = "list", nElements)
  abundanceCombinations <- massCombinations
  nIsotopicCombinations <- 1
  el <- 0
  for (i in non0Elements) {
    elementalCombinations <- SUB_COMB(massAbundanceList[[i]], MoleFormVec[i], UFA_IP_memeory_variables[2])
    el <- el + 1
    massCombinations[[el]] <- elementalCombinations[[1]]
    abundanceCombinations[[el]] <- elementalCombinations[[2]]
    nIsotopicCombinations <- nIsotopicCombinations * length(massCombinations[[el]])
  }
  ##
  ##############################################################################
  ##
  MZ <- rep(0, nIsotopicCombinations)
  listIndex <- 1
  sum2 <- 0
  counterIsotopicCombination <- 0
  RecursiveMass <- function(massCombinations, nElements, listIndex, sum2) {
    if (listIndex > nElements) {
      counterIsotopicCombination <<- counterIsotopicCombination + 1       # A global variable
      MZ[counterIsotopicCombination] <<- sum2                             # A global variable
    } else {
      for (c in massCombinations[[listIndex]]) {
        RecursiveMass(massCombinations, nElements, listIndex + 1, sum2 + c)
      }
    }
  }
  RecursiveMass(massCombinations, nElements, listIndex, sum2)
  ##
  ##############################################################################
  ##
  if (nIsotopicCombinations > 1) {
    ##
    RA <- rep(1, nIsotopicCombinations) # RA indicates the abundance of isotopic combinations
    listIndex <- 1
    prod2 <- 1
    counterIsotopicCombination <- 0
    RecursiveAbundance <- function(abundanceCombinations, nElements, listIndex, prod2) {
      if (listIndex > nElements) {
        counterIsotopicCombination <<- counterIsotopicCombination + 1     # A global variable
        RA[counterIsotopicCombination] <<- prod2                          # A global variable
      } else {
        for (c in abundanceCombinations[[listIndex]]) {
          RecursiveAbundance(abundanceCombinations, nElements, listIndex + 1, prod2*c)
        }
      }
    }
    RecursiveAbundance(abundanceCombinations, nElements, listIndex, prod2)
    ##
    ############################################################################
    ##
    xRAe12 <- which(RA > UFA_IP_memeory_variables[2])
    MZ <- MZ[xRAe12]
    RA <- RA[xRAe12]
    nIsotopicCombinations <- length(xRAe12)
    ##
    orderMZ <- order(MZ, decreasing = FALSE)
    MZ <- MZ[orderMZ]
    RA <- RA[orderMZ]
    if (peak_spacing > 0) {
      B <- cbind(MZ, RA)
      xDIffMZ <- c(0, which(diff(B[, 1]) > peak_spacing), nIsotopicCombinations)
      LxDIffMZ <- length(xDIffMZ) - 1
      ##
      MZ <- rep(0, nIsotopicCombinations)
      RA <- MZ
      nB <- 0
      ##
      for (q in 1:LxDIffMZ) {
        xQ <- seq((xDIffMZ[q] + 1), xDIffMZ[q + 1], 1)
        ##
        nQ <- xDIffMZ[q + 1] - xDIffMZ[q]
        if (nQ == 1) {
          nB <- nB + 1
          RA[nB] <- B[xQ, 2]
          MZ[nB] <- B[xQ, 1]
          B[xQ, ] <- 0
        } else {
          xQ <- xQ[order(B[xQ, 2], decreasing = TRUE)]
          ##
          for (i in xQ) {
            if (B[i, 1] != 0) {
              x <- which(abs(B[xQ, 1] - B[i, 1]) <= peak_spacing)
              x <- xQ[x]
              nB <- nB + 1
              if (length(x) == 1) {
                RA[nB] <- B[x, 2]
                MZ[nB] <- B[x, 1]
              } else {
                RA[nB] <- sum(B[x, 2])
                MZ[nB] <- sum(B[x, 1]*B[x, 2])/RA[nB]
              }
              ##
              B[x, ] <- 0
            }
          }
        }
      }
      RA <- RA[1:nB]
      MZ <- MZ[1:nB]
    }
    ##
    Intensity <- RA/max(RA)*100 # intensity or relative abundance
    ##
    if (intensity_cutoff > 0) {
      INDEX <- which(Intensity >= intensity_cutoff)
      MZ <- MZ[INDEX]
      Intensity <- Intensity[INDEX]
    }
  } else {
    Intensity <- 100
  }
  ##
  ##############################################################################
  ##
  IsotopicProfile <- matrix(cbind(MZ, Intensity), ncol = 2)
  ##
  return(IsotopicProfile)
}
