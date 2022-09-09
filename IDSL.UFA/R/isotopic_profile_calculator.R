isotopic_profile_calculator <- function(MoleFormVec, Elements_mass_abundance, peak_spacing, intensity_cutoff,
                                        UFA_IP_memeory_variables = c(1e30, 1e-12, 10)) {
  ##############################################################################
  on.exit({
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  })
  #
  setTimeLimit(elapsed = UFA_IP_memeory_variables[3], transient = TRUE)
  ##############################################################################
  combination_formula <- function(n, k) {
    factorial(n + k - 1)/factorial(n - 1)/factorial(k)
  }
  Non0Elements <- which(MoleFormVec > 0)
  max_memeory_variables <- 1
  for (i in Non0Elements) {
    N_iso <- length(Elements_mass_abundance[[i]][[1]])
    Rep_ncr <- combination_formula(MoleFormVec[i], N_iso)
    max_memeory_variables <- max_memeory_variables * Rep_ncr
    if (is.nan(max_memeory_variables) | (max_memeory_variables > UFA_IP_memeory_variables[1])) {
      IsotopicProfile <- matrix(c(Inf, 100), ncol = 2)
      return(IsotopicProfile)
    }
  }
  ##
  nmultichoosek <- function(vec, k) {
    n <- length(vec)
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
    Am <- do.call(cbind, lapply(1:k, function(y) { vec[A[, y]]}))
    return(Am)
  }
  ##
  SUB_ISO_PRO <- function(NumberAtoms, IsotopeMasses, IsotopeCompositions) {
    IsotopeMasses_Combination <- nmultichoosek(IsotopeMasses, NumberAtoms)
    L_Ai <- dim(IsotopeMasses_Combination)[1]
    Ai <- rep(0, L_Ai)
    for (k in 1:L_Ai) {
      ElementCombination <- IsotopeMasses_Combination[k, ]
      FactValue <- 1
      R <- 1
      for (i in 1:length(IsotopeMasses)) {
        X <- length(which(ElementCombination == IsotopeMasses[i]))
        FactValue <- FactValue*factorial(X)
        R <- R*(IsotopeCompositions[i]^X)
      }
      Ai[k] <- factorial(NumberAtoms)/FactValue*R
    }
    SumIsotopeMasses_Combination <- rowSums(IsotopeMasses_Combination)
    X <- list(SumIsotopeMasses_Combination, Ai)
    return(X)
  }
  ## Subfunction for element and isotopic combinations
  SUB_COMB <- function(NumberAtoms, IsotopeMasses, IsotopeCompositions, UFA_IP_memeory_variables2) {
    if (length(IsotopeCompositions) == 1) {
      element_combination <- IsotopeMasses*NumberAtoms
      RAe_element <- IsotopeCompositions
    } else {
      X <- SUB_ISO_PRO(NumberAtoms, IsotopeMasses, IsotopeCompositions)
      element_combination <- X[[1]]
      RAe_element <- X[[2]]
      x_c <- which(RAe_element > UFA_IP_memeory_variables2)
      RAe_element <- RAe_element[x_c]
      element_combination <- element_combination[x_c]
    }
    Combo <- list (element_combination, RAe_element)
    return(Combo)
  }
  ##
  N_elements <- length(Non0Elements)
  El_Mass <- vector(mode = "list", N_elements)
  RAel_El <- El_Mass
  Combination_Size <- 1
  el <- 0
  for (i in Non0Elements) {
    IsotopeMasses <- Elements_mass_abundance[[i]][[1]]
    IsotopeCompositions <- Elements_mass_abundance[[i]][[2]]
    Element_Combo <- SUB_COMB(MoleFormVec[i], IsotopeMasses, IsotopeCompositions, UFA_IP_memeory_variables[2])
    el <- el + 1
    El_Mass[[el]] <- Element_Combo[[1]]
    RAel_El[[el]] <- Element_Combo[[2]]
    Combination_Size <- Combination_Size * length(El_Mass[[el]])
  }
  MW <- rep(0, Combination_Size)
  RA <- rep(1, Combination_Size) # RA indicates the abundance of isotopic combinations
  counter <- 0
  ##
  listIndex <- 1
  sum2 <- 0
  CounterIsotopicProfile <- 0
  RecursiveMass <- function(El_Mass, N_elements, listIndex, sum2) {
    if (listIndex == (N_elements + 1)) {
      CounterIsotopicProfile <<- CounterIsotopicProfile + 1 # A global variable
      MW[CounterIsotopicProfile] <<- sum2 # A global variable
    } else {
      for (c in El_Mass[[listIndex]]) {
        RecursiveMass(El_Mass, N_elements, listIndex + 1, sum2 + c)
      }
    }
  }
  RecursiveMass(El_Mass, N_elements, listIndex, sum2)
  #
  listIndex <- 1
  prod2 <- 1
  CounterIsotopicProfile <- 0
  RecursiveAbundance <- function(RAel_El, N_elements, listIndex, prod2) {
    if (listIndex == (N_elements + 1)) {
      CounterIsotopicProfile <<- CounterIsotopicProfile + 1 # A global variable
      RA[CounterIsotopicProfile] <<- prod2 # A global variable
    } else {
      for (c in RAel_El[[listIndex]]) {
        RecursiveAbundance(RAel_El, N_elements, listIndex + 1, prod2*c)
      }
    }
  }
  RecursiveAbundance(RAel_El, N_elements, listIndex, prod2)
  ##
  B <- cbind(MW, RA)
  L_B <- length(MW)
  if (L_B == 1) {
    MolWeight <- B[1]
    RA2 <- B[2]
  } else {
    B <- B[which(B[, 2] > 1e-12), ]
    B <- B[order(B[, 2], decreasing = TRUE), ]
    if (peak_spacing != 0) {
      MolWeight <- rep(0, L_B)
      RA2 <- MolWeight
      Counter <- 0
      for (i in 1:dim(B)[1]) {
        if (B[i, 1] != 0) {
          x <- which(abs(B[, 1] - B[i, 1]) <= peak_spacing)
          Counter <- Counter + 1
          RA2[Counter] <- sum(B[x, 2])
          MolWeight[Counter] <- sum(B[x, 1]*B[x, 2])/RA2[Counter]
          B[x, ] <- 0
        }
      }
      RA2 <- RA2[1:Counter]
      MolWeight <- MolWeight[1:Counter]
    } else {
      RA2 <- B[, 2]
      MolWeight <- B[, 1]
    }
  }
  Intensity <- RA2/max(RA2)*100 # intensity or relative abundance
  IsotopicProfile <- matrix(cbind(MolWeight, Intensity), ncol = 2)
  IsotopicProfile <- matrix(IsotopicProfile[order(IsotopicProfile[, 1], decreasing = FALSE), ], ncol = 2)
  INDEX <- which(IsotopicProfile[, 2] >= intensity_cutoff)
  IsotopicProfile <- matrix(IsotopicProfile[INDEX, ], ncol = 2)
  return(IsotopicProfile)
}
