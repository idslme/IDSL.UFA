ionization_pathway_deconvoluter <- function(IonPathways, Elements) {
  ##
  IonPathways <- gsub(" ", "", IonPathways, fixed = TRUE)
  IonPathways <- gsub("*", "", IonPathways, fixed = TRUE)
  ##
  elementSortList <- element_sorter(ElementList = Elements, alphabeticalOrder = TRUE)
  ELementsDeconvolution <- elementSortList[["Elements"]]
  LElements <- length(ELementsDeconvolution)
  x_alpha <- do.call(c, lapply(Elements, function(i) {
    which(ELementsDeconvolution == i)
  }))
  ##
  ##############################################################################
  ##
  xH <- which(Elements == "H")
  xPosChargedElements <- c(xH, which((Elements == "Li") | (Elements == "Na") | (Elements == "K")))
  xNegChargedElements <- which((Elements == "Br") | (Elements == "Cl") | (Elements == "I"))
  xNeutralElements <- setdiff(seq(1, LElements, 1), c(xPosChargedElements, xNegChargedElements))
  ##
  xN <- which(Elements == "N")
  xNonChargedNH4 <- setdiff(seq(1, LElements, 1), c(xN, xH))
  ##
  functionalGroupChargeCalculator <- function(moleFVec, numberAdducts, xH, xN, xPosChargedElements, xNegChargedElements, xNeutralElements, xNonChargedNH4) {
    ##
    pwFunctionalGroupCharge <- 0
    ##
    if (length(which(moleFVec[xNeutralElements] > 0)) == 0) {
      if ((sum(moleFVec[xPosChargedElements]) == 1) & (sum(moleFVec[xNegChargedElements]) == 0)) {
        pwFunctionalGroupCharge <- +1 # H/Li/Na/K
      } else if ((sum(moleFVec[xPosChargedElements]) == 0) & (sum(moleFVec[xNegChargedElements]) == 1)) {
        if (numberAdducts == 2) {
          pwFunctionalGroupCharge <- -1 # Br/Cl/I
        }
      }
    } else if (length(which(moleFVec[xNonChargedNH4] > 0)) == 0) {
      if ((moleFVec[xN] == 1) & (moleFVec[xH] == 4)) { # NH4
        pwFunctionalGroupCharge <- +1
      }
    }
    ##
    return(pwFunctionalGroupCharge)
  }
  ##
  ##############################################################################
  ##
  IonPW_DC <- lapply(IonPathways, function(pw) {
    ##
    adductTypeList <- NULL
    ncharIPW <- nchar(pw)
    ##
    if (ncharIPW > 2) {
      coeff <- 0
      coeff2 <- substr(pw, 2, 2)
      ## coeff = 1
      if (coeff2 == "M") {
        coeff <- 1
      } else {## coeff > 1
        coeff3 <- substr(pw, 3, 3)
        if (coeff3 == "M") {
          coeff <- tryCatch(as.numeric(coeff2), warning = function(w) {0})
        }
      }
      ##
      if (coeff > 0) {
        x_bracket <- UFA_locate_regex(pw, "[]]")
        if (!is.null(x_bracket)) {
          x_bracket <- x_bracket[nrow(x_bracket), 1]
          ##
          xp <- UFA_locate_regex(pw, "[+]")
          if (!is.null(xp)) {
            xp <- xp[, 1]
          }
          ##
          xm <- UFA_locate_regex(pw, "-")
          if (!is.null(xm)) {
            xm <- xm[, 1]
          }
          ##
          xpm <- sort(unique(c(xp, xm, x_bracket)), decreasing = FALSE)
          xpm <- xpm[which(xpm <= x_bracket)]
          ##
          molfp <- rep(0, LElements)
          molfm <- molfp
          Lxpm <- length(xpm)
          ##
          pwCharge <- 0
          ##
          if (Lxpm > 1) {
            for (j in 1:(Lxpm - 1)) {
              if (substr(pw, xpm[j], xpm[j]) == "+") {
                molfp1 <- substr(pw, (xpm[j] + 1), (xpm[j + 1] - 1))
                ##
                if (grepl("[[:digit:]]", substr(molfp1, 1, 1))) {
                  coeffMolfp1 <- as.numeric(substr(molfp1, 1, 1))
                  molfp1 <- substr(molfp1, 2, nchar(molfp1))
                } else {
                  coeffMolfp1 <- 1
                }
                ##
                molfp1_vec <- formula_vector_generator(molfp1, Elements, LElements, allowedRedundantElements = TRUE)
                ##
                pwCharge <- pwCharge + coeffMolfp1*functionalGroupChargeCalculator(molfp1_vec, Lxpm, xH, xN, xPosChargedElements, xNegChargedElements, xNeutralElements, xNonChargedNH4)
                ##
                molfp <- molfp + coeffMolfp1*molfp1_vec
              } else if (substr(pw, xpm[j], xpm[j]) == "-") {
                molfm1 <- substr(pw, (xpm[j] + 1), (xpm[j + 1] - 1))
                ##
                if (grepl("[[:digit:]]", substr(molfm1, 1, 1))) {
                  coeffMolfm1 <- as.numeric(substr(molfm1, 1, 1))
                  molfm1 <- substr(molfm1, 2, nchar(molfm1))
                } else {
                  coeffMolfm1 <- 1
                }
                ##
                molfm1_vec <- formula_vector_generator(molfm1, Elements, LElements, allowedRedundantElements = TRUE)
                ##
                pwCharge <- pwCharge - coeffMolfm1*functionalGroupChargeCalculator(molfm1_vec, Lxpm, xH, xN, xPosChargedElements, xNegChargedElements, xNeutralElements, xNonChargedNH4)
                ##
                molfm <- molfm + coeffMolfm1*molfm1_vec
              }
            }
          }
          ##
          ######################################################################
          ##
          if (!is.infinite(molfp[1]) & !is.infinite(molfm[1])) {
            molfpm <- molfp - molfm
            ##
            if (pwCharge != 0) {
              signCharge <- if (pwCharge > 0) {"+"} else {"-"}
              absCharge <- abs(pwCharge)
              if (absCharge == 1) {
                chargeIPW <- signCharge
              } else {
                chargeIPW <- paste0(absCharge, signCharge)
              }
            } else {
              chargeIPW <- ""
            }
            ##
            if (ncharIPW >= (x_bracket + 1)) {
              originalChargeIPW <- substr(pw, (x_bracket + 1), ncharIPW)
              originalChargeIPW <- gsub(',| |[.]|;|"', '', originalChargeIPW)
              if (nchar(originalChargeIPW) >= nchar(chargeIPW)) {
                chargeIPW <- originalChargeIPW
              }
            }
            ##
            molfpm <- molfpm[x_alpha]
            ##
            adductTypeList <- list(coeff, molfpm, chargeIPW)
          }
        }
      }
    }
    ##
    return(adductTypeList)
  })
  ##
  names(IonPW_DC) <- IonPathways
  ##
  return(IonPW_DC)
}
