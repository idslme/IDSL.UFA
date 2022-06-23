ionization_pathway_deconvoluter <- function(IonPathways, Elements) {
  ##
  IonPathways <- gsub(" ", "", IonPathways, fixed = TRUE)
  IonPathways <- gsub("*", "", IonPathways, fixed = TRUE)
  ##
  EL <- element_sorter(ElementList = Elements)
  EL_DC <- EL[[1]]
  L_Elements <- length(EL_DC)
  x_alpha <- sapply(1:L_Elements, function(i) {
    which(EL_DC == Elements[i])
  })
  ##
  IonPW_DC <- lapply(1:length(IonPathways), function(i) {
    b <- gsub(" ", "", IonPathways[i])
    xp <- UFA_locate_regex(b, "[+]")[, 1]
    xm <- UFA_locate_regex(b, "-")[, 1]
    x_bracket <- UFA_locate_regex(b, "[]]")[, 1]
    xpm <- sort(unique(c(xp, xm)))
    xpm <- xpm[xpm < x_bracket]
    xpm <- c(xpm, x_bracket)
    #
    coeff <- 1
    if (substr(b, 2, 2) != "M") {
      coeff <- as.numeric(substr(b, 2, 2))
    }
    #
    molfp <- rep(0, L_Elements)
    molfm <- rep(0, L_Elements)
    if (length(xpm) > 1) {
      for (j in 1:(length(xpm) - 1)) {
        if (substr(b, xpm[j], xpm[j]) == "+") {
          molfp1 <- substr(b, (xpm[j] + 1), (xpm[j + 1] - 1))
          molfp1_vec <- formula_vector_generator(molfp1, EL_DC, L_Elements)
          molfp <- molfp + molfp1_vec
        }
        if (substr(b, xpm[j], xpm[j]) == "-") {
          molfm1 <- substr(b, (xpm[j] + 1), (xpm[j + 1] - 1))
          molfm1_vec <- formula_vector_generator(molfm1, EL_DC, L_Elements)
          molfm <- molfm + molfm1_vec
        }
      }
    }
    molfpm <- molfp - molfm
    molfpm <- molfpm[x_alpha]
    list(coeff, molfpm)
  })
  return(IonPW_DC)
}
