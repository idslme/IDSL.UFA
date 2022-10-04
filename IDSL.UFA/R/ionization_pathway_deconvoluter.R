ionization_pathway_deconvoluter <- function(IonPathways, Elements) {
  ##
  IonPathways <- gsub(" ", "", IonPathways, fixed = TRUE)
  IonPathways <- gsub("*", "", IonPathways, fixed = TRUE)
  ##
  EL_DC <- element_sorter(ElementList = Elements)[["Elements"]]
  L_Elements <- length(EL_DC)
  x_alpha <- sapply(1:L_Elements, function(i) {
    which(EL_DC == Elements[i])
  })
  ##
  IonPW_DC <- lapply(IonPathways, function(pw) {
    xp <- UFA_locate_regex(pw, "[+]")[, 1]
    xm <- UFA_locate_regex(pw, "-")[, 1]
    x_bracket <- UFA_locate_regex(pw, "[]]")[, 1]
    xpm <- sort(unique(c(xp, xm)))
    xpm <- xpm[xpm < x_bracket]
    xpm <- c(xpm, x_bracket)
    #
    coeff <- 1
    if (substr(pw, 2, 2) != "M") {
      coeff <- tryCatch(as.numeric(substr(pw, 2, 2)), warning = function(w){
        message(paste0("WARNING!!! Problem with ", pw))
        return(1)
      })
    }
    #
    molfp <- rep(0, L_Elements)
    molfm <- molfp
    L_xpm <- length(xpm)
    if (L_xpm > 1) {
      for (j in 1:(L_xpm - 1)) {
        if (substr(pw, xpm[j], xpm[j]) == "+") {
          molfp1 <- substr(pw, (xpm[j] + 1), (xpm[j + 1] - 1))
          molfp1_vec <- formula_vector_generator(molfp1, EL_DC, L_Elements)
          molfp <- molfp + molfp1_vec
        } else if (substr(pw, xpm[j], xpm[j]) == "-") {
          molfm1 <- substr(pw, (xpm[j] + 1), (xpm[j + 1] - 1))
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
