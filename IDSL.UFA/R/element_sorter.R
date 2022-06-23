element_sorter <- function(ElementList = "all", ElementOrder = "alphabetical") {
  IDSL.UFA::IUPAC_Isotopes
  NonUniElements <- IUPAC_Isotopes$element
  Valence_EL <- IUPAC_Isotopes$valence
  if (tolower(ElementList[1]) == "all") {
    Elements <- unique(NonUniElements)
  } else {
    Elements <- unique(ElementList)
  }
  L_Elements <- length(Elements)
  if (tolower(ElementOrder) == "alphabetical") {
    Elements <- sort(Elements, decreasing = FALSE)
    X <- sapply(1:L_Elements, function(i) {nchar(Elements[i])})
    Elements <- cbind(Elements, X)
    Elements <- Elements[order(Elements[, 2], decreasing = TRUE), ]
    OutputElements <- Elements[, 1]
  } else if(tolower(ElementOrder) == "same") {
    OutputElements <- Elements
  }
  Elements_mass_abundance <- vector(mode = "list", L_Elements)
  valence <- rep(0, L_Elements)
  for (i in 1:L_Elements) {
    x_el <- which(NonUniElements == Elements[i])
    Elements_mass_abundance[[i]] <- list(IUPAC_Isotopes$mass[x_el], IUPAC_Isotopes$abundance[x_el])
    valence[i] <- Valence_EL[x_el[1]]
  }
  ElementList <- list(OutputElements, Elements_mass_abundance, valence)
  return(ElementList)
}
