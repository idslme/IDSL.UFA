element_sorter <- function(ElementList = "all", ElementOrder = "alphabetical") {
  ##
  IDSL.UFA::IUPAC_Isotopes
  ##
  NonUniElements <- IUPAC_Isotopes$element
  Valence_EL <- IUPAC_Isotopes$valence
  if (tolower(ElementList[1]) == "all") {
    Elements <- unique(NonUniElements)
  } else {
    Elements <- unique(ElementList)
  }
  ##
  if (tolower(ElementOrder) == "alphabetical") {
    Elements <- sort(Elements, decreasing = FALSE)
    X <- do.call(c, lapply(Elements, nchar))
    Elements <- cbind(Elements, X)
    Elements <- Elements[order(Elements[, 2], decreasing = TRUE), ]
    Elements <- Elements[, 1]
  }
  ##
  L_Elements <- length(Elements)
  Elements_mass_abundance <- vector(mode = "list", L_Elements)
  valence <- rep(0, L_Elements)
  for (i in 1:L_Elements) {
    x_el <- which(NonUniElements == Elements[i])
    Elements_mass_abundance[[i]] <- list(IUPAC_Isotopes$mass[x_el], IUPAC_Isotopes$abundance[x_el])
    valence[i] <- Valence_EL[x_el[1]]
  }
  ##
  names(Elements_mass_abundance) <- Elements
  ##
  ElementList <- list(Elements, Elements_mass_abundance, valence)
  names(ElementList) <- c("Elements", "massAbundanceList", "Valence")
  return(ElementList)
}
