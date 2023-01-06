element_sorter <- function(ElementList = "all", alphabeticalOrder = TRUE) {
  ##
  if (!is.null(ElementList) & (ElementList[1] != "")) {
    IUPAC_Isotopes <- NULL
    load(paste0(system.file('extdata', package = 'IDSL.UFA'), "/IUPAC_Isotopes.rda"))
    ##
    if (tolower(ElementList[1]) == "all") {
      Elements <- unique(IUPAC_Isotopes$element)
    } else {
      Elements <- unique(ElementList)
    }
    ##
    if (alphabeticalOrder) {
      Elements <- sort(Elements, decreasing = FALSE)
      nCharElements <- do.call(c, lapply(Elements, nchar))
      Elements <- Elements[order(nCharElements, decreasing = TRUE)]
    }
    ##
    LElements <- length(Elements)
    massAbundanceList <- vector(mode = "list", LElements)
    valence <- rep(0, LElements)
    listIndexElements <- base::tapply(seq(1, length(IUPAC_Isotopes$element), 1), IUPAC_Isotopes$element, 'c', simplify = FALSE)
    ##
    for (i in 1:LElements) {
      x_el <- listIndexElements[[Elements[i]]]
      nStableIsotopes <- length(x_el)
      if (nStableIsotopes > 0) {
        seqSatbleIsotopes <- seq(1, nStableIsotopes, 1)
        ##
        massStableIsotopes <- IUPAC_Isotopes$mass[x_el]
        names(massStableIsotopes) <- seqSatbleIsotopes
        ##
        abundanceStableIsotopes <- IUPAC_Isotopes$abundance[x_el]
        names(abundanceStableIsotopes) <- seqSatbleIsotopes
        ##
        massAbundanceList[[i]] <- list(massStableIsotopes, abundanceStableIsotopes, nStableIsotopes)
        valence[i] <- IUPAC_Isotopes$valence[x_el[1]]
      } else {
        massAbundanceList[[i]] <- list(0, 1, 0)
      }
    }
    ##
    names(massAbundanceList) <- Elements
    ##
  } else {
    Elements <- NULL
    massAbundanceList <- NULL
    valence <- NULL
  }
  ##
  ElementList <- list(Elements, massAbundanceList, valence)
  names(ElementList) <- c("Elements", "massAbundanceList", "Valence")
  ##
  return(ElementList)
}
