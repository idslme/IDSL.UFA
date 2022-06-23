formula_vector_generator <- function(molecular_formula, Elements, L_Elements = length(Elements)) {
  MolecularFormulaVector <- rep(0, L_Elements)
  for (i in 1:L_Elements) {
    Z <- UFA_locate_regex(molecular_formula, Elements[i])
    if (!is.na(Z[1, 1])) {
      if (dim(Z)[1] > 1) {
        return(rep(-Inf, L_Elements))
      }
      z <- c(Z[1, 1], Z[1, 2])
      Charmolecular_formula <- substr(molecular_formula, z[2] + 1, z[2] + 1)
      if (!grepl("[[:digit:]]", Charmolecular_formula)) {
        Number <- 1
      } else {
        Number <- ""
        NextChar <- TRUE
        while (NextChar) {
          Charmolecular_formula <- substr(molecular_formula, z[2] + 1, z[2] + 1)
          if (grepl("[[:digit:]]", Charmolecular_formula)) {
            Number <- paste0(Number, Charmolecular_formula)
            z[2] <- z[2] + 1
          } else {
            NextChar <- FALSE
            Number <- as.numeric(Number)
          }
        }
      }
      MolecularFormulaVector[i] <- Number
      molecular_formula <- paste0(substr(molecular_formula, 1, z[1] - 1), substr(molecular_formula, z[2] + 1, nchar(molecular_formula)))
    }
  }
  if (molecular_formula != "") {
    MolecularFormulaVector <- rep(-Inf, L_Elements)
  }
  return(MolecularFormulaVector)
}
