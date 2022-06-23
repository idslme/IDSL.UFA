molecular_formula_library_generator <- function(entire_molecular_formulas) {
  entire_molecular_formulas <- sort(entire_molecular_formulas, decreasing = FALSE)
  A <- table(entire_molecular_formulas)
  MolFs <- names(A)
  x_UIC <- which(grepl("[[:punct:]]", MolFs) == FALSE)
  MolFs <- MolFs[x_UIC]
  Freq_Molfs <- as.numeric(A[x_UIC])
  molecular_formula_freq_database <- Freq_Molfs
  names(molecular_formula_freq_database) <- MolFs
  return(molecular_formula_freq_database)
}
