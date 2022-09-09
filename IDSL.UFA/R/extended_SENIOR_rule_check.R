extended_SENIOR_rule_check <- function(mol_vec, valence_vec, ionization_correction = 0) {
  rule2 <- FALSE
  sumVal <- ionization_correction + sum(valence_vec*mol_vec)
  if (sumVal %% 2 == 0) {
    max_non0_el <- max(valence_vec[mol_vec > 0])
    if (sumVal >= 2*max_non0_el) {
      if (sumVal >= (2*(sum(mol_vec) - 1))) {
        rule2 <- TRUE
      }
    }
  }
  return(rule2)
}
