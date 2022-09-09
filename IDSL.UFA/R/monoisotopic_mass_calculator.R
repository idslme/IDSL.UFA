monoisotopic_mass_calculator <- function(MoleFormVec, Elements_mass_abundance) {
  Non0Elements <- which(MoleFormVec != 0)
  ##
  if (length(Non0Elements) > 0) {
    MImass <- do.call(sum, lapply(Non0Elements, function(i) {
      IsotopeMonoMass <- Elements_mass_abundance[[i]][[1]][1]
      MoleFormVec[i]*IsotopeMonoMass
    }))
  } else {
    MImass <- 0
  }
  return(MImass)
}
