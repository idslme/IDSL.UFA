monoisotopicMassCalculator <- function(MoleFormVec, massAbundanceList, LElements = length(massAbundanceList)) {
  monoIsoMass <- do.call(sum, lapply(1:LElements, function(i) {
    MoleFormVec[i]*massAbundanceList[[i]][[1]][['1']]
  }))
  ##
  return(monoIsoMass)
}
