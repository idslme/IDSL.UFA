UFA_IPdbMerger <- function(path, vecIPDB) {
  ##
  MassMAIso <- NULL
  MolecularFormulaVec <- NULL
  IsotopicProfileList <- NULL
  IP_R13C <- NULL
  IndexMAIso <- NULL
  IPsize <- NULL
  ##
  for (i in vecIPDB) {
    IPDB <- IDSL.IPA::loadRdata(paste0(path, "/", i))
    ##
    MassMAIso <- c(MassMAIso, IPDB[["MassMAIso"]])
    MolecularFormulaVec <- c(MolecularFormulaVec, IPDB[["MolecularFormula"]])
    IsotopicProfileList <- c(IsotopicProfileList, IPDB[["IsotopicProfile"]])
    IP_R13C <- c(IP_R13C, IPDB[["R13C"]])
    IndexMAIso <- c(IndexMAIso, IPDB[["IndexMAIso"]])
    IPsize <- c(IPsize, IPDB[["IPsize"]])
  }
  IPDB <- NULL
  ##
  ##############################################################################
  ##
  AggregatedList <- aggregatedIPdbListGenerator(MassMAIso)
  ##
  ##############################################################################
  ##
  IPDB <- list(AggregatedList, MassMAIso, MolecularFormulaVec, IsotopicProfileList, IP_R13C, IndexMAIso, IPsize)
  names(IPDB) <- c("AggregatedList", "MassMAIso", "MolecularFormula", "IsotopicProfile", "R13C", "IndexMAIso", "IPsize")
  ##
  return(IPDB)
}
