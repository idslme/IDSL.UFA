aggregatedIPdbListGenerator <- function(MassMAIso) {
  LMassMAIso <- length(MassMAIso)
  if (LMassMAIso > 1) {
    IDroundMass <- cbind(round(MassMAIso, digits = 2), seq(1, LMassMAIso, 1))
    IDroundMass <- IDroundMass[order(IDroundMass[, 1], decreasing = FALSE), ]
    xDiff <- c(0, which(abs(diff(IDroundMass[, 1])) > 0), LMassMAIso)
    LxDiff <- length(xDiff) - 1
    AggregatedList <- lapply(1:LxDiff, function(i) {
      IDroundMass[(xDiff[i] + 1):xDiff[i + 1], 2]
    })
    names(AggregatedList) <- as.character(IDroundMass[(xDiff[1:LxDiff] + 1), 1])
  } else {
    AggregatedList <- list(1)
    names(AggregatedList) <- as.character(round(MassMAIso, digits = 2))
  }
  return(AggregatedList)
}
