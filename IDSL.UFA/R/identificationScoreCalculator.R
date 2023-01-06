# Multiplicative identification score
identificationScoreCalculator <- function(scoreCoefficients, nIisotopologues, PCS, RCS, NEME, R13C_PL, R13C_IP) {
  S <- (nIisotopologues^scoreCoefficients[1])*((PCS/1000)^scoreCoefficients[2])*((RCS/100)^scoreCoefficients[3])/(NEME^scoreCoefficients[4])/(exp(abs(log(R13C_PL/R13C_IP)))^scoreCoefficients[5])
  return(S)
}
