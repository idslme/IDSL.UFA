# Multiplicative identification score
identification_score <- function(Score_coefficients, N_isotopologues, PCS, RCS, NEME, maxNEME, R13C_PL, R13C_IP) {
  S <- (N_isotopologues^Score_coefficients[1])*((PCS/1000)^Score_coefficients[2])*((RCS/100)^Score_coefficients[3])/((NEME/maxNEME)^Score_coefficients[4])/(exp(abs(log(R13C_PL/R13C_IP)))^Score_coefficients[5])
  return(S)
}
