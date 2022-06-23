UFA_score_coefficient_workflow <- function(spreadsheet) {
  print("Initiated testing the UFA score coefficient spreadsheet consistency!")
  PARAM_SFT <- UFA_score_function_optimization_xlsxAnalyzer(spreadsheet)
  if (length(PARAM_SFT) > 1) {
    print("The UFA score coefficient spreadsheet is consistent with the score optimization workflow!")
    ##
    x0001 <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0001'), 2]
    if (tolower(x0001) == "yes") {
      zero_score_function(PARAM_SFT)
    }
    ##
    x0003 <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0003'), 2]
    if (tolower(x0003) == "yes") {
      score_coefficients_optimization(PARAM_SFT)
    }
    ##
    x0004 <- PARAM_SFT[which(PARAM_SFT[, 1] == 'SFT0004'), 2]
    if (tolower(x0004) == "yes") {
      score_coefficient_evaluation(PARAM_SFT)
    }
  }
}
