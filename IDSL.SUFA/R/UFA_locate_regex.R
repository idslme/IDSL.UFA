UFA_locate_regex <- function(string, pattern) {
  UFAgreg <- gregexpr(pattern, string)[[1]]
  if (UFAgreg[1] > 0) {
    UFAgreg_lengthchar <- attributes(UFAgreg)$match.length
    #
    loc_mat <- do.call(rbind, lapply(1:length(UFAgreg_lengthchar), function(i) {
      c(UFAgreg[i], (UFAgreg[i] + UFAgreg_lengthchar[i] - 1))
    }))
  } else {
    loc_mat <- matrix(c(NA, NA), nrow = 1)
  }
  #
  return(loc_mat)
}
