UFA_locate_regex <- function(string, pattern, ignore.case = FALSE, perl = FALSE,
                             fixed = FALSE, useBytes = FALSE) {
  ##
  UFAgreg <- gregexpr(pattern, string, ignore.case, perl, fixed, useBytes)[[1]]
  if (UFAgreg[1] > 0) {
    UFAgreg_lengthchar <- attributes(UFAgreg)$match.length
    #
    loc_mat <- do.call(rbind, lapply(1:length(UFAgreg_lengthchar), function(i) {
      c(UFAgreg[i], (UFAgreg[i] + UFAgreg_lengthchar[i] - 1))
    }))
  } else {
    loc_mat <- NULL
  }
  #
  return(loc_mat)
}
