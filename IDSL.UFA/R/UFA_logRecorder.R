UFA_logRecorder <- function(messageQuote, printMessage = TRUE) {
  if (printMessage) {
    print(messageQuote)
  }
  ##
  if (exists('logUFA')) {
    if (typeof(messageQuote) == "list") {
      namesMessageQuote <- names(messageQuote)
      for (i in 1:length(messageQuote)) {
        write(paste0(i, ": In ", paste0(deparse(messageQuote[[i]]), collapse = "\n")), file = logUFA, append = TRUE, sep = "\n")
        write(paste0("  ", namesMessageQuote[i]), file = logUFA, append = TRUE, sep = "\n")
      }
      ##
    } else {
      write(messageQuote, file = logUFA, append = TRUE, sep = "\n")
    }
  }
}
