UFA_PubChem_formula_extraction <- function(path) {
  ##
  ##############################################################################
  ##
  if (!dir.exists(path)) {
    tryCatch(dir.create(path, recursive = TRUE), warning = function(w) {stop(paste0("Can't create `", path, "`!"))})
  }
  ##
  ##############################################################################
  ##
  IPA_message("This module may require a few hours to complete downloading and aggregating molecular formulas!", failedMessage = FALSE)
  sdfPubChemURL <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
  ##
  sdfPath <- paste0(path, "/SDF")
  if (!dir.exists(sdfPath)) {
    dir.create(sdfPath, recursive = TRUE)
  }
  ##
  sdfURLs <- paste0(sdfPath, "/sdfURLs.txt")
  utils::download.file(url = sdfPubChemURL, destfile = sdfURLs, quiet = TRUE, mode = "wb")
  sdfURLtxt <- readLines(sdfURLs)
  tryCatch(unlink(sdfURLs, recursive = TRUE), error = function(e) {warning(paste0("Can't delete `", sdfURLs, "`!"))})
  ##
  sdfList <- do.call(c, lapply(1:length(sdfURLtxt), function(i) {
    sdfURLloc <- UFA_locate_regex(sdfURLtxt[i], 'Compound_.*?.sdf.gz')
    if (!is.null(sdfURLloc)) {
      if (dim(sdfURLloc)[2] >= 2) {
        paste0(substr(sdfURLtxt[i], sdfURLloc[2, 1], sdfURLloc[2, 2]))
      }
    }
  }))
  sdfList <- unique(sdfList)
  LsdfList <- length(sdfList)
  ##
  progressBARboundaries <- txtProgressBar(min = 0, max = LsdfList, initial = 0, style = 3)
  i <- 0
  while (i < LsdfList) {
    i <- i + 1
    Sys.sleep(2) ## To wait for two seconds between downloads
    tryCatch(utils::download.file(url = paste0(sdfPubChemURL, sdfList[i]), destfile = paste0(sdfPath, "/", sdfList[i]), quiet = TRUE, mode = "wb"),
             error = function(e) {
               i <- i - 1
               Sys.sleep(180) ## To wait for three minutes and re-start downloading files
             })
    ##
    setTxtProgressBar(progressBARboundaries, i)
  }
  close(progressBARboundaries)
  ##
  IPA_message("Completed downloading `SDF` data!", failedMessage = FALSE)
  ##
  IPA_message("Initiated aggregating molecular formulas from `SDF` data!", failedMessage = FALSE)
  ##
  progressBARboundaries <- txtProgressBar(min = 0, max = LsdfList, initial = 0, style = 3)
  ##
  PubChemFormulas <- do.call(c, lapply(1:LsdfList, function(i) {
    setTxtProgressBar(progressBARboundaries, i)
    ##
    sdftext <- tryCatch(readLines(paste0(sdfPath, "/", sdfList[i]), warn = FALSE), error = function(e) {""})
    xmf <- grep("PUBCHEM_MOLECULAR_FORMULA", sdftext)
    if (length(xmf) > 0) {
      sdftext[xmf + 1]
    }
  }))
  close(progressBARboundaries)
  IPA_message("Completed aggregating molecular formulas from `SDF` data!", failedMessage = FALSE)
  ##
  PubChemFormulaFileName <- gsub("-", "_", paste0(path, "/PubChemFormula_", Sys.Date(), ".txt"), fixed = TRUE)
  write.table(PubChemFormulas, file = PubChemFormulaFileName, quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
  IPA_message(paste0("Stored molecular formulas as `", PubChemFormulaFileName,"` in the `", path, "` directory!"), failedMessage = FALSE)
  ##
  return()
}
