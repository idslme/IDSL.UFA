molecular_formula_annotator <- function(IPDB, spectraList, peaklist, selectedIPApeaks, massAccuracy, maxNEME, minPCS,
                                        minNDCS, minRCS, scoreCoefficients, RTtolerance = NA, correctedRTpeaklist = NULL,
                                        exportSpectraParameters = NULL, number_processing_threads = 1) {
  ##
  MolecularFormulaAnnotationTable <- NULL
  ##
  ######################### retention time check ###############################
  ##
  if (is.na(RTtolerance)) {
    retentionTimeCheck <- FALSE
  } else {
    retentionTimeCheck <- TRUE
    if (is.null(correctedRTpeaklist)) {
      correctedRetentionTimeCheck <- FALSE
    } else {
      correctedRetentionTimeCheck <- TRUE
    }
  }
  ##
  ##############################################################################
  ##
  if (is.null(exportSpectraParameters)) {
    exportSpectraCheck <- FALSE
  } else {
    exportSpectraCheck <- TRUE
    maxAllowedNumberHits <- as.numeric(exportSpectraParameters[1])
    HRMSfileName <- exportSpectraParameters[2]
    msPolarity <- exportSpectraParameters[3]
    outputProfileSpectraHRMS <- paste0(exportSpectraParameters[4], "/", HRMSfileName)
    ##
    if (!dir.exists(outputProfileSpectraHRMS)) {
      dir.create(outputProfileSpectraHRMS, recursive = TRUE)
    }
    ##
    ############################################################################
    ##
    namesAnnotation <- c('IDSL.IPA_PeakID', 'IonFormulaID', 'isotopologueCount', 'IonFormula',
                         'mzTheoretical', 'mzPeaklist', 'retentionTime(min)', 'PeakHeight',
                         'normEucMassError(mDa)', 'Similaity(per-mille)', 'R13C Peaklist(%)',
                         'R13C Theoretical(%)', 'NDCS @80%', 'RatioChromScan(%) @80', 'rankOrder')
  }
  ##
  ##############################################################################
  ##
  call_molecular_formula_annotator <- function(i) {
    ##
    peakMolFormMatrix <- NULL
    ##
    mz12CTarget <- peaklist[i, 8]
    massPLround <- round(mz12CTarget, digits = 2)
    IDjIPDB <- do.call(c, lapply(c(-0.01, 0, 0.01), function(j) {
      IPDB[["AggregatedList"]][[as.character(massPLround + j)]]
    }))
    ##
    if (!is.null(IDjIPDB)) {
      rtTargetPeaklist <- peaklist[i, 3]
      ##########################################################################
      ######################## Retention Time Check ############################
      ##########################################################################
      if (retentionTimeCheck) {
        if (correctedRetentionTimeCheck) {
          xIDjIPDB <- which(abs(correctedRTpeaklist[i] - IPDB[["Retention Time"]][IDjIPDB]) <= RTtolerance)
        } else {
          xIDjIPDB <- which(abs(rtTargetPeaklist - IPDB[["Retention Time"]][IDjIPDB]) <= RTtolerance)
        }
        if (length(xIDjIPDB) > 0) {
          IDjIPDB <- IDjIPDB[xIDjIPDB]
        } else {
          IDjIPDB <- NULL
        }
      }
      ##########################################################################
      ##########################################################################
      ##########################################################################
      if (!is.null(IDjIPDB)) {
        xIDjIPDB <- which(abs(mz12CTarget - IPDB[["MassMAIso"]][IDjIPDB]) <= massAccuracy)
        if (length(xIDjIPDB) > 0) {
          IDjIPDB <- IDjIPDB[xIDjIPDB]
          ##
          int12CTarget <- peaklist[i, 4]
          PL_R13C <- peaklist[i, 11]
          scanNumberStart <- peaklist[i, 1]
          scanNumberEnd <- peaklist[i, 2]
          nScans <- scanNumberEnd - scanNumberStart + 1
          ##
          mzAnnotationEIClist <- lapply(IDjIPDB, function(IDj) {
            ##
            annotationEIClist <- NULL
            ##
            IsotopicProfile <- IPDB[["IsotopicProfile"]][[IDj]]
            IPsize <- IPDB[["IPsize"]][IDj]
            experimentalMZ <- matrix(rep(0, IPsize*nScans), ncol = nScans)
            experimentalINT <- experimentalMZ
            counterScan <- 0
            for (sc in scanNumberStart:scanNumberEnd) {
              counterScan <- counterScan + 1
              specScan <- spectraList[[sc]]
              for (Iso in 1:IPsize) {
                xIso <- which(abs(specScan[, 1] - IsotopicProfile[Iso, 1]) <= massAccuracy)
                LxIso <- length(xIso)
                if (LxIso > 0) {
                  if (LxIso > 1) {
                    xIsoMin <- which.min(abs(specScan[xIso, 1] - IsotopicProfile[Iso, 1]))
                    xIso <- xIso[xIsoMin[1]]
                  }
                  experimentalMZ[Iso, counterScan] <- specScan[xIso, 1]
                  experimentalINT[Iso, counterScan] <- specScan[xIso, 2]
                }
              }
            }
            integratedExperimentalINT <- rowSums(experimentalINT)
            xInt0 <- which(integratedExperimentalINT == 0)
            if (length(xInt0) == 0) {
              PCS <- sum(integratedExperimentalINT*IsotopicProfile[, 2])/sqrt(sum(integratedExperimentalINT^2)*sum(IsotopicProfile[, 2]^2))*1000 # in per-mille
              if (PCS >= minPCS) {
                integratedExperimentalMZ <- rowSums(experimentalMZ*experimentalINT)/integratedExperimentalINT
                NEME <- sqrt(sum((integratedExperimentalMZ - IsotopicProfile[, 1])^2)/IPsize)*1000 # in mDa
                if (NEME <= maxNEME) {
                  experimentalMZ[which(experimentalMZ > 0)] <- 1
                  nd <- colSums(experimentalMZ)
                  x100 <- IPDB[["IndexMAIso"]][IDj]
                  Int100 <- experimentalINT[x100, ]
                  intensityHeight <- max(Int100)
                  if (abs(intensityHeight - int12CTarget) < 2) { # To eliminate hits that do not match the IDSL.IPA peak height
                    x80 <- which(Int100/intensityHeight >= 0.2)
                    NDCS <- length(which(nd[x80] == IPsize))
                    if (NDCS >= minNDCS) {
                      L80 <- x80[length(x80)] - x80[1] + 1
                      RCS <- NDCS/L80*100
                      if (RCS >= minRCS) {
                        IP_R13C <- IPDB[["R13C"]][IDj]
                        IdentificationScore <- identificationScoreCalculator(scoreCoefficients, IPsize, PCS, RCS, NEME, PL_R13C, IP_R13C)
                        annotation <- c(i, IDj, IPsize,
                                        round(IsotopicProfile[x100, 1], digits = 5),
                                        round(mz12CTarget, digits = 5),
                                        round(rtTargetPeaklist, digits = 3),
                                        round(intensityHeight, digits = 0),
                                        round(NEME, digits = 2),
                                        round(PCS, digits = 3),
                                        round(PL_R13C, digits = 2),
                                        round(IP_R13C, digits = 2),
                                        NDCS,
                                        round(RCS, digits = 2),
                                        IdentificationScore)
                        ##
                        if (exportSpectraCheck) {
                          EIClist <- cbind(integratedExperimentalMZ, integratedExperimentalINT/integratedExperimentalINT[x100]*100)
                        } else {
                          EIClist <- NULL
                        }
                        ##
                        annotationEIClist <- list(annotation, EIClist)
                        ##
                      }
                    }
                  }
                }
              }
            }
            return(annotationEIClist)
          })
          ##
          xNonNULL <- which(do.call(c, lapply(mzAnnotationEIClist, function(j) {!is.null(j)})))
          numberMatchedMolForm <- length(xNonNULL)
          if (numberMatchedMolForm > 0) {
            peakMolFormMatrix <- do.call(rbind, lapply(xNonNULL, function(j) {
              mzAnnotationEIClist[[j]][[1]]
            }))
            ##
            if (numberMatchedMolForm > 1) {
              peakMolFormMatrix <- peakMolFormMatrix[order(peakMolFormMatrix[, 14], decreasing = TRUE), ]
            }
            ##
            molecularFormulaVec <- IPDB[["MolecularFormula"]][peakMolFormMatrix[, 2]]
            ##
            if (numberMatchedMolForm > 1) {
              peakMolFormMatrix <- cbind(peakMolFormMatrix[, 1:3], molecularFormulaVec, peakMolFormMatrix[, 4:13], seq(1, numberMatchedMolForm, 1))
            } else {
              peakMolFormMatrix <- matrix(c(peakMolFormMatrix[, 1:3], molecularFormulaVec, peakMolFormMatrix[, 4:13], 1), nrow = 1)
            }
            ##
            if (exportSpectraCheck) {
              names(mzAnnotationEIClist) <- as.character(IDjIPDB)
              ##
              ##################################################################
              ##
              maxNumberHits <- min(c(maxAllowedNumberHits, numberMatchedMolForm))
              ##
              for (R in 1:maxNumberHits) {
                ##
                PeakID <- peakMolFormMatrix[R, 1]
                characterIDj <- peakMolFormMatrix[R, 2]
                molecularFormulaIDj <- peakMolFormMatrix[R, 4]
                experimentalProfile <- mzAnnotationEIClist[[characterIDj]][[2]]
                ##
                IDj <- as.numeric(characterIDj)
                IsotopicProfile <- IPDB[["IsotopicProfile"]][[IDj]]
                IPsize <- IPDB[["IPsize"]][IDj]
                ##
                if (retentionTimeCheck) {
                  rtIPDB <- IPDB[["Retention Time"]][IDj]
                  spectraFilename <- paste0(outputProfileSpectraHRMS, "/", HRMSfileName, "_PeakID_", PeakID, "_", R, "_", molecularFormulaIDj, "_", round(rtIPDB, 3), "_.png")
                } else {
                  spectraFilename <- paste0(outputProfileSpectraHRMS, "/", HRMSfileName, "_PeakID_", PeakID, "_", R, "_", molecularFormulaIDj, "_.png")
                }
                fileCreateRCheck <- file.create(file = spectraFilename, showWarnings = FALSE)
                if (fileCreateRCheck) {
                  ##
                  annotationLabel <- do.call(c, lapply(1:15, function(j) {paste0(namesAnnotation[j], " = ", peakMolFormMatrix[R, j])}))
                  ##
                  lablelSpectra <- cbind(round(experimentalProfile[, 1], 5), do.call(c, lapply(1:IPsize, function(j) {max(c(experimentalProfile[j, 2], IsotopicProfile[j, 2]))[1] + 3.5})))
                  ##
                  png(spectraFilename, width = 20, height = 10, units = "in", res = 100)
                  ##
                  layout(matrix(c(1, 2), ncol = 2), widths = c(2, 1))
                  ##
                  plot(IsotopicProfile[, 1], IsotopicProfile[, 2], type = "h", xlim = c((min(IsotopicProfile[, 1]) - 1), (max(IsotopicProfile[, 1]) + 1)), ylim = c(0, 115),
                       lwd = 16, lend = 2, col = "blue", xlab = "", ylab = "", yaxt = "n", yaxs = "i")
                  ##
                  lines(experimentalProfile[, 1], experimentalProfile[, 2], type = "h", lwd = 4, lend = 2, col = "red")
                  ##
                  text(x = lablelSpectra[, 1], y = lablelSpectra[, 2], cex = 1.25, label = lablelSpectra[, 1])
                  text(x = (IsotopicProfile[1, 1] + IsotopicProfile[IPsize, 1])/2, y = 110, cex = 1.4, label = paste0("[", molecularFormulaIDj, "]", msPolarity))
                  ##
                  mtext(HRMSfileName, side = 3, adj = 0, line = 0.25, cex = 1.2)
                  mtext("m/z", side = 1, adj = 0.5, line = 2, cex = 1.35)
                  mtext("Intensity (%)", side = 2, adj = 0.5, line = 0.25, cex = 1.35)
                  mtext(text = paste0("+/- ", massAccuracy, " Da"), side = 3, adj = 1, line = 0.25, cex = 1.0)
                  ##
                  legend(x = "topright", legend = c("Theoretical", "Experimental"),
                         col = c("blue", "red"),
                         lwd = c(8, 4),
                         cex = 1.5, seg.len = 1, x.intersp = 0.5, y.intersp = 1)
                  ##
                  plot.new()
                  legend(x = "center", legend = annotationLabel,
                         cex = 1.40, bty = "n", x.intersp = 0.05, y.intersp = 1.3, seg.len = 0)
                  ##
                  dev.off()
                } else {
                  IPA_logRecorder(paste0("WARNING!!! Figure can not be created for `", molecularFormulaIDj, "` due to character length limit for `", HRMSfileName, "`!"))
                }
              }
            }
          }
        }
      }
    }
    ##
    return(peakMolFormMatrix)
  }
  ##
  if (length(selectedIPApeaks) > 0) {
    if (max(selectedIPApeaks) > dim(peaklist)[1]) {
      stop("Incorrect peak IDs from peaklist!")
    }
    ##
    ############################################################################
    ##
    if (number_processing_threads == 1) {
      MolecularFormulaAnnotationTable <- do.call(rbind, lapply(selectedIPApeaks, function(i) {
        call_molecular_formula_annotator(i)
      }))
    } else {
      osType <- Sys.info()[['sysname']]
      ##
      ##########################################################################
      ##
      if (osType == "Linux") {
        MolecularFormulaAnnotationTable <- do.call(rbind, mclapply(selectedIPApeaks, function(i) {
          call_molecular_formula_annotator(i)
        }, mc.cores = number_processing_threads))
        ##
        closeAllConnections()
        ##
        ########################################################################
        ##
      } else if (osType == "Windows") {
        clust <- makeCluster(number_processing_threads)
        registerDoParallel(clust)
        ##
        MolecularFormulaAnnotationTable <- foreach(i = selectedIPApeaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
          call_molecular_formula_annotator(i)
        }
        ##
        stopCluster(clust)
        ##
      }
    }
    ##
    ############################################################################
    ##
    if (!is.null(MolecularFormulaAnnotationTable)) {
      MolecularFormulaAnnotationTable <- data.frame(MolecularFormulaAnnotationTable) ## This must be a dataframe!
      rownames(MolecularFormulaAnnotationTable) <- NULL
      colnames(MolecularFormulaAnnotationTable) <- c('IDSL.IPA_PeakID', 'IonFormulaID', 'isotopologueCount', 'IonFormula', 'mzTheoretical', 'mzPeaklist', 'retentionTime',
                                                     'PeakHeight', 'normEucMassError', 'PCS', 'r13cPeaklist', 'r13cTheoretical', 'NDCS', 'ratioChromScan(%)', 'rank')
    }
  }
  ##
  return(MolecularFormulaAnnotationTable)
}
