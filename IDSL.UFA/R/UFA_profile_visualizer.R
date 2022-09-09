UFA_profile_visualizer <- function(PARAM_SA) {
  print("Initiated producing mass spectra!")
  input_path_hrms <- PARAM_SA[which(PARAM_SA[, 1] == 'PARAM0010'), 2]
  if (tolower(PARAM_SA[which(PARAM_SA[, 1] == 'PARAM0011'), 2]) == "all") {
    file_name_hrms <- dir(path = input_path_hrms)
    file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM_SA[which(PARAM_SA[, 1] == 'PARAM0012'), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
  } else {
    samples_string <- PARAM_SA[which(PARAM_SA[, 1] == 'PARAM0011'), 2]
    file_name_hrms <- strsplit(samples_string, ";")[[1]]
  }
  ##
  input_path_pl <- PARAM_SA[which(PARAM_SA[, 1] == 'PARAM0013'), 2]
  file_names_peaklist1 <- dir(path = input_path_pl, pattern = ".Rdata")
  file_names_peaklist2 <- dir(path = input_path_pl, pattern = "peaklist_")
  file_names_peaklist <- file_names_peaklist1[file_names_peaklist1 %in% file_names_peaklist2]
  #
  file_names_peaklist_hrms1 <- gsub(".Rdata", "", file_names_peaklist)
  file_names_peaklist_hrms2 <- gsub("peaklist_", "", file_names_peaklist_hrms1)
  file_names_peaklist_hrms <- file_name_hrms %in% file_names_peaklist_hrms2
  L_PL <- length(which(file_names_peaklist_hrms == TRUE))
  if (length(file_name_hrms) != L_PL) {
    stop("Error!!! peaklist files are not available for the entire selected HRMS files!")
  }
  ##
  output_path <- PARAM_SA[which(PARAM_SA[, 1] == 'PARAM0014'), 2]
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
    print("Created output directory!")
  }
  ##
  output_path_spectra <- paste0(output_path, "/UFA_spectra/")
  if (!dir.exists(output_path_spectra)) {
    dir.create(output_path_spectra, recursive = TRUE)
    print("Created UFA_spectra directory!")
  }
  ##
  molecular_formula <- eval(parse(text = paste0("c(", PARAM_SA[which(PARAM_SA[, 1] == 'SA0001'), 2], ")")))
  RT_target <- eval(parse(text = paste0("c(", PARAM_SA[which(PARAM_SA[, 1] == 'SA0002'), 2], ")")))
  delta_rt <- as.numeric(PARAM_SA[which(PARAM_SA[, 1] == 'SA0003'), 2])
  IonPathways <- eval(parse(text = paste0("c(", PARAM_SA[which(PARAM_SA[, 1] == 'SA0004'), 2], ")")))
  peak_spacing <- as.numeric(PARAM_SA[which(PARAM_SA[, 1] == 'SA0005'), 2])
  intensity_cutoff_str <- PARAM_SA[which(PARAM_SA[, 1] == 'SA0006'), 2]
  UFA_IP_memeory_variables <- eval(parse(text = paste0("c(", PARAM_SA[which(PARAM_SA[, 1] == "SA0007"), 2], ")")))
  mass_accuracy <- as.numeric(PARAM_SA[which(PARAM_SA[, 1] == 'SA0008'), 2])
  number_processing_threads <- as.numeric(PARAM_SA[which(PARAM_SA[, 1] == 'SA0009'), 2])
  exportSpectra <- if (tolower(PARAM_SA[which(PARAM_SA[, 1] == 'SA0010'), 2]) == "yes") {TRUE} else {FALSE}
  exportedAnnotatedSpectraTable <- if (tolower(PARAM_SA[which(PARAM_SA[, 1] == 'SA0011'), 2]) == "yes") {TRUE} else {FALSE}
  if (exportSpectra == TRUE | exportedAnnotatedSpectraTable == TRUE) {
    ##
    EL <- element_sorter()
    Elements <- EL[[1]]
    Elements_mass_abundance <- EL[[2]]
    L_Elements <- length(Elements)
    ##
    x_el_c <- which(Elements == "C")
    x_el_b <- which(Elements == "B")
    x_el_br <- which(Elements == "Br")
    x_el_cl <- which(Elements == "Cl")
    x_el_k <- which(Elements == "K")
    x_el_s <- which(Elements == "S")
    x_el_se <- which(Elements == "Se")
    x_el_si <- which(Elements == "Si")
    ##
    IonPW_DC <- ionization_pathway_deconvoluter(IonPathways, Elements)
    ##
    L_MolF <- length(molecular_formula)
    RT_target_ion <- NULL
    MoleFormVecMat <- do.call(rbind, lapply(1:L_MolF, function(i_molf) {
      FormulaVector <- formula_vector_generator(molecular_formula[i_molf], Elements, L_Elements)
      rt1 <- RT_target[i_molf]
      molf_deconvoluter_ipw <- do.call(rbind, lapply(IonPW_DC, function(IonPW) {
        molv_ipw <- NULL
        Ion_coeff <- IonPW[[1]]
        Ion_adduct <- IonPW[[2]]
        MoleFormVec <- Ion_coeff*FormulaVector + Ion_adduct
        x_neg <- which(MoleFormVec < 0)
        if (length(x_neg) == 0) {
          RT_target_ion <<- c(RT_target_ion, rt1)
          molv_ipw <- MoleFormVec
        }
        molv_ipw
      }))
      molf_deconvoluter_ipw
    }))
    ##
    L_MoleFormVecMat <- dim(MoleFormVecMat)[1]
    molecular_formula_hill <- hill_molecular_formula_printer(Elements, MoleFormVecMat, number_processing_threads)
    ##
    IP_calculator <- "IP_calculator <- function(i_mat) {
      ##
      c <- MoleFormVecMat[i_mat, x_el_c]
      b <- MoleFormVecMat[i_mat, x_el_b]
      br <- MoleFormVecMat[i_mat, x_el_br]
      cl <- MoleFormVecMat[i_mat, x_el_cl]
      k <- MoleFormVecMat[i_mat, x_el_k]
      s <- MoleFormVecMat[i_mat, x_el_s]
      se <- MoleFormVecMat[i_mat, x_el_se]
      si <- MoleFormVecMat[i_mat, x_el_si]
      ##
      intensity_cutoff <- intensity_cutoff_str
      isotopic_profile_calculator(MoleFormVecMat[i_mat, ], Elements_mass_abundance, peak_spacing, intensity_cutoff, UFA_IP_memeory_variables)
    }"
    IP_calculator <- gsub("intensity_cutoff_str", intensity_cutoff_str, IP_calculator)
    eval(parse(text = IP_calculator))
    ##
    ip_db_function <- function(i) {
      IPP <- IsotopicProfile_DataBase[[i]]
      x_100 <- which.max(IPP[, 2])
      L_IPP <- length(IPP[, 2])
      ##
      r13c_ip <- 0
      if (L_IPP > x_100) {
        M13C <- abs(IPP[, 1] - IPP[x_100, 1] - 1.00335484)
        M13C <- M13C[(x_100 + 1):L_IPP]
        x_101 <- which.min(M13C)[1]
        if (M13C[x_101] <= 0.015) {
          x_101 <- x_101 + x_100
          r13c_ip <- IPP[x_101, 2]/IPP[x_100, 2]*100
        }
      }
      c(IPP[x_100, 1], r13c_ip, x_100, L_IPP)
    }
    ##
    SpectraAnalysis_call <- function (i_pl) {
      peaklist <- loadRdata(paste0(input_path_pl, "/peaklist_", file_name_hrms[i_pl], ".Rdata"))
      ##
      outputer <- IDSL.IPA::IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i_pl])
      spectraList <- outputer[["spectraList"]]
      MS_polarity <- outputer[["MS_polarity"]]
      ##
      mzList.m <- do.call(rbind, lapply(1:L_MoleFormVecMat, function(j) {
        Annotation <- NULL
        ##
        x_pl <-mzRTindexer(peaklist[, 8], peaklist[, 3], mz_DataBase[j], RT_target_ion[j], mass_accuracy, delta_rt)
        if (!is.null(x_pl)) {
          ##
          R13C_PL <- peaklist[x_pl, 11]
          RangeScan <- peaklist[x_pl, 1]:peaklist[x_pl, 2]
          NumberScans <- length(RangeScan)
          ##
          IsotopicProfile <- IsotopicProfile_DataBase[[j]]
          size_IP <- SizeIP_IsotopicProfile_DataBase[j]
          R13C_IP <- R13C_DataBase[j]
          x_100 <- MAIso_IsotopicProfile_DataBase[j]
          MW_exp <- matrix(rep(0, size_IP*NumberScans), ncol = NumberScans)
          INT_exp <- MW_exp
          for (sc in 1:NumberScans) {
            PEAKS <- spectraList[[RangeScan[sc]]]
            for (Iso in 1:size_IP) {
              x_Iso <- which(abs(PEAKS[, 1] - IsotopicProfile[Iso, 1]) <= mass_accuracy)
              if (length(x_Iso) > 0) {
                if (length(x_Iso) > 1) {
                  x_Iso_min <- which.min(abs(PEAKS[x_Iso, 1] - IsotopicProfile[Iso, 1]))
                  x_Iso <- x_Iso[x_Iso_min[1]]
                }
                MW_exp[Iso, sc] <- PEAKS[x_Iso, 1]
                INT_exp[Iso, sc] <- PEAKS[x_Iso, 2]
              }
            }
          }
          sum_INT_exp <- rowSums(INT_exp)
          Ave_MW_exp <- rowSums(MW_exp*INT_exp)/sum_INT_exp
          PCS <- sum(sum_INT_exp*IsotopicProfile[, 2])/sqrt(sum(sum_INT_exp^2)*sum(IsotopicProfile[, 2]^2))*1000 # in per-mille
          NEME <- sqrt(sum((Ave_MW_exp - IsotopicProfile[, 1])^2)/size_IP)*1000 # in mDa
          MW_exp1 <- MW_exp
          MW_exp1[which(MW_exp1 > 0)] <- 1
          nd <- colSums(MW_exp1)
          Int_100 <- INT_exp[x_100, ]
          max_Int <- max(Int_100)
          x_80 <- which(Int_100/max_Int > 0.2)
          NDCS <- length(which(nd[x_80] == size_IP))
          L_80 <- x_80[length(x_80)] - x_80[1] + 1
          RCS <- NDCS/L_80*100
          Annotation <- c(x_pl, j, size_IP, molecular_formula_hill[j],
                          round(IsotopicProfile[x_100, 1], 5),
                          round(peaklist[x_pl, 8], 5),
                          mass_accuracy,
                          peaklist[x_pl, 3],
                          round(max_Int, 0),
                          round(NEME, 2),
                          round(PCS, 2),
                          round(R13C_PL, 2),
                          round(R13C_IP, 2),
                          NDCS,
                          round(RCS, 2))
          Annotation <- data.frame(Annotation)
          rownames(Annotation) <- c("PeakID", "ID_IonFormula", "sizeIP", "IonFormula", "m/z theoretical", "m/z peaklist", "Mass accuracy (Da)", "RetentionTime(min)", "PeakHeight", "NEME (mDa)", "PCS (per-mille)", "R13C peaklist (%)", "R13C theoretical (%)", "NDCS @ 80%", "RCS (%) @ 80%")
          if (exportSpectra == TRUE) {
            exp_spectra <- data.frame(matrix(cbind(Ave_MW_exp, sum_INT_exp/sum_INT_exp[x_100]*100), ncol = 2))
            colnames(exp_spectra) <- c("mz", "int")
            theo_spectra <- data.frame(IsotopicProfile)
            colnames(theo_spectra) <- c("mz", "int")
            lablel_spectra <- data.frame(cbind(round(exp_spectra[, 1], 5), sapply(1:size_IP, function(la_i) {max(c(exp_spectra[la_i, 2], theo_spectra[la_i, 2])) + 2})))
            max_int <- 1.15*max(c(theo_spectra[, 2], exp_spectra[, 2]))
            colnames(lablel_spectra) <- c("mz", "int")
            spectra_figure <- ggplot(data=exp_spectra, aes(x=mz, y=int)) +
              geom_segment(data=theo_spectra, aes(x=mz, xend=mz, y=0, yend=int, col = "Theoretical"), size = 6) +
              geom_segment(data=exp_spectra, aes(x=mz, xend=mz, y=0, yend=int, col = "Experimental"), size = 2) +
              scale_color_discrete(name = NULL) +
              xlab("m/z") + ylab("Intensity (%)") +
              scale_x_continuous(limits = c(theo_spectra[1, 1]-1, theo_spectra[size_IP, 1]+1), expand = c(0, 0)) +
              scale_y_continuous(limits = c(0, max_int), expand = c(0, 0)) +
              geom_text(data=lablel_spectra, aes(x=mz, y=int, label = round(mz, 5)), cex = 5) +
              annotate("text", x = (theo_spectra[1, 1]+theo_spectra[size_IP, 1])/2, y = max_int/1.05, label = paste0("[", molecular_formula_hill[j], "]", MS_polarity), size = 5) +
              labs(title = file_name_hrms[i_pl]) +
              theme_bw() + theme(legend.position = c(0.83, 0.95), legend.margin = margin(0, 0, 0, 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 text = element_text(size = 24), plot.title = element_text(size = 14))
            my_table <- tableGrob(Annotation, theme = ttheme_default()) ## create insert table
            spectra_figure <- grid.arrange(spectra_figure, my_table, ncol = 2) ## final result
            ggsave(filename=paste0("/UFA_spectra_", file_name_hrms[i_pl], "_", j, "_", molecular_formula_hill[j], "_",round(RT_target[j], 2), ".png"),
                   plot = spectra_figure,
                   device = "png",
                   path = output_path_spectra,
                   scale = 1,
                   width = 16,
                   height = 8,
                   units = "in",
                   dpi = 100)
          }
          Annotation <- t(rbind(file_name_hrms[i_pl], Annotation))
        }
        Annotation
      }))
      mzList.m
    }
    ##
    osType <- Sys.info()[['sysname']]
    if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      IsotopicProfile_DataBase <- foreach(counter = 1:L_MoleFormVecMat, .verbose = FALSE) %dopar% {
        IP_calculator(counter)
      }
      ip_db_mat <- foreach(counter = 1:L_MoleFormVecMat, .combine = "rbind", .verbose = FALSE) %dopar% {
        ip_db_function(counter)
      }
      ##
      mz_DataBase <- ip_db_mat[, 1]
      R13C_DataBase <- ip_db_mat[, 2]
      MAIso_IsotopicProfile_DataBase <- ip_db_mat[, 3]
      SizeIP_IsotopicProfile_DataBase <- ip_db_mat[, 4]
      ##
      AnnotatedSpectraTable <- foreach(counter = 1:L_PL, .combine = "rbind", .verbose = FALSE) %dopar% {
        SpectraAnalysis_call(counter)
      }
      stopCluster(clust)
      ##
    } else if (osType == "Linux") {
      ##
      IsotopicProfile_DataBase <- mclapply(1:L_MoleFormVecMat, function (counter) {
        IP_calculator(counter)
      }, mc.cores = number_processing_threads)
      ip_db_mat <- do.call(rbind, mclapply(1:L_MoleFormVecMat, function (counter) {
        ip_db_function(counter)
      }, mc.cores = number_processing_threads))
      ##
      mz_DataBase <- ip_db_mat[, 1]
      R13C_DataBase <- ip_db_mat[, 2]
      MAIso_IsotopicProfile_DataBase <- ip_db_mat[, 3]
      SizeIP_IsotopicProfile_DataBase <- ip_db_mat[, 4]
      ##
      AnnotatedSpectraTable <- do.call(rbind, mclapply(1:L_PL, function (counter) {
        SpectraAnalysis_call(counter)
      }, mc.cores = number_processing_threads))
      closeAllConnections()
    }
    ##
    print("Completed producing mass spectra!")
    ##
    if (exportedAnnotatedSpectraTable == TRUE) {
      AnnotatedSpectraTable <- data.frame(AnnotatedSpectraTable)
      rownames(AnnotatedSpectraTable) <- NULL
      colnames(AnnotatedSpectraTable) <- c("Filename", "PeakID", "ID_IonFormula", "sizeIP", "IonFormula", "m/z theoretical", "m/z peaklist", "Mass accuracy (Da)", "RetentionTime(min)", "PeakHeight", "NEME (mDa)", "PCS (per-mille)", "R13C peaklist (%)", "R13C theoretical (%)", "NDCS @ 80%", "RCS (%) @ 80%")
      ##
      return(AnnotatedSpectraTable)
    } else {
      return()
    }
  }
}
