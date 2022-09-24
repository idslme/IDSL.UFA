# IDSL.UFA<img src='UFA_educational_files/Figures/IDSL.UFA-logo.PNG' width="250px" align="right" />

<!-- badges: start -->
[![Maintainer](https://img.shields.io/badge/maintainer-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.UFA)](https://cran.r-project.org/package=IDSL.UFA)
![](http://cranlogs.r-pkg.org/badges/IDSL.UFA?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.UFA?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.UFA)](https://cran.r-project.org/package=IDSL.UFA)


[![DOI](https://zenodo.org/badge/140601694.svg)](https://zenodo.org/record/7065107#.YxvObaHMK70)
<!-- badges: end -->

[**United Formula Annotation (UFA)**](https://ufa.idsl.me/) by the [**Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME)**](https://www.idsl.me/) is an R pipeline to annotate peaklists from the [**IDSL.IPA**](https://github.com/idslme/IDSL.IPA) package with molecular formula of a prioritized chemical space using an isotopic profile matching approach. The IDSL.UFA pipeline only requires MS1 for formula annotation.

	install.packages("IDSL.UFA")

## <img src='UFA_educational_files/Figures/IDSL.UFA-TOC_Art.png' align="right" />

## Workflow
To annotate your mass spectrometry data (**mzXML**, **mzML**, **netCDF**), mass spectrometry data should be processed using the [IDSL.IPA](https://github.com/idslme/IDSL.IPA) workflow to acquire chromatographic information of the peaks (***m/z-RT***). When the chromatographic information of individual and aggregated aligned peaklists were generated using the [IDSL.IPA](https://github.com/idslme/IDSL.IPA) workflow, download the [UFA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.UFA/main/UFA_parameters.xlsx) and select the parameters accordingly and then use this spreadsheet as the input for the IDSL.UFA workflow:

	library(IDSL.UFA)
	UFA_workflow("Address of the UFA parameter spreadsheet")

Visit [**wiki**](https://github.com/idslme/IDSL.UFA/wiki) for the detailed documentations and tutorials on specific functions for [**Isotopic Profile DataBase (IPDB)**](https://github.com/idslme/IDSL.UFA/wiki/Isotopic-Profile-DataBase-(IPDB)), [**PubChem molecular formula database for IDSL.UFA**](https://github.com/idslme/IDSL.UFA/wiki/PubChem-molecular-formula-database-for-IDSL.UFA), and [**Molecular formula class detection**](https://github.com/idslme/IDSL.UFA/wiki/Molecular-formula-class-detection).

Visit https://ufa.idsl.me/ for the detailed documentation and tutorial.

## Citation
Fakouri Baygi, S., Banerjee S. K., Chakraborty P., Kumar, Y. Barupal, D.K. [IDSL.UFA assigns high confidence molecular formula annotations for untargeted LC/HRMS datasets in metabolomics and exposomics](https://pubs.acs.org/doi/10.1021/acs.analchem.2c00563). *Analytical Chemistry*, **2022**.


Fakouri Baygi, S., Kumar, Y. Barupal, D.K. [IDSL. IPA characterizes the organic chemical space in untargeted LC/HRMS datasets](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00120). *Journal of proteome research*, **2022**, *21(6)*, 1485-1494.
