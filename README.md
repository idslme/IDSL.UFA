# IDSL.UFA <img src='UFA_educational_files/Figures/IDSL.UFA-logo.PNG' width="250px" align="right" />

<!-- badges: start -->
[![Maintainer](https://img.shields.io/badge/maintainer-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.UFA)](https://cran.r-project.org/package=IDSL.UFA)
![](http://cranlogs.r-pkg.org/badges/IDSL.UFA?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.UFA?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.UFA)](https://cran.r-project.org/package=IDSL.UFA)


[![DOI](https://zenodo.org/badge/140601694.svg)](https://zenodo.org/record/7512923#.Y7m-DdXMJPY)
<!-- badges: end -->

[**United Formula Annotation (UFA)**](https://ufa.idsl.me/) by the [**Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME)**](https://www.idsl.me/) is a light-weight R package to annotate peaklists from the [**IDSL.IPA**](https://github.com/idslme/IDSL.IPA) package with molecular formula of a prioritized chemical space using an isotopic profile matching approach. The IDSL.UFA pipeline only requires MS1 for formula annotation.

	install.packages("IDSL.UFA")

## <img src='UFA_educational_files/Figures/IDSL.UFA-TOC_Art.png' align="right" />

## Features of IDSL.UFA

1) Parameter selection through a well-described [parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.UFA/main/UFA_parameters.xlsx)
2) Generating comprehensive *in-silico* theoretical libraries using natural isotopic distribution profiles
3) Annotating high-throughput large scale chromatographic peaklists
4) Aggregating annotated molecular formulas on the aligned peak table
5) Generating batch untargeted isotopic profile match figures
6) Compatibility with parallel processing in Windows and Linux environments

## Workflow
To annotate your mass spectrometry data (**mzXML**, **mzML**, **netCDF**), mass spectrometry data should be processed using the [IDSL.IPA](https://github.com/idslme/IDSL.IPA) workflow to acquire chromatographic information of the peaks (***m/z-RT***). When the chromatographic information of individual and aggregated aligned peaklists were generated using the [IDSL.IPA](https://github.com/idslme/IDSL.IPA) workflow, download the [UFA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.UFA/main/UFA_parameters.xlsx) and select the parameters accordingly and then use this spreadsheet as the input for the IDSL.UFA workflow:

	library(IDSL.UFA)
	UFA_workflow("Address of the UFA parameter spreadsheet")

Visit [**wiki**](https://github.com/idslme/IDSL.UFA/wiki) for detailed documentations and tutorials for the [**list of consistent labeled isotopes**](https://github.com/idslme/IDSL.UFA/wiki/Consistent-Labeled-Isotopes), [**Standard Adduct Type**](https://github.com/idslme/IDSL.UFA/wiki/Standard-Adduct-Type), [**Definitions of Peak Spacing and Intensity Cutoff**](https://github.com/idslme/IDSL.UFA/wiki/Peak-Spacing-and-Intensity-Cutoff), [**Isotopic Profile DataBase (IPDB)**](https://github.com/idslme/IDSL.UFA/wiki/Isotopic-Profile-DataBase-(IPDB)), [**PubChem molecular formula database for IDSL.UFA**](https://github.com/idslme/IDSL.UFA/wiki/PubChem-molecular-formula-database-for-IDSL.UFA), [**Definitions of NDCS and RCS**](https://github.com/idslme/IDSL.UFA/wiki/NDCS-RCS), and [**Molecular formula class detection**](https://github.com/idslme/IDSL.UFA/wiki/Molecular-formula-class-detection).

## Note
The IDSL.UFA pipeline originally was developed to annotate IDSL.IPA peaklists that were generated using <sup>12</sup>C/<sup>13</sup>C isotopologue pairs by the IDSL.IPA pipeline. Nevertheless, IDSL.UFA can still annotate IDSL.IPA peaklists with non-carbon ion pairs when the fifth coefficient in the score function in the equation 6 in the main manuscript is zero to neutralize interferences of R<sup>13</sup>C values. The score coefficients can be adjusted through **PARAM0023** in the [UFA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.UFA/main/UFA_parameters.xlsx).

## Citation
Fakouri Baygi, S., Banerjee S. K., Chakraborty P., Kumar, Y. Barupal, D.K. [IDSL.UFA assigns high confidence molecular formula annotations for untargeted LC/HRMS datasets in metabolomics and exposomics](https://pubs.acs.org/doi/10.1021/acs.analchem.2c00563). *Analytical Chemistry*, **2022**, *94(39)*, 13315�13322.


Fakouri Baygi, S., Kumar, Y. Barupal, D.K. [IDSL. IPA characterizes the organic chemical space in untargeted LC/HRMS datasets](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00120). *Journal of proteome research*, **2022**, *21(6)*, 1485-1494.
