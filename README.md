# IDSL.UFA <img src='UFA_educational_files/Figures/IDSL.UFA-logo.PNG' width="250px" align="right" />

<!-- badges: start -->
[![Developed-by](https://img.shields.io/badge/Developed_by-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.UFA)](https://cran.r-project.org/package=IDSL.UFA)
![](http://cranlogs.r-pkg.org/badges/IDSL.UFA?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.UFA?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.UFA)](https://cran.r-project.org/package=IDSL.UFA)


[![DOI](https://zenodo.org/badge/140601694.svg)](https://zenodo.org/record/7512923#.Y7m-DdXMJPY)
<!-- badges: end -->

**United Formula Annotation (UFA)** by the [**Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME)**](https://www.idsl.me/) is a light-weight R package to annotate peaklists from the [**IDSL.IPA**](https://github.com/idslme/IDSL.IPA) package with molecular formula of a prioritized chemical space using an isotopic profile matching approach. The IDSL.UFA pipeline only requires MS1 for molecular formula annotation.

## <img src='UFA_educational_files/Figures/IDSL.UFA-TOC_Art.png' align="right" />

## Table of Contents

- [Background](https://github.com/idslme/IDSL.UFA#background)
- [Features of IDSL.UFA](https://github.com/idslme/IDSL.UFA#features-of-idslufa)
- [Installation](https://github.com/idslme/IDSL.UFA#installation)
- [Workflow](https://github.com/idslme/IDSL.UFA#workflow)
- [Quick Batch Example](https://github.com/idslme/IDSL.UFA#quick-batch-example)
- [Wiki](https://github.com/idslme/IDSL.UFA#wiki)
- [Citation](https://github.com/idslme/IDSL.UFA#citation)

## Background

***molecular formulas*** are fundamental property of chemical compounds and represent their elemental compositions. Assigning molecular formulas to peaks in data generated using untargeted LC/HRMS can help in gaining biological insights from metabolomics and exposomics datasets. Molecular formula annotation can also complement the peak annotation pipelines that need MS2 spectra to assign a structural identity to a peak. Formulas can be assigned using only MS1 spectral data which is available for every sample analyzed using a LC/HRMS instrument in metabolomics and exposomics studies.  

Because of the naturally occuring isotope atoms for each element, MS1 spectral data have more than one mass to charge ratio (m/z) values observed for an ionized species. The isotopic pattern for a chemical structure can be accurately predicted using a set of combinatorial rules that uses atomic mass tables provided by the [International Union of Pure and Applied Chemistry (IUPAC)](https://www.isotopesmatter.com). To assign a molecular formula, the theoretical isotopic profile of carbon-containing compounds can be queried against the MS1 spectral data using a set of matching criteria and ranking system. The universality of molecular formula assignment can allow almost all commercial and academic software to process untargeted LC/HRMS datasets to search a single or list of molecular formulas against the raw MS1 data. Community guidelines for peak annotation also recommend performing the molecular formula assignment step on untargeted LC/HRMS datasets. 

While existing solutions offer a straightforward solution to match theoretical isotopic patterns against the MS1 spectral data, there is still unmet needs to improve the workflow for larger studies and various sources of molecular formula. This is important for exposomics studies where we do expect to see many more compounds from formula sources other than common metabolite databases.

## Features of IDSL.UFA

1) Parameter selection through a user-friendly and well-described [parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.UFA/main/UFA_parameters.xlsx)
2) Analyzing population size untargeted studies (n > 500)
3) Generating comprehensive *in-silico* theoretical libraries (known as [IPDB](https://github.com/idslme/IDSL.UFA/wiki/Isotopic-Profile-DataBase-(IPDB))) using natural isotopic distribution profiles
4) Aggregating annotated molecular formulas on the aligned peak table. This is a very unique feature that is only presented by IDSL.UFA. To familiarize with this statistical mass spectrometry feature, try **PARAM0006** in the `parameters` tab in the [UFA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.UFA/main/UFA_parameters.xlsx)
5) Ranking candidate molecular formulas for a peak
6) Generating batch untargeted isotopic profile match figures
7) Parallel processing in Windows and Linux environments

## Installation

	install.packages("IDSL.UFA")

## Workflow

To annotate your mass spectrometry data (**mzXML**, **mzML**, **netCDF**), mass spectrometry data should be processed using the [IDSL.IPA](https://github.com/idslme/IDSL.IPA) workflow to acquire chromatographic information of the peaks (***m/z-RT***). When the chromatographic information of individual and aggregated aligned peaklists were generated using the [IDSL.IPA](https://github.com/idslme/IDSL.IPA) workflow, download the [UFA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.UFA/main/UFA_parameters.xlsx) and select the parameters accordingly and then use this spreadsheet as the input for the IDSL.UFA workflow:

	library(IDSL.UFA)
	UFA_workflow("Address of the UFA parameter spreadsheet")

## Quick Batch Example

Follow these steps for a quick case study (n = 33) [ST002263](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST002263&DataMode=AllData&ResultType=1) which has Thermo Q Exactive HF hybrid Orbitrap data collected in the HILIC-ESI-POS/NEG modes. 

1. Process raw mass spectrometry data to generate chromatographics information using the method described by [IDSL.IPA](https://github.com/idslme/IDSL.IPA#quick-batch-example) for this study.

2. Download these pre-calculated [IPDBs](https://zenodo.org/record/7512923/preview/IPDB_v1.8.zip#tree_item16) and use positive or negative mode IPDB from RefMetDB folder according to the IDSL.IPA folder results. RefMet represents a [**Ref**erence list of **Met**abolite names](https://www.metabolomicsworkbench.org/databases/refmet/).

3. IDSL.UFA requires 30 parameters distributed into 4 separate sections for a full scale analysis. For this study, use default parameter values presented in the [UFA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.UFA/main/UFA_parameters.xlsx). Next, provide information for 
	
	3.1. **PARAM0004** for the *Address of the IPDB (.Rdata)*
	
	3.2. **PARAM0005** and **PARAM0006** should be **YES**
	
	3.3. **PARAM0009** for *HRMS data location address (MS1 level HRMS data)*
	
	3.4. **PARAM0011** for *Address of the `peaklists` directory generated by the IDSL.IPA workflow*
	
	3.5. **PARAM0012** for *Address of the `peak_alignment` directory generated by the IDSL.IPA workflow*
	
	3.6. **PARAM0014** for *Output location address (MS1 processed data)*
	
	3.7. You may also increase the number of processing threads using **PARAM0008** according to your computational power

4. Run this command in the R/Rstudio console or terminal

```
library(IDSL.UFA)
UFA_workflow("Address of the UFA parameter spreadsheet")
```

5. You see the results in the address you provided for **PARAM0014** including:

	5.1. Individual annotated peaklists with molecular formulas for each HRMS file in the *annotated_mf_tables* directory in the *.Rdata* and *.csv* formats
	
	5.2. Aligned molecular formula table in the *aligned_molecular_formula_table* directory in the *.Rdata* and *.csv* formats. We strongly recommend to familiarize yourself with the structure of this table to find the most probable candidate molecular formulas.
	
	5.3. If you had selected numbers greater than **0** for **PARAM0024**, match spectra figures are presented in the *UFA_spectra* folder.

## [**Wiki**](https://github.com/idslme/IDSL.UFA/wiki)

1. [**A population size study with 499 individual mass spectrometry file**](https://github.com/idslme/IDSL.UFA/wiki/IDSL.UFA-for-MTBLS1684-study)
2. [**List of consistent labeled isotopes**](https://github.com/idslme/IDSL.UFA/wiki/Consistent-Labeled-Isotopes)
3. [**Standard Adduct Type**](https://github.com/idslme/IDSL.UFA/wiki/Standard-Adduct-Type)
4. [**Definitions of Peak Spacing and Intensity Cutoff**](https://github.com/idslme/IDSL.UFA/wiki/Peak-Spacing-and-Intensity-Cutoff)
5. [**Isotopic Profile DataBase (IPDB)**](https://github.com/idslme/IDSL.UFA/wiki/Isotopic-Profile-DataBase-(IPDB))
6. [**PubChem molecular formula database for IDSL.UFA**](https://github.com/idslme/IDSL.UFA/wiki/PubChem-molecular-formula-database-for-IDSL.UFA)
7. [**Definitions of NEME and PCS**](https://github.com/idslme/IDSL.UFA/wiki/NEME-PCS)
8. [**Definitions of NDCS and RCS**](https://github.com/idslme/IDSL.UFA/wiki/NDCS-RCS)
9. [**Score Coefficients Optimization**](https://github.com/idslme/IDSL.UFA/wiki/Score-Coefficients-Optimization)
10. [**Molecular formula class detection**](https://github.com/idslme/IDSL.UFA/wiki/Molecular-formula-class-detection)

## Citation

[1] Fakouri Baygi, S., Banerjee S. K., Chakraborty P., Kumar, Y. Barupal, D.K. [IDSL.UFA assigns high confidence molecular formula annotations for untargeted LC/HRMS datasets in metabolomics and exposomics](https://pubs.acs.org/doi/10.1021/acs.analchem.2c00563). *Analytical Chemistry*, **2022**, *94(39)*, 13315-13322.

[2] Fakouri Baygi, S., Kumar, Y. Barupal, D.K. [IDSL. IPA characterizes the organic chemical space in untargeted LC/HRMS datasets](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00120). *Journal of proteome research*, **2022**, *21(6)*, 1485-1494.
