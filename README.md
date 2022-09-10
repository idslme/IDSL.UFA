# IDSL.UFA<img src='UFA_educational_files/Figures/IDSL.UFA-logo.PNG' width="250px" align="right" />

<!-- badges: start -->
[![Maintainer](https://img.shields.io/badge/maintainer-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.UFA)](https://cran.r-project.org/package=IDSL.UFA)
![](http://cranlogs.r-pkg.org/badges/IDSL.UFA?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.UFA?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.UFA)](https://cran.r-project.org/package=IDSL.UFA)


[![DOI](https://zenodo.org/badge/140601694.svg)](https://zenodo.org/record/7065107#.YxvObaHMK70)
<!-- badges: end -->

[**United Formula Annotation (UFA)**](https://ufa.idsl.me/) by the [**Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME)**](https://www.idsl.me/) is an R pipeline to annotate peaklists from the [**IDSL.IPA**](https://github.com/idslme/IDSL.IPA) package with molecular formula using an isotopic profile matching approach. The IDSL.UFA pipeline is especially beneficial when MS/MS data is not available.

	install.packages("IDSL.UFA")

## <img src='UFA_educational_files/Figures/IDSL.UFA-TOC_Art.png' align="right" />

## Workflow
To annotate your mass spectrometry data (**mzXML**, **mzML**, **netCDF**), mass spectrometry data should be processed using the [IDSL.IPA](https://github.com/idslme/IDSL.IPA) workflow to acquire chromatographic information of the peaks (***m/z-RT***). When the chromatographic information of individual and aggregated aligned peaklists were generated using the [IDSL.IPA](https://github.com/idslme/IDSL.IPA) workflow, download the [UFA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.UFA/main/UFA_parameters.xlsx) and select the parameters accordingly and then use this spreadsheet as the input for the IDSL.UFA workflow:

	library(IDSL.UFA)
	UFA_workflow("Address of the UFA parameter spreadsheet")

Visit https://ufa.idsl.me/ for the detailed documentation and tutorial.

## Isotopic Profile DataBase (IPDB)
IDSL.UFA requires a library of calculated isotopic profiles (**IPDB**) to efficiently annotate chromatographic peaks with molecular formulas. **IPDB** libraries can be saved and re-used in similar workflows which is also necessary for consistency in population-size studies.
Generally, **IPDB** objects are R lists consisting of eight primary objects including:

**logIPDB**: Parameters used to create the **IPDB** object

**AggregatedList**: A list of rounded mass and IDs

**MassMAIso**: A vector of mass of the most abundant isotopologues

**MolecularFormulaDB**: A list of 1) a vector elements used in the molecular ions and 2) a matrix of stacked vectorized molecular ions

**IsotopicProfile**: A list of theoretical isotopic profiles

**R<sup>13</sup>C**: A vector of theoretical R<sup>13</sup>C values

**IndexMAIso**: A vector of indices of the most abundant isotopologues in the isotopic profiles

**IPsize**: A vector of number of isotopologues in the isotopic profiles

#
We updated the **IPDB** structure for faster performance after v1.5 and older IPDBs do not work with IDSL.UFA version >= 1.5. We generated IPDBs consistent with IDSL.UFA >= 1.5 for [**blood exposome**](https://bloodexposome.org/#/dashboard), [**EPA CompTox chemicals dashboard**](https://comptox.epa.gov/dashboard/), [**FDA substance registry**](https://www.fda.gov/industry/fda-data-standards-advisory-board/structured-product-labeling-resources), [**IDSL.Exposome**](https://chemcor.idsl.site/exposomechemicals/), [**LIPID MAPS**](https://www.lipidmaps.org/), [**RefMet**](https://metabolomicsworkbench.org/databases/refmet/index.php) and [**PubChem**](https://pubchem.ncbi.nlm.nih.gov/) databases in positive and negative modes. These **IPDB** libraries can be accessed using this [**link**](https://zenodo.org/record/7065107/files/IPDB_v1.5.zip?download=1) and this link for [**PubChem IPDB**](https://zenodo.org/record/7065107/files/PubChem07092021.zip?download=1). These **IPDB** libraries were generated presuming occurance of c("**[M+H]<sup>+</sup>**", "**[M+Na]<sup>+</sup>**", "**[M-H<sub>2</sub>O+H]<sup>+</sup>**") and c("**[M-H]<sup>-</sup>**", "**[M-H<sub>2</sub>O-H]<sup>-</sup>**") ionization pathways in positive and negative modes, respectively. Therefore, numbers of molecular formula ions in IPDBs are roughly a factor of the number of ionization pathways multiplied by the number of intact molecular formulas. Non-carbon-containing compounds are excluded from IPDBs since IDSL.IPA cannot detect non-carbon-containing compounds in the first place.


**Tip:** In case, molecular formula ions used in an **IPDB** are needed the following command can be used:

	 Elements <- IPDB[["MolecularFormulaDB"]][["Elements"]]
	 MolVecMat <- IPDB[["MolecularFormulaDB"]][["MolecularFormulaMatrix"]]
	 Molecular_Formulas_IPDB <- IDSL.UFA::hill_molecular_formula_printer(Elements, MolVecMat, number_processing_threads = 1)

## [**PubChem**](https://pubchem.ncbi.nlm.nih.gov/) molecular formula database for IDSL.UFA
In addition to [**positive and negative IPDBs from the PubChem molecular formula database**](https://zenodo.org/record/7065107/files/PubChem07092021.zip?download=1), a database of frequency of molecular formulas in the [**PubChem**](https://pubchem.ncbi.nlm.nih.gov/) molecular formula database are also provided in this [**link**](https://zenodo.org/record/7065107/files/PubChem_MolecularFormula_Freq_Database.Rdata?download=1). When the ionization pathway of the analytical platform is known, intact molecular formulas can be obtained and then check its presence and its frequency in the [**PubChem**](https://pubchem.ncbi.nlm.nih.gov/) molecular formula database. Molecular formula database search can be activated through PARAM0022-PARAM0024 in the [UFA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.UFA/main/UFA_parameters.xlsx) for individual file annotation.

## Molecular formula class detection
In many instances, the molecular formula enumeration method usually results with many molecular formulas. We develoepd a module to sort these molecular formulas based on their classes to facilitate similar molecular formulas. The aligned annotated molecular formula tables are the best source to find these related molecular formulas in a study.

**From manuscript:** Many compounds belong to a chemical class with a distinct sub-structure pattern such as polychlorinated biphenyl (PCBs), polybrominated diphenyl ethers (PBDEs), polycyclic aromatic hydrocarbons (PAHs), perfluoroalkyl substances (PFAS), lipids and phthalates, etc. The formula annotations generated via the enumerated chemical space (ECS) approach were processed to detect such classes within a list of formulas. The IDSL.UFA function *detect_formula_sets* was used to detect 1) constant &Delta;H/&Delta;C ratios for polymeric (&Delta;H/&Delta;C = 2) and cyclic (&Delta;H/&Delta;C = 1/2) chain progressions within polymeric and cyclic classes (Table S.2- S.4) and 2) a constant number of carbons and fixed summation of hydrogens and halogens (&Sigma;(H+Br+Cl+F+I)) representing classes similar to PCBs, PBDEs (Table S.5).

	detect_formula_sets(molecular_formulas, ratio_delta_HBrClFI_C,
	mixed.HBrClFI.allowed, min_molecular_formula_class, max_number_formula_class,
	number_processing_threads = 1)

**molecular_formulas:** a vector of molecular formulas

**ratio_delta_HBrClFI_C:** c(**2**, **1/2**, **0**). **2** to detect structures with linear carbon chains such as PFAS, lipids, chlorinated paraffins, etc. **1/2** to detect structures with cyclic chains such as PAHs. **0** to detect molecular formulas with a fixed structures but changing H/Br/Cl/F/I atoms similar to PCBs, PBDEs, etc.

**mixed.HBrClFI.allowed:** c(`TRUE`, `FALSE`). Select `FALSE` to detect halogenated-saturated compounds similar to PFOS or select `TRUE` to detect mixed halogenated compounds with hydrogen.

**min_molecular_formula_class:** minimum number of molecular formulas in each class. This number should be greater than or equal to **2**.

**max_number_formula_class:** maximum number of molecular formulas in each class

**number_processing_threads:** Number of processing threads for multi-threaded computations


	library(IDSL.UFA)
	##
	molecular_formulas <- c("C3F7O3S", "C4F9O3S", "C5F11O3S", "C6F9O3S", "C8F17O3S",
	"C9F19O3S", "C10F21O3S", "C7ClF14O4", "C10ClF20O4", "C11ClF22O4", "C11Cl2F21O4",
	"C12ClF24O4")
	##
	ratio_delta_HBrClFI_C <- 2 # to aggregate polymeric classes
	mixed.HBrClFI.allowed <- FALSE # To detect only halogen saturated classes
	min_molecular_formula_class <- 2
	max_number_formula_class <- 20
	##
	classes <- detect_formula_sets(molecular_formulas, ratio_delta_HBrClFI_C,
	mixed.HBrClFI.allowed, min_molecular_formula_class, max_number_formula_class,
	number_processing_threads = 1)


## Citation
Fakouri Baygi, S., Banerjee S. K., Chakraborty P., Kumar, Y. Barupal, D.K. [IDSL.UFA assigns high confidence molecular formula annotations for untargeted LC/HRMS datasets in metabolomics and exposomics](https://www.biorxiv.org/content/10.1101/2022.02.02.478834v2). *Analytical Chemistry*, **2022**.


Fakouri Baygi, S., Kumar, Y. Barupal, D.K. [IDSL. IPA characterizes the organic chemical space in untargeted LC/HRMS datasets](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00120). *Journal of proteome research*, **2022**, *21(6)*, 1485-1494.
