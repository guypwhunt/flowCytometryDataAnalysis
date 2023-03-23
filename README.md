# Specialised Pro-Resolving Mediators and Resolvin Receptors Analysis

This GitHub Repository contains the R scripts to run the differential analysis on specialised pro-resolvin mediators as well as the clustering, differential expression analysis and survival analysis.

## Contacts
If you have an issues running the code please feel free to raise an issue on the GitHub or contact either:

- Guy Hunt <guy.hunt@kcl.ac.uk>
- Alfredo Iacoangeli <alfredo.iacoangeli@kcl.ac.uk>

## Folder Descriptions
The "R" folder contains all the R scripts to rerun the analysis. Within the R folder there are several R scripts and folders which are described below:

- 01_specialised_pro_resolving_mediators_analysis contains the scripts to run differential expression analysis on the specialised pro-resolving mediators
- 02_GPR18_isotypes_analysis contains the scripts to run the isotype analysis on GPR18 and Chem23
- 03_GPR32_isotypes_analysis contains the scripts to run the isotype analysis on GPR32 and FPRL1
- Scripts 04 - 11 contain the R scripts to pre-process, perform clustering and differential expression analysis for B Cells, Monocytes, T Cells and Senescent T Cells for both GPR32 and GPR18
- 12_Visulise_the_Cell_Populations contains the scripts to visulise the UMAP analysis
- 13_Visulise_the_Cell_Populations_Stability contains the scripts to visulise the cell population stability analysis
- 14_Visulise_the_Differential_Expression_Results contains the scripts to visulise the differential expression analysis
- 15_Visulise_the_Clinical_Data contains the scripts to visulise the umap analysis of the clinical data
- 16_Survival_Analysis contains the scripts to perform and visulise the survival analysis on GPR32 and GPR18
- 00_datasets.R contains the small variables that are used in multiple analyses
- 01_functions.R contains the functions that are used in multiple analyses
