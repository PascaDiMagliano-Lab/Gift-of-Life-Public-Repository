This repository contains code used for single-cell and spatial transcriptomics data analysis in the [Analysis of donor pancreases defines the transcriptomic signature and microenvironment of early pre-neoplastic pancreatic lesions](link) paper.

## How to use this repo to reproduce the analysis done in the paper?

* `data` folder contains instructions for downloading data from GEO and samples metadata spreadsheet, along with other low-size input data.
* `scripts` folder contains sub-folders for each section of the analysis.
* `results` folder will be the output directory for the scripts. 
* `doc` contains reports for the results of the analysis. 
* `src` folder contains utility functions that are called within other scripts for a cleaner code.

*** 

## Data Acquisition

Instructions for data acquisition is described under [`Data_Acquisition`](scripts/Data_Acquisition) folder. Input data will be stored under `data` folder.

## Data Processing

Fastq files were aligned to hg38 using CellRanger v6. Ambient RNA correction was done using [SoupX](https://github.com/constantAmateur/SoupX) as described in [`01_ambient_RNA_correction.R`](scripts/Data_Processing/01_ambient_RNA_correction.R). Samples underwent QC to remove low quality cells using [`02_data_procesisng_QC.R`](scripts/Data_Processing/02_data_procesisng_QC.R). Samples were then integrated using the recommeded [Seurat rPCA MNN integration](https://satijalab.org/seurat/articles/integration_rpca.html) using the [`03_data_integration_script.R`](scripts/Data_Processing/03_data_integration_script.R). Clusters were then labeled using previously publised markers using [`04_clusters_annotation.R`](scripts/Data_Processing/04_clusters_annotation.R) 

## Exploratory Analysis

Differential abundance analysis was done using [miloR R Package](https://github.com/MarioniLab/miloR) as described in [`01_differential_abundance_with_miloR.R`](scripts/Exploratory_Analysis/01_differential_abundance_with_miloR.R). Pseudobulk analysis was done using [`02_Pseudobulk.R`](scripts/Exploratory_Analysis/02_Pseudobulk.R). [CellChat](https://github.com/sqjin/CellChat) was used to perform Cell-cell communication analysis as described in [`03_CellChat_analysis.R`](scripts/Exploratory_Analysis/03_CellChat_analysis.R). Interactome circos plots were generated using [`04_interactome_analysis.R`](scripts/Exploratory_Analysis/04_interactome_analysis.R) script. KRAS and previously published gene set scoring was done utilizing [AUCell](https://github.com/aertslab/AUCell) as described in [`05_geneset_scoring.R`](scripts/Exploratory_Analysis/05_geneset_scoring.R) script. Differential gene expression analysis and UMAP plotting was done using Seurat v4 functions; [`FindMarkers()`](https://satijalab.org/seurat/reference/findmarkers), [`DimPlot()`](https://satijalab.org/seurat/reference/dimplot) [`FeaturePlot()`](https://satijalab.org/seurat/reference/featureplot) and [`VlnPlot()`](https://satijalab.org/seurat/reference/vlnplot)

## Spatial Transcriptomics

Spatial transcriptomics was done using the Nanostring GeoMx platform. Segments QC was done to remove the low quality segments and genes with low detection level as recommeded by Nanostring using [GeoMxWorkflows R Package](https://github.com/Nanostring-Biostats/GeoMxWorkflows), then, differential gene expression was done in [`01_DGE_analysis.R`](scripts/Spatial_Transcriptomics_Analysis/01_DGE_analysis.R). Mapping of marker gene sets derived from spatial transcriptomics onto the single cell data was done using [`02_mapping_markers_to_SC.R`](scripts/Spatial_Transcriptomics_Analysis/02_mapping_markers_to_SC.R)
 
