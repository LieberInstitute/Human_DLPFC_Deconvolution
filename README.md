# Human_DLPFC_Deconvolution

## Internal

JHPCE location: `/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution`

## File locations

This is an integrative, multi-assay project including individual-matched data generated from human dorsolateral prefrontal cortex (DLPFC) using the RNAscope, 10X Chromium, and Illumina HiSeq platforms. 

The main directory path on JHPCE is located at:

`/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/`

Note that integrative analyses and content for the deconvolution methods paper are stored in the subdirectory `Human_DLPFC_Deconvolution`.

### 10X Chromium single-nucleus RNA-seq datasets

The snRNA-seq datasets for this project currently include 9 donors and 10 samples. For the deconvolution method paper, we access these data as a `SingleCellExperiment` object. The latest (final?) version is a file called `sce_DLPFC.Rdata` located at the subdirectory:

`DLPFC_snRNAseq/processed-data/sce/`

### Bulk RNA-seq datasets

The bulk RNA-seq datasets include 6 samples groups comprised of mRNA isolated from either nucleus, cytoplasm, or bulk, and prepped with either polyA selection or rRNA depletion strategies. Data come from 19 samples, meaning 113 total data points. 

These bulk RNA-seq data are currently being QC'd, but they will be accessible as a `SummarizedExperiment` object located at the subdirectory:

`Human_DLPFC_Deconvolution/processed-data/01_SPEAQeasy/`

## Variables and data dictionaries

The deconvolution project makes use of a number of metadata attributes and variables in the results files mentioned about. This section describes the key terms and definitions of these attributes and variables for the deconvolution method paper.

### Cell type labels

Cell type labels for `snRNA-seq` datasets are determined from the variable `cellType_broad_hc`. This can be accessed from the `sce` object in various ways such as `colData(sce)$cellType_broad_hc` or `sce[["cellType_broad_hc"]]`. Note that the deconvolution methods paper focuses on just 6 cell types of interest, and these are identified from among the cell type labels in the `cellType_broad_hc` variable.

Cell type labels for RNAscope experiments are obtained from the image analysis outputs produced by the HALO software. In brief, outputs each contain a series of columns corresponding to the cell type makers. Since each row in these outputs corresponds to an individual detected nucleus, we simply look at which marker is positive for that nucleus to determine its cell type. Cell type proportions and abundances are then calculated from these outputs.

Cell type labels aren't available for the bulk RNA-seq and other datasets produced for this project.

### Data dictionaries for marker tables, results, etc.

## Terminology

Throughout the deconvolution project we try to use standard terms to refer to key project entities. Here are some of the
key terms to be aware of when using these project files and understanding analysis outputs.

### The deconvolution equation

We take deconvolution to be the prediction of cell type amounts in a mixed sample by leveraging data from a non-mixed sample. For the deconvolution method paper, we focus on predicting either cell abundances or proportions for each of 6 cell types in bulk RNA-seq data by leveraging a single-nucleus RNA-seq reference dataset. The deconvolution equation looks like:

$Y = \pi * Z$

Where $Y$ is a matrix of mixed signals having dimensions $[G,J]$, $\pi$ is a vector of cell type amounts (proportions or abundances), and $Z$ is a matrix of non-mixed signals having dimensions $[G,K]$. The dimensions of these objects include the $G$ marker gene set, $J$ samples with mixed signal, and $K$ cell types.

### Strict deconvolution

Strict deconvolution refers to solving for $\pi$ given a set of matrices $[Y,Z]$. In other words, if we assume we have a set of matrices $Y$ and $Z$, we then solve for $\pi$ based on these.

### Deconvolution preparation

This simply refers to steps taken to prepare $Y$ and $Z$ prior to performing deconvolution. This may include differential weighting of the marker genes, bias corrections, data transformations and rescaling, etc.

### Marker genes

The $G$ marker genes represent markers of cell types which we can use in deconvolution. Note, this is the "row" dimension in the $Y$ and $Z$ matrices.

### Samples

The $J$ samples are the samples having mixed signal, for which we predict cell type amounts. Note, this is the "column" dimension in the matrix $Y$.

### Cell types

The $K$ cell types represent the cell types for which we deconvolute signal and obtain amount estimates. Note, this is the "column" dimension in matrix $Y$ and the length of the vector represented by $\pi$.

## Funding sources

* 1R01MH123183-01 ([grantome_link](https://grantome.com/grant/NIH/R01-MH123183-01))