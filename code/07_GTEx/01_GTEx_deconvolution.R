
library("recount3")
# library(jaffelab)
# library(SummarizedExperiment)
library("SingleCellExperiment")
library("xbioc")
library("MuSiC")
library("BisqueRNA")
library("tidyverse")
# library(reshape2)
# library(compositions)
library("here")
library("sessioninfo")
# library(DeconvoBuddies)

#### Load GTEx data with recount3 ####

## Find all available human projects
human_projects <- available_projects()
## Find the project you are interested in
proj_info <- subset(
  human_projects,
  project == "BRAIN" & file_source == "gtex")

proj_info
# project organism file_source      project_home project_type n_samples
# BRAIN    human        gtex data_sources/gtex data_sources      2931

## Create RSE (defaults to Gencode v26, could also be Gencode v29)
rse_gene_brain_gtex <- create_rse(proj_info)
## Compute read counts (not needed if you want to go to RPKMs directly)
# assay(rse_gene_brain_gtex, "counts") <- compute_read_counts(rse_gene_brain_gtex)

## Drop bad GTEx samples
table(rse_gene_brain_gtex$gtex.smafrze)
# EXCLUDE  RNASEQ 
# 261    2670
rse_gene_brain_gtex <- rse_gene_brain_gtex[, rse_gene_brain_gtex$gtex.smafrze != "EXCLUDE"]

## Compute RPKMs
assay(rse_gene_brain_gtex, "counts") <- transform_counts(rse_gene_brain_gtex)
# assay(rse_gene_brain_gtex, "RPKM") <- recount::getRPKM(rse_gene_brain_gtex, length_var = "score")

## Brain region info
table(rse_gene_brain_gtex$gtex.smtsd)
# Brain - Amygdala                          Brain - Anterior cingulate cortex (BA24) 
# 154                                       178 
# Brain - Caudate (basal ganglia)           Brain - Cerebellar Hemisphere 
# 246                                       217 
# Brain - Cerebellum                        Brain - Cortex 
# 245                                       257 
# Brain - Frontal Cortex (BA9)              Brain - Hippocampus 
# 209                                       203 
# Brain - Hypothalamus                      Brain - Nucleus accumbens (basal ganglia) 
# 204                                       250 
# Brain - Putamen (basal ganglia)           Brain - Spinal cord (cervical c-1) 
# 207                                       159 
# Brain - Substantia nigra 
# 141

rowData(rse_gene_brain_gtex)$ensembl <- gsub("\\..+$", "",rowData(rse_gene_brain_gtex)$gene_id)

dim(rse_gene_brain_gtex)
# [1] 63856  2670
length(unique(rowData(rse_gene_brain_gtex)$ensembl))
# [1] 63811

rownames(rse_gene_brain_gtex) <- rowData(rse_gene_brain_gtex)$ensembl

#### sce data ####
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)
# sce

rowData(sce)

common_genes <- intersect(rowData(sce)$gene_id, rowData(rse_gene_brain_gtex)$ensembl)
length(common_genes)
# [1] 34057

rownames(sce) <- rowData(sce)$gene_id

## load marker genes
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)
# marker_stats

marker_stats |>
  filter(gene %in% common_genes,
         rank_ratio <= 25) |>
  count(cellType.target)

# cellType.target     n
# <fct>           <int>
# 1 Astro              25
# 2 EndoMural          25
# 3 Micro              25
# 4 Oligo              23
# 5 OPC                23
# 6 Excit              24
# 7 Inhib              25

markers <- marker_stats |> 
  filter(gene %in% common_genes, rank_ratio <= 25) |>
  pull(gene)
               
#### Build Expression sets ####

exp_set_bulk <- ExpressionSet(assayData = assays(rse_gene_brain_gtex)$counts[markers,],
                              phenoData=AnnotatedDataFrame(
                                as.data.frame(colData(rse_gene_brain_gtex))[c("gtex.subjid")]))

exp_set_sce <- ExpressionSet(assayData = as.matrix(assays(sce)$counts[markers,]),
                             phenoData=AnnotatedDataFrame(
                               as.data.frame(colData(sce)[,c("key","Sample","BrNum", "cellType_broad_hc", "cellType_hc")])))


#### Run MuSiC ####

est_prop_music <- music_prop(bulk.mtx = assays(rse_gene_brain_gtex)$counts,
                             sc.sce = sce,
                             markers = markers,
                             clusters = "cellType_broad_hc",
                             samples = "Sample")

### run Bisque ####
exp_set_sce_temp <- exp_set_sce[,]
zero_cell_filter <- colSums(exprs(exp_set_sce_temp)) != 0
message("Exclude ",sum(!zero_cell_filter), " cells")
exp_set_sce_temp <- exp_set_sce_temp[,zero_cell_filter]
  
  est_prop_bisque <- ReferenceBasedDecomposition(bulk.eset = exp_set_bulk[m,],
                                                 sc.eset = exp_set_sce_temp,
                                                 cell.types = "cellType.Broad",
                                                 subject.names = "donor",
                                                 use.overlap = FALSE)
