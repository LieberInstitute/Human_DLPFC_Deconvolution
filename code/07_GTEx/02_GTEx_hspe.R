
library("recount3")
library("SingleCellExperiment")
library("hspe")
library("tidyverse")
library("here")
library("sessioninfo")
# library(DeconvoBuddies)

data_dir <- here("processed-data",  "07_GTEx", "02_GTEx_hspe")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

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

## Compute logtransformed counts consistent with paired data
message(Sys.time(), " - get logcounts")
# assay(rse_gene_brain_gtex, "counts") <- transform_counts(rse_gene_brain_gtex)
assay(rse_gene_brain_gtex, "counts") <- compute_read_counts(rse_gene_brain_gtex)
assay(rse_gene_brain_gtex, "RPKM") <- recount::getRPKM(rse_gene_brain_gtex)
assay(rse_gene_brain_gtex, "logcounts") <- log2(assays(rse_gene_brain_gtex)$RPKM +1)

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

#### pseudobulk sce data ####
sce_pb <- readRDS(here("processed-data", "sce","sce_broad_pseudobulk.rds"))
table(sce_pb$cellType_broad_hc)

## use ensemblIDs
rownames(sce_pb) <- rowData(sce_pb)$gene_id

## find common genes
common_genes <- intersect(rowData(sce_pb)$gene_id, rowData(rse_gene_brain_gtex)$ensembl)
message("common genes: ", length(common_genes))
# 14809

#### prep hspe input ####
pure_samples = rafalib::splitit(sce_pb$cellType_broad_hc)
# hspe assumes log2 transformed expressions
mixture_samples = t(assays(rse_gene_brain_gtex)$logcounts[common_genes,])
reference_samples = t(assays(sce_pb)$logcounts[common_genes,])

stopifnot(ncol(mixture_samples) == ncol(reference_samples))

#### marker genes ####
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)

marker_tab <- marker_stats |> 
  dplyr::filter(gene %in% common_genes, 
                rank_ratio <= 25)
               
marker_genes <- purrr::map(rafalib::splitit(marker_tab$cellType.target), ~marker_tab$gene[.x])

marker_genes <- marker_genes[names(pure_samples)]

message(Sys.time(), "- hspe w/ Mean Ratio top25 markers")
purrr::map_int(marker_genes, length)

est_prop_hspe = hspe(Y = mixture_samples,
                     reference = reference_samples,
                     pure_samples = pure_samples,
                     markers = marker_genes,
                     seed = 10524)


save(est_prop_hspe, file = here(data_dir, "GTEx_est_prop_hspe_MeanRatio_top25_rc.Rdata"))

# slurmjobs::job_single('02_GTEx_hspe_rc', create_shell = TRUE, memory = '25G', command = "Rscript 02_GTEx_hspe.R")
# slurmjobs::job_single('02_GTEx_hspe', create_shell = TRUE, memory = '25G', command = "Rscript 02_GTEx_hspe.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

