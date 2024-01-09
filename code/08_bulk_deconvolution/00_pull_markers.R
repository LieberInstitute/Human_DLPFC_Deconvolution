
library("SingleCellExperiment")
library("BisqueRNA")
library("here")
library("sessioninfo")

#### load data ####
## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

rownames(sce) <- rowData(sce)$gene_id

## find common genes
common_genes <- intersect(rowData(sce)$gene_id, rowData(rse_gene)$ensemblID)
length(common_genes)
# [1] 17804

## load marker gene data
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)
# marker_stats

# cellType.target     n
# <fct>           <int>
# 1 Astro              24
# 2 EndoMural          24
# 3 Micro              17
# 4 Oligo              25
# 5 OPC                18
# 6 Excit              23
# 7 Inhib              20

markers_top <- marker_stats |> 
  dplyr::filter(gene %in% common_genes, rank_ratio <= 25) 

markers_top |>
  dplyr::count(cellType.target)

write.csv(markers_top, file = here("processed-data","08_bulk_deconvolution", "markers_top25.csv"))

markers_top25 <- marker_stats |> 
  dplyr::filter(gene %in% common_genes, rank_ratio <= 25) |>
  dplyr::pull(gene)

save(markers_top25, file = here("processed-data","08_bulk_deconvolution", "markers_top25.Rdata"))

# slurmjobs::job_single('00_pull_markerse', create_shell = TRUE, memory = '25G', command = "Rscript 00_pull_markers.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


