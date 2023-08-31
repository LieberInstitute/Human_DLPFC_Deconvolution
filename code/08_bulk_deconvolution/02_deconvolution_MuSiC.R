
library("SingleCellExperiment")
library("MuSiC")
library("here")
library("sessioninfo")

#### load data ####
## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

rownames(sce) <- rowData(sce)$gene_id

## find common genes
common_genes <- intersect(rowData(sce)$gene_id, rowData(rse_gene)$ensemblID)
length(common_genes)

## load marker gene data
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)
# marker_stats

marker_stats |>
  dplyr::filter(gene %in% common_genes,
         rank_ratio <= 25) |>
  dplyr::count(cellType.target)


markers <- marker_stats |> 
  dplyr::filter(gene %in% common_genes, rank_ratio <= 25) |>
  dplyr::pull(gene)

#### Run MuSiC ####
message(Sys.time(), " - MuSiC deconvolution")
est_prop_music <- music_prop(bulk.mtx = assays(rse_gene)$counts,
                             sc.sce = sce,
                             markers = markers,
                             clusters = "cellType_broad_hc",
                             samples = "Sample")

save(est_prop_music, file = here("processed-data","08_bulk_deconvolution","est_prop_music.Rdata"))

# sgejobs::job_single('02_deconvolution_MuSiC', create_shell = TRUE, memory = '25G', command = "Rscript 02_deconvolution_MuSiC.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
