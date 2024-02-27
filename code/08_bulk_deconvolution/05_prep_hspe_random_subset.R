#   This script was interactively run to produce a lightweight SCE object,
#   subset to the top-25 mean-ratio markers and using a dense counts assay. The
#   goal is to make randomly subsetting and pseudobulking as fast as possible in
#   the '05_deconvolution_hpse_random_subset.*' array job, since it gets
#   performed 1000 times

library("SingleCellExperiment")
library("jaffelab")
library("tidyverse")
library("here")

marker_label <- "MeanRatio_top25"
sce_path = here("processed-data", "sce", "sce_DLPFC.Rdata")
bulk_path = here("processed-data", "rse", "rse_gene.Rdata")
marker_stats_path = here(
    "processed-data", "06_marker_genes", "03_find_markers_broad",
    "marker_stats_broad.Rdata"
)
sce_out_path = here(
    "processed-data", "08_bulk_deconvolution", "05_deconvolution_hspe",
    "sce_filtered.rds"
)

#### load data ####
## load bulk data
load(bulk_path, verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
dim(rse_gene)
# [1] 21745   110

load(sce_path, verbose = TRUE)
rownames(sce) <- rowData(sce)$gene_id
colnames(sce) = sce$key

#   Drop ambiguous cells
sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

## find common genes
common_genes <- intersect(rowData(sce)$gene_id, rowData(rse_gene)$ensemblID)
message("common genes: ", length(common_genes))
# common genes: 17804 

## load marker gene data
load(marker_stats_path, verbose = TRUE) 
marker_tab <- marker_stats |> 
    filter(gene %in% common_genes, rank_ratio <= 25)

marker_genes <- purrr::map(
    rafalib::splitit(marker_tab$cellType.target), ~marker_tab$gene[.x]
)
marker_genes <- marker_genes[levels(sce$cellType_broad_hc)]

#   Subset to markers, drop unused assays and pull counts assay into memory.
#   Drop reducedDims to keep the object as small as possible
sce = sce[unlist(marker_genes), ]
assays(sce) = list(counts = as.matrix(assays(sce)$counts))
reducedDim(sce) = NULL

saveRDS(sce, sce_out_path)
