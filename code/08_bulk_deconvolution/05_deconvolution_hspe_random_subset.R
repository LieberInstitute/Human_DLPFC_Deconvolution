library("hspe")
library("SingleCellExperiment")
library("jaffelab")
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")

marker_label <- "MeanRatio_top25"
sce_path = here("processed-data", "sce", "sce_DLPFC.Rdata")
bulk_path = here("processed-data", "rse", "rse_gene.Rdata")
marker_stats_path = here(
    "processed-data", "06_marker_genes", "03_find_markers_broad",
    "marker_stats_broad.Rdata"
)

#### data output folder ####
data_dir <- here("processed-data","08_bulk_deconvolution", "05_deconvolution_hspe")
if(!dir.exists(data_dir)) dir.create(data_dir)

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

#   The number of cells present for the cell type with the least cells
min_n_cells = colData(sce) |>
    as_tibble() |>
    group_by(cellType_broad_hc) |>
    summarize(n = n()) |>
    pull(n) |>
    min()

#   Randomly subset cells of the sce such that each cell type is equally
#   represented
subset_keys = colData(sce) |>
    as_tibble() |>
    group_by(cellType_broad_hc) |>
    slice_sample(n = min_n_cells) |>
    pull(key)
sce_sub = sce[, subset_keys]

sce_pb = registration_pseudobulk(
    sce_sub, var_registration = "cellType_broad_hc", var_sample_id = "Sample"
)
table(sce_pb$cellType_broad_hc)

## find common genes
common_genes <- intersect(rowData(sce_pb)$gene_id, rowData(rse_gene)$ensemblID)
message("common genes: ", length(common_genes))
# [1] 13311 

# we can instead explicitly pass a list of markers to hspe specifying the marker genes
# elements of the list correspond one to each cell type in the same order specified either in elements of pure_samples

pure_samples = rafalib::splitit(sce_pb$cellType_broad_hc)
# hspe assumes log2 transformed expressions
mixture_samples = t(assays(rse_gene)$logcounts[common_genes,])
reference_samples = t(assays(sce_pb)$logcounts[common_genes,])

stopifnot(ncol(mixture_samples) == ncol(reference_samples))

est_prop_hspe <- NULL

## load marker gene data
load(marker_stats_path, verbose = TRUE) 
marker_tab <- marker_stats |> 
    dplyr::filter(gene %in% common_genes, rank_ratio <= 25)

marker_genes <- purrr::map(
    rafalib::splitit(marker_tab$cellType.target), ~marker_tab$gene[.x]
)
marker_genes <- marker_genes[names(pure_samples)]
  
message(Sys.time(), "- hspe w/ markers ", marker_label)
purrr::map_int(marker_genes, length)

est_prop_hspe = hspe(
    Y = mixture_samples,
    reference = reference_samples,
    pure_samples = pure_samples,
    markers = marker_genes
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
