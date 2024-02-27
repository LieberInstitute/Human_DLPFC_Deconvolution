library("hspe")
library("SingleCellExperiment")
library("jaffelab")
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")

task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(task_id)

marker_label <- "MeanRatio_top25"
sce_path = here(
    "processed-data", "08_bulk_deconvolution", "05_deconvolution_hspe",
    "sce_filtered.rds"
)
bulk_path = here("processed-data", "rse", "rse_gene.Rdata")
marker_stats_path = here(
    "processed-data", "06_marker_genes", "03_find_markers_broad",
    "marker_stats_broad.Rdata"
)
out_path = here(
    "processed-data", "08_bulk_deconvolution", "05_deconvolution_hspe",
    "random_subset_runs", sprintf("est_prop_%s.csv", task_id)
)

dir.create(dirname(out_path), showWarnings = FALSE)

################################################################################
#   Load bulk data, single-cell data, and markers
################################################################################

#### load data ####
## load bulk data
load(bulk_path, verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
dim(rse_gene)
# [1] 21745   110

sce = readRDS(sce_path)
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

################################################################################
#   Randomly subset single-cell data then pseudobulk
################################################################################

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

#   Pseudobulking performs gene filtering, but markers should all be
#   well-expressed
stopifnot(all(unlist(marker_genes) %in% rownames(sce_pb)))

################################################################################
#   Run hspe and export estimates to CSV
################################################################################

# we can instead explicitly pass a list of markers to hspe specifying the marker genes
# elements of the list correspond one to each cell type in the same order specified either in elements of pure_samples

pure_samples = rafalib::splitit(sce_pb$cellType_broad_hc)
# hspe assumes log2 transformed expressions
mixture_samples = t(assays(rse_gene)$logcounts[unlist(marker_genes),])
reference_samples = t(assays(sce_pb)$logcounts)

stopifnot(ncol(mixture_samples) == ncol(reference_samples))
  
message(Sys.time(), "- hspe w/ markers ", marker_label)
purrr::map_int(marker_genes, length)

est_prop_hspe = hspe(
    Y = mixture_samples,
    reference = reference_samples,
    pure_samples = pure_samples,
    markers = marker_genes,
    seed = task_id
)

#   Tidy up and export to CSV
est_prop_hspe$estimates |>
    as.data.frame() |>
    rownames_to_column('sample_id') |>
    as_tibble() |>
    pivot_longer(
        cols = !matches('sample_id'),
        names_to = 'cell_type',
        values_to = 'prop'
    ) |>
    mutate(subset_run = task_id) |>
    write_csv(out_path)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
