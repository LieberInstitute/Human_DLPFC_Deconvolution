#   This script was run to produce a lightweight SCE object,
#   subset to the top-25 mean-ratio markers and using a dense counts assay. The
#   goal is to make randomly subsetting and pseudobulking as fast as possible in
#   the following hspe and bisque scripts, which perform these operations up
#   to 1000 times

library("SingleCellExperiment")
library("tidyverse")
library("here")
library(sessioninfo)

marker_stats_path <- here(
    'processed-data', '13_PEC_deconvolution', 'CMC_marker_stats.csv'
)
sce_in_path = here("processed-data", "13_PEC_deconvolution", "sce_CMC.rds")
sce_out_path = here(
    "processed-data", "13_PEC_deconvolution", "sce_CMC_subset.rds"
)
n_markers_per_type = 25

sce = readRDS(sce_in_path)
marker_stats = read_csv(marker_stats_path, show_col_types = FALSE)

#   Take up to 25 markers per cell type, provided all ratios exceed 1
filtered_stats = marker_stats |>
    group_by(cellType.target) |>
    arrange(desc(ratio)) |>
    filter(ratio > 1) |>
    slice_head(n = n_markers_per_type) |>
    ungroup()

#   Subset to markers, cell types for which markers exist, and make counts dense
sce = sce[
    filtered_stats$gene,
    sce$subclass %in% unique(filtered_stats$cellType.target)
]
assays(sce) = list(counts = as.matrix(assays(sce)$counts))
saveRDS(sce, sce_out_path)

session_info()
