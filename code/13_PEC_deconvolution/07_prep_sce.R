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
ct_table_path = here(
    "processed-data", "13_PEC_deconvolution", "PEC_cell_types.csv"
)
n_markers_per_type = 25

sce = readRDS(sce_in_path)

#   Remove logcounts, which are large and no longer needed
assays(sce)$logcounts = NULL
gc()

#   Add cell type at broad resolution. Drop certain broad cell types as
#   consistent with the other PEC analysis
ct_table = read_csv(ct_table_path, show_col_types = FALSE)
stopifnot(all(sce$subclass %in% ct_table$subclass))
sce$cell_type_broad = ct_table$cellType_broad[
    match(sce$subclass, ct_table$subclass)
]
sce = sce[, sce$cell_type_broad != "drop"]

marker_stats = read_csv(marker_stats_path, show_col_types = FALSE)

#   Take up to 25 markers per cell type, provided all ratios exceed 1
filtered_stats = marker_stats |>
    group_by(cellType.target) |>
    arrange(desc(ratio)) |>
    filter(ratio > 1) |>
    slice_head(n = n_markers_per_type) |>
    ungroup()

#   Subset to markers and make counts dense
sce = sce[filtered_stats$gene,]
assays(sce) = list(counts = as.matrix(assays(sce)$counts))
saveRDS(sce, sce_out_path)

session_info()
