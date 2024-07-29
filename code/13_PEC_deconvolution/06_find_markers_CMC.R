library("DeconvoBuddies")
library("SingleCellExperiment")
library("tidyverse")
library("here")
library("sessioninfo")

sce_path = here("processed-data", "13_PEC_deconvolution", "sce_CMC.rds")
stats_out_path = here(
    "processed-data", "13_PEC_deconvolution", "CMC_marker_stats.csv"
)
markers_out_path = here(
    "processed-data", "13_PEC_deconvolution", "CMC_markers_MeanRatio_top25.txt"
)
ct_table_path = here(
    "processed-data", "13_PEC_deconvolution", "PEC_cell_types.csv"
)
n_markers_per_type = 25

sce = readRDS(sce_path)

#   Remove counts in this script to reduce memory usage
assays(sce)$counts = NULL
gc()

#   Add cell type at broad resolution. Drop certain broad cell types as
#   consistent with the other PEC analysis
ct_table = read_csv(ct_table_path, show_col_types = FALSE)
stopifnot(all(sce$subclass %in% ct_table$subclass))
sce$cell_type_broad = ct_table$cellType_broad[
    match(sce$subclass, ct_table$subclass)
]
sce = sce[, sce$cell_type_broad != "drop"]

#   Find markers at broad resolution
mean_ratio <- get_mean_ratio2(
    sce, cellType_col = "cell_type_broad", assay_name = "logcounts",
    add_symbol = TRUE
)

write_csv(mean_ratio, stats_out_path)

#   Take up to 25 markers per cell type, provided all ratios exceed 1
filtered_stats = mean_ratio |>
    group_by(cellType.target) |>
    arrange(desc(ratio)) |>
    filter(ratio > 1) |>
    slice_head(n = n_markers_per_type) |>
    ungroup()

message(
    sprintf(
        "Number of markers found by cell type (%s total):", nrow(filtered_stats)
    )
)
table(filtered_stats$cellType.target)

missing_types = setdiff(sce$cell_type_broad, filtered_stats$cellType.target)
if (length(missing_types) > 0) {
    warning(
        sprintf(
            "No markers were found for %s cell type(s)",
            paste(missing_types, collapse = ", ")
        )
    )
}

#   Write one marker per line to a text file
writeLines(filtered_stats$gene, markers_out_path)

session_info()
