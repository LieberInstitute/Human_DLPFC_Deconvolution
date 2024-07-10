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
n_markers_per_type = 25

sce = readRDS(sce_path)

mean_ratio <- get_mean_ratio2(
    sce, cellType_col = "subclass", assay_name = "X",
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

#   Write one marker per line to a text file
writeLines(filtered_stats$gene, markers_out_path)

session_info()
