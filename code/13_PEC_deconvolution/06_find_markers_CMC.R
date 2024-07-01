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
    sce, cellType_col = "cellType_broad", assay_name = "counts",
    add_symbol = TRUE
)

write_csv(mean_ratio, stats_out_path)

#   Write one marker per line to a text file
mean_ratio |>
    filter(rank_ratio <= n_markers_per_type) |>
    pull(gene) |>
    writeLines(markers_out_path)

session_info()
