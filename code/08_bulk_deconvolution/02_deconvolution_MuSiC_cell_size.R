library("SingleCellExperiment")
library("MuSiC")
library("here")
library("sessioninfo")
library("HDF5Array")
library("tidyverse")

#   Later controls the argument to 'cell_size' parameter for music_prop
cell_size_opt = c("nuc_area", "akt3", "nuc_area_akt3")[
    as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
]

marker_label <- 'MeanRatio_top25'
marker_file <- here(
    'processed-data', '08_bulk_deconvolution', 'markers_MeanRatio_top25.txt'
)
sce_path = here("processed-data", "sce", "sce_DLPFC.Rdata")
bulk_path = here("processed-data","rse", "rse_gene.Rdata")
halo_path = here("processed-data", "03_HALO", "halo_all.Rdata")
out_path <- here(
    "processed-data", "08_bulk_deconvolution", "02_deconvolution_MuSiC",
    sprintf("est_prop_bisque_%s_%s_random_subsets.csv", marker_label, n_runs)
)

message("Using ", marker_label," marker genes from: ", marker_file)
message(sprintf("Using cell size option '%s'.", cell_size_opt))

################################################################################
#   Load bulk data, single-cell data, and top-25 mean-ratio markers
################################################################################

## load bulk data
load(bulk_path, verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
dim(rse_gene)
# [1] 21745   110

## sce data
load(sce_path, verbose = TRUE)
rownames(sce) <- rowData(sce)$gene_id
colnames(sce) = sce$key

#   Drop Ambiguous cell type
sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

## find common genes
common_genes <- intersect(rowData(sce)$gene_id, rowData(rse_gene)$ensemblID)
length(common_genes)
# [1] 17804

message("Input Markers:")
markers <- scan(marker_file, what="", sep="\n")
if(!all(markers %in% common_genes)) {
    warning(
        "Markers missing from common genes: ",
        paste(setdiff(markers, common_genes), collapse = ", ")
    )
}
markers <- intersect(markers, common_genes)

load(halo_path, verbose = TRUE)

################################################################################
#   Load HALO data and produce a table of "cell sizes" for each cell type
################################################################################

#   Keep only cell types that match the single-cell data; treat "OligoOPC"
#   as if it's just oligos (later we'll use define OPCs to have the same cell
#   sizes as oligos, since this category was merged in the HALO data and not
#   recoverable into distinct categories)
halo_all = halo_all |>
    mutate(cell_type = ifelse(cell_type == "OligoOPC", "Oligo", cell_type)) |>
    filter(cell_type != "Other")
stopifnot(all(halo_all$cell_type %in% unique(sce$cellType_broad_hc)))

#   For each cell type, compute an average "cell size" based on some combination
#   of nuclear area and ATK3 puncta
if (cell_size_opt == "nuc_area") {
    cell_size_df = halo_all |>
        group_by(cell_type) |>
        summarize(cell_size = mean(Nucleus_Area))
} else if (cell_size_opt == "akt3") {
    cell_size_df = halo_all |>
        group_by(cell_type) |>
        summarize(cell_size = mean(DAPI_AKT3))
} else if (cell_size_opt == "nuc_area_akt3") {
    cell_size_df = halo_all |>
        group_by(cell_type) |>
        summarize(cell_size = mean(Nucleus_Area * DAPI_AKT3))
} else {
    stop("Invalid cell-size option.")
}

#   Add a row for OPC that matches the row for Oligo
cell_size_df = rbind(
    cell_size_df,
    tibble(
        cell_type = 'OPC',
        cell_size = cell_size_df |>
            filter(cell_type == "Oligo") |>
            pull(cell_size)
    )
) |> as.data.frame()

################################################################################
#   Run MuSiC and save
################################################################################

message(Sys.time(), " - MuSiC deconvolution")
est_prop_music <- music_prop(
    bulk.mtx = assays(rse_gene)$counts,
    sc.sce = sce,
    markers = markers,
    clusters = "cellType_broad_hc",
    samples = "Sample",
    cell_size = cell_size_df
)

save(est_prop_music, file = here(data_dir,paste0("est_prop_music-", marker_label, ".Rdata")))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
