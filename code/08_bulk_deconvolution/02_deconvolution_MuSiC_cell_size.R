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
halo_path = here("processed-data", "03_HALO", "12_HALO_cell_size", "MuSiC_cell_size.csv")
out_path = here(
    "processed-data", "08_bulk_deconvolution", "02_deconvolution_MuSiC",
    sprintf("est_prop_cell_size_%s.csv", cell_size_opt)
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


################################################################################
#   Load HALO "cell sizes" for each cell type
################################################################################

#   Keep only cell types that match the single-cell data; treat "OligoOPC"
#   as if it's just oligos (later we'll use define OPCs to have the same cell
#   sizes as oligos, since this category was merged in the HALO data and not
#   recoverable into distinct categories)

MuSiC_cell_size <- read.csv(halo_path) |>
    mutate(cell_type = ifelse(cell_type == "OligoOPC", "Oligo", cell_type)) |>
    filter(cell_type != "Other")

stopifnot(all(halo_all$cell_type %in% unique(sce$cellType_broad_hc)))

#   For each cell type, compute an average "cell size" based on some combination
#   of nuclear area and ATK3 puncta
if (cell_size_opt == "nuc_area") {
    cell_size_df = MuSiC_cell_size |>
        select(cell_type, cell_size = nuc_area)
} else if (cell_size_opt == "akt3") {
  cell_size_df = MuSiC_cell_size |>
    select(cell_type, cell_size = akt3)
} else if (cell_size_opt == "nuc_area_akt3") {
  cell_size_df = MuSiC_cell_size |>
    select(cell_type, cell_size = nuc_area_akt3)
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

#   Tidy up and export to CSV
est_prop_music$Est.prop.weighted |>
    as.data.frame() |>
    rownames_to_column('sample_id') |>
    as_tibble() |>
    pivot_longer(
        cols = !matches('sample_id'),
        names_to = 'cell_type',
        values_to = 'prop'
    ) |>
    mutate(cell_size_opt = {{ cell_size_opt }}) |>
    write_csv(out_path)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
