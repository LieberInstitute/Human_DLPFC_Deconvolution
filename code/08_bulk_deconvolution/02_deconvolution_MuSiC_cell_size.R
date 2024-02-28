marker_label <- 'MeanRatio_top25'
marker_file <- here(
    'processed-data', '08_bulk_deconvolution', 'markers_MeanRatio_top25.txt'
)
sce_path = here("processed-data", "sce", "sce_DLPFC.Rdata")
bulk_path = here("processed-data","rse", "rse_gene.Rdata")
out_path <- here(
    "processed-data", "08_bulk_deconvolution", "02_deconvolution_MuSiC",
    sprintf("est_prop_bisque_%s_%s_random_subsets.csv", marker_label, n_runs)
)


message("Using ", marker_label," marker genes from:", marker_file)

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
#   Run MuSiC and save
################################################################################

message(Sys.time(), " - MuSiC deconvolution")
est_prop_music <- music_prop(
    bulk.mtx = assays(rse_gene)$counts,
    sc.sce = sce,
    markers = markers,
    clusters = "cellType_broad_hc",
    samples = "Sample",
    cell_size = NULL
)

save(est_prop_music, file = here(data_dir,paste0("est_prop_music-", marker_label, ".Rdata")))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
