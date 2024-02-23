library("SingleCellExperiment")
library("BisqueRNA")
library("here")
library("tidyverse")
library("sessioninfo")

n_cores = 8
n_runs = 1000
marker_label <- 'MeanRatio_top25'
marker_file <- here(
    'processed-data', '08_bulk_deconvolution', 'markers_MeanRatio_top25.txt'
)
sce_path = here("processed-data", "sce", "sce_DLPFC.Rdata")
bulk_path = here("processed-data","rse", "rse_gene.Rdata")
out_path <- here(
    "processed-data", "08_bulk_deconvolution", "01_deconvolution_Bisque",
    sprintf("est_prop_bisque_%s_%s_random_subsets.csv", marker_label, n_runs)
)
sce_coldata_cols = c(
    "key", "Sample", "BrNum", "cellType_broad_hc", "cellType_hc"
)

message("Using ", marker_label," marker genes from:", marker_file)

################################################################################
#   Function definitions
################################################################################

run_bisque = function(i, sce, exp_set_bulk) {
    #   Randomly subset cells of the sce such that each cell type is equally
    #   represented
    subset_keys = colData(sce) |>
        as_tibble() |>
        group_by(cellType_broad_hc) |>
        slice_sample(n = min_n_cells) |>
        pull(key)
    sce_sub = sce[, subset_keys]

    #   Convert to ExpressionSet as required for Bisque
    exp_set_sce <- ExpressionSet(
        assayData = assays(sce_sub)$counts,
        phenoData = AnnotatedDataFrame(
            as.data.frame(colData(sce_sub)[,sce_coldata_cols])
        )
    )

    #   Deconvolve bulk using Bisque
    est_prop_bisque <- ReferenceBasedDecomposition(
        bulk.eset = exp_set_bulk,
        sc.eset = exp_set_sce,
        cell.types = "cellType_broad_hc",
        subject.names = "Sample",
        use.overlap = FALSE
    )

    #   Re-format estimated bulk cell-type proportions to be long and tidy
    props_df = est_prop_bisque$bulk.props |>
        as.data.frame() |>
        rownames_to_column('cell_type') |>
        as_tibble() |>
        pivot_longer(
            cols = !matches('cell_type'),
            names_to = 'sample_id',
            values_to = 'prop'
        ) |>
        mutate(subset_run = i)
    
    return(props_df)
}

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
#   Subset and prepare data for Bisque
################################################################################

#-------------------------------------------------------------------------------
#   Bulk data
#-------------------------------------------------------------------------------

exp_set_bulk <- ExpressionSet(
    assayData = assays(rse_gene)$counts[markers,],
    phenoData = AnnotatedDataFrame(
        as.data.frame(colData(rse_gene))[c("SAMPLE_ID")]
    )
)

#-------------------------------------------------------------------------------
#   Single-cell data
#-------------------------------------------------------------------------------

#   Drop ambiguous cells
sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

#   Subset to markers and cells with at least some gene expression
nonempty_cells = colSums(assays(sce)$counts[markers,]) > 0
sce = sce[markers, nonempty_cells]
message("Excluding ", sum(!nonempty_cells), " zero-expression cells")

#   Bring into memory to greatly speed up random subsetting later
assays(sce)$counts = as.matrix(assays(sce)$counts)

#   The number of cells present for the cell type with the least cells
min_n_cells = colData(sce) |>
    as_tibble() |>
    group_by(cellType_broad_hc) |>
    summarize(n = n()) |>
    pull(n) |>
    min()

################################################################################
#   Run Bisque and export results
################################################################################

#   Run Bisque on the many different subsets in parallel
props_df_list = bplapply(
    1:n_runs,
    run_bisque,
    sce = sce,
    exp_set_bulk = exp_set_bulk,
    BPPARAM = MulticoreParam(n_cores)
)

#   Export a combined CSV of all runs
do.call(rbind, props_df_list) |>
    write_csv(out_path)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
