library("SingleCellExperiment")
library("BisqueRNA")
library("here")
library("tidyverse")
library("sessioninfo")
library("BiocParallel")

n_runs = 100
n_total_donors = 52
n_donors = round(1.4**(1:10) * n_total_donors / 1.4**10)[
    as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
]

marker_file <- here(
    "processed-data", "13_PEC_deconvolution", "CMC_markers_MeanRatio_top25.txt"
)
sce_path = here("processed-data", "13_PEC_deconvolution", "sce_CMC_subset.rds")
bulk_path = here("processed-data","rse", "rse_gene.Rdata")
out_path <- here(
    "processed-data", "13_PEC_deconvolution", "bisque_donor_subset",
    sprintf("est_prop_bisque_n%s_%s_random_subsets.csv", n_donors, n_runs)
)

sce_assay_name = "counts"

#   colnames in colData(sce)
sce_cell_type_col = "cell_type_broad"
sce_individual_col = "individualID"

n_cores = as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))

dir.create(dirname(out_path), showWarnings = FALSE)

################################################################################
#   Function definitions
################################################################################

run_bisque = function(i, sce, exp_set_bulk, unique_donors) {
    #   Randomly take [n_donors] donors from the set of [n_total_donors]
    sce_sub = sce[
        , sce[[sce_individual_col]] %in%
            unique_donors[sample(seq(n_total_donors), n_donors)]
    ]

    #   Convert to ExpressionSet as required for Bisque
    exp_set_sce <- ExpressionSet(
        assayData = assays(sce_sub)[[sce_assay_name]],
        phenoData = AnnotatedDataFrame(
            as.data.frame(
                colData(sce_sub)[,c(sce_cell_type_col, sce_individual_col)]
            )
        )
    )

    #   Deconvolve bulk using Bisque
    est_prop_bisque <- ReferenceBasedDecomposition(
        bulk.eset = exp_set_bulk,
        sc.eset = exp_set_sce,
        cell.types = sce_cell_type_col,
        subject.names = sce_individual_col,
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
        mutate(subset_run = i, num_donors = n_donors)

    return(props_df)
}

################################################################################
#   Load bulk data, single-cell data, and top-25 mean-ratio markers
################################################################################

markers <- readLines(marker_file)

#-------------------------------------------------------------------------------
#   Bulk data
#-------------------------------------------------------------------------------

load(bulk_path, verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
stopifnot(all(markers %in% rownames(rse_gene)))

exp_set_bulk <- ExpressionSet(
    assayData = assays(rse_gene)$counts[markers,],
    phenoData = AnnotatedDataFrame(
        as.data.frame(colData(rse_gene))[c("SAMPLE_ID")]
    )
)

#-------------------------------------------------------------------------------
#   Single-cell data
#-------------------------------------------------------------------------------

## sce data and markers
sce = readRDS(sce_path)

unique_donors = unique(sce$individualID)
stopifnot(length(unique_donors) == n_total_donors)
stopifnot(setequal(markers, rownames(sce)))

#   Subset to cells with at least some gene expression
nonempty_cells = colSums(assays(sce)[[sce_assay_name]]) > 0
sce = sce[, nonempty_cells]
message("Excluding ", sum(!nonempty_cells), " zero-expression cells")

################################################################################
#   Run Bisque and export results
################################################################################

#   Run Bisque on the many different subsets in parallel
props_df_list = bplapply(
    1:n_runs,
    run_bisque,
    sce = sce,
    exp_set_bulk = exp_set_bulk,
    unique_donors = unique_donors,
    BPPARAM = MulticoreParam(
        n_cores, RNGseed = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
    )
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
