library("SingleCellExperiment")
library("BisqueRNA")
library("here")
library("tidyverse")
library("sessioninfo")
library("BiocParallel")

n_runs = 100
n_total_donors = 104
n_donors = round(1.5**(1:10) * n_total_donors / 1.5**10)[
    as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
]

marker_file <- here(
    "processed-data", "13_PEC_deconvolution", "CMC_markers_MeanRatio_top25.txt"
)
sce_path = here("processed-data", "13_PEC_deconvolution", "sce_CMC.rds")
bulk_path = here("processed-data","rse", "rse_gene.Rdata")
out_path <- here(
    "processed-data", "13_PEC_deconvolution", "bisque_donor_subset",
    sprintf("est_prop_bisque_n%s_%s_random_subsets.csv", n_donors, n_runs)
)
sce_coldata_cols = c(
    "key", "Sample", "BrNum", "cellType_broad_hc", "cellType_hc"
)

n_cores = as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))

dir.create(dirname(out_path), showWarnings = FALSE)

################################################################################
#   Function definitions
################################################################################

unique_donors = unique(sce$donor)
stopifnot(length(unique_donors) == n_total_donors)

run_bisque = function(i, sce, exp_set_bulk) {
    #   Randomly take [n_donors] donors from the set of [n_total_donors]
    sce_sub = sce[
        , sce$donor %in% unique_donors[sample(seq(n_total_donors), n_donors)]
    ]

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
        mutate(subset_run = i, num_donors = n_donors)
    
    return(props_df)
}

################################################################################
#   Load bulk data, single-cell data, and top-25 mean-ratio markers
################################################################################

## load bulk data
load(bulk_path, verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce data and markers
sce = readRDS(sce_path)
markers <- readLines(marker_file)
