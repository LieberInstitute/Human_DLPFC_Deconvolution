library(getopt)
library(sessioninfo)
library("hspe")
library("SingleCellExperiment")
library("jaffelab")
library("tidyverse")
library("here")
library("spatialLIBD")

task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(task_id)

# Import command-line parameters
spec <- matrix(
    c(
        c("n_donors", "run_num"),
        c("n", "r"),
        rep("1", 2),
        rep("character", 2),
        rep("Add variable description here", 2)
    ),
    ncol = 5
)
opt <- getopt(spec)

print("Using the following parameters:")
print(opt)

marker_stats_path <- here(
    'processed-data', '13_PEC_deconvolution', 'CMC_marker_stats.csv'
)
sce_path = here("processed-data", "13_PEC_deconvolution", "sce_CMC_subset.rds")
bulk_path = here("processed-data", "rse", "rse_gene.Rdata")
out_path <- here(
    "processed-data", "13_PEC_deconvolution", "hspe_donor_subset",
    sprintf(
        "est_prop_hspe_n%s_%s_random_subsets.csv", opt$n_donors, opt$run_num
    )
)
n_markers_per_type = 25
n_total_donors = 52

dir.create(dirname(out_path), showWarnings = FALSE)

#   colnames in colData(sce)
sce_cell_type_col = "cell_type_broad"
sce_individual_col = "individualID"

#   Randomly take [opt$n_donors] donors from the set of [n_total_donors]
sce = readRDS(sce_path)
sce[[sce_cell_type_col]] = factor(sce[[sce_cell_type_col]])

#   Take up to 25 markers per cell type, provided all ratios exceed 1
filtered_stats = read_csv(marker_stats_path, show_col_types = FALSE) |>
    group_by(cellType.target) |>
    arrange(desc(ratio)) |>
    filter(ratio > 1) |>
    slice_head(n = n_markers_per_type) |>
    ungroup()

#   The SCE should be filtered to markers and contain the expected number of
#   donors
unique_donors = unique(sce$individualID)
stopifnot(setequal(filtered_stats$gene, rownames(sce)))
stopifnot(length(unique_donors) == n_total_donors)

sce_sub = sce[
    , sce[[sce_individual_col]] %in%
        unique_donors[sample(seq(n_total_donors), opt$n_donors)]
]

sce_pb = registration_pseudobulk(
    sce_sub,
    var_registration = sce_cell_type_col,
    var_sample_id = sce_individual_col
)
sce_pb[[sce_cell_type_col]] = droplevels(sce_pb[[sce_cell_type_col]])
table(sce_pb[[sce_cell_type_col]])

#   Pseudobulking performs gene filtering, which can drop markers
num_dropped_markers = length(which(!filtered_stats$gene %in% rownames(sce_pb)))
if (num_dropped_markers > 0) {
    warning(
        sprintf(
            "Pseudobulking dropped %s lowly expressed markers.",
            num_dropped_markers
        )
    )
}

#   Both subsetting and pseudobulking can drop cells and therefore cell types
dropped_types = setdiff(
    filtered_stats$cellType.target, sce_pb[[sce_cell_type_col]]
)
if (length(dropped_types) > 0) {
    warning(
        sprintf(
            "Subsetting and pseudobulking dropped %s cell types.",
            paste(dropped_types, collapse = ", ")
        )
    )
}

#   Since both markers and cells can be dropped during pseudobulking, ensure
#   we only keep markers of the remaining types that are present in the pb
#   object
filtered_stats = filtered_stats |>
    filter(
        gene %in% rownames(sce_pb),
        cellType.target %in% sce_pb[[sce_cell_type_col]]
    )

#   Form a named list with names as cell types and values as markers
marker_genes = split(filtered_stats, f = filtered_stats$cellType.target) |>
    map(~.x$gene)

#   In case markers were dropped
sce_pb = sce_pb[unlist(marker_genes),]

load(bulk_path, verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

################################################################################
#   Run hspe and export estimates to CSV
################################################################################

# we can instead explicitly pass a list of markers to hspe specifying the marker genes
# elements of the list correspond one to each cell type in the same order specified either in elements of pure_samples

pure_samples = rafalib::splitit(sce_pb[[sce_cell_type_col]])
# hspe assumes log2 transformed expressions
mixture_samples = t(assays(rse_gene)$logcounts[unlist(marker_genes),])
reference_samples = t(assays(sce_pb)$logcounts)

stopifnot(ncol(mixture_samples) == ncol(reference_samples))
stopifnot(setequal(names(pure_samples), names(marker_genes)))
marker_genes = marker_genes[names(pure_samples)]

est_prop_hspe = hspe(
    Y = mixture_samples,
    reference = reference_samples,
    pure_samples = pure_samples,
    markers = marker_genes,
    seed = task_id
)

#   Tidy up and export to CSV
est_prop_hspe$estimates |>
    as.data.frame() |>
    rownames_to_column('sample_id') |>
    as_tibble() |>
    pivot_longer(
        cols = !matches('sample_id'),
        names_to = 'cell_type',
        values_to = 'prop'
    ) |>
    mutate(subset_run = opt$run_num, num_donors = opt$n_donors) |>
    write_csv(out_path)

gc()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/
