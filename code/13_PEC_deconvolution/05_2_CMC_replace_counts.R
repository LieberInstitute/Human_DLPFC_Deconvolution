library(here)
library(tidyverse)
library(SingleCellExperiment)
library(sessioninfo)
library(data.table)
library(Matrix)

rse_gene_path = here("processed-data", "rse", "rse_gene.Rdata")
sce_in_path = here(
    "processed-data", "13_PEC_deconvolution", "sce_CMC_initial.rds"
)
sce_out_path = here(
    "processed-data", "13_PEC_deconvolution", "sce_CMC.rds"
)
raw_count_dir = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version6/CMC_update_Feb2023/Count_matrices"

sce = readRDS(sce_in_path)

load(rse_gene_path, verbose = TRUE)

#   Ensembl IDs for rownames
sce$symbol = rownames(sce)
rownames(sce) <- rowData(sce)$featureid
rownames(rse_gene) = rowData(rse_gene)$ensemblID

#   Subset to genes present in bulk data
gene_match = rownames(sce) %in% rownames(rse_gene)
message('Genes from sce present in bulk data?')
table(gene_match)
sce = sce[gene_match,]

#   Raw counts are stored in separate files, one per individual. Form a list by
#   reading one by one, then combine into one sparse matrix
raw_counts_list = list()
for (ind in unique(sce$individualID)) {
    counts_path = file.path(
        raw_count_dir, paste0(ind, '-annotated_matrix.txt.gz')
    )
    a = fread(counts_path)

    #   Subset (and order) counts matrix to the genes present in the SCE
    stopifnot(all(sce$symbol %in% a$featureid))
    a = a[match(sce$symbol, a$featureid), - match('featureid', colnames(a))]

    #   Columns in the count matrix are named by their cell types. Ensure they
    #   line up, then use the colnames in the SCE object instead
    stopifnot(identical(sce$subclass[sce$individualID == ind], colnames(a)))
    colnames(a) = colnames(sce)[sce$individualID == ind]

    #   Save memory by representing as a sparse matrix
    raw_counts_list[[ind]] = as(a, "dgCMatrix")
}
raw_counts = do.call(cbind, raw_counts_list)
raw_counts = raw_counts[, colnames(sce)]

#   Replace logcounts with counts
assays(sce) = list('counts' = raw_counts)
