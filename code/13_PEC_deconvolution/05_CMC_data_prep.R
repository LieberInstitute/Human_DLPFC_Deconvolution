library(zellkonverter)
library(here)
library(tidyverse)
library(SingleCellExperiment)
library(sessioninfo)

sce_ad_path = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version5/CMC_annotated.h5ad"
meta_path = '/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/sample_metadata_v5/CMC_metadata.csv'
sce_out_path = here("processed-data", "13_PEC_deconvolution", "sce_CMC.rds")
rse_gene_path = here("processed-data", "rse", "rse_gene.Rdata")

sce = readH5AD(sce_ad_path)

#   Reduce memory footprint by removing parts of the object we won't need
metadata(sce) = list()
reducedDims(sce) = list()

#   Subset to controls as established from the external metadata
meta = read_csv(meta_path, show_col_types = FALSE)
stopifnot(setequal(meta$individualID, sce$individualID))
sce$primaryDiagnosis = meta$primaryDiagnosis[
    match(sce$individualID, meta$individualID)
]
sce = sce[,sce$primaryDiagnosis == "control"]

load(rse_gene_path, verbose = TRUE)

#   Ensembl IDs for rownames
rownames(sce) <- rowData(sce)$featureid
rownames(rse_gene) = rowData(rse_gene)$ensemblID

#   Subset to genes present in bulk data
gene_match = rownames(sce) %in% rownames(rse_gene)
message('Genes from sce present in bulk data?')
table(gene_match)
sce = sce[gene_match,]

saveRDS(sce, sce_out_path)

session_info()
