library(zellkonverter)
library(here)
library(tidyverse)
library(SingleCellExperiment)
library(sessioninfo)

sce_ad_path = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version6/CMC_update_Feb2023/CMC_annotated.h5ad"
meta_path = '/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/sample_metadata_v5/CMC_metadata.csv'
sce_out_path = here("processed-data", "13_PEC_deconvolution", "sce_CMC_initial.rds")

sce = readH5AD(sce_ad_path)

#   Reduce memory footprint by removing parts of the object we won't need
metadata(sce) = list()
reducedDims(sce) = list()
rowData(sce)$varm = NULL

#   Subset to controls as established from the external metadata
meta = read_csv(meta_path, show_col_types = FALSE)
stopifnot(all(sce$individualID %in% meta$individualID))
sce$primaryDiagnosis = meta$primaryDiagnosis[
    match(sce$individualID, meta$individualID)
]
sce = sce[,sce$primaryDiagnosis == "control"]

saveRDS(sce, sce_out_path)

session_info()
