
library("SingleCellExperiment")
library("MuSiC")
library("here")
library("sessioninfo")

## get args
args = commandArgs(trailingOnly=TRUE)
marker_label <- args[1]
marker_file <- NULL

if(marker_label == "FULL"){
  message("Using FULL gene-set")
} else {
  marker_file <- args[2]
  stopifnot(file.exists(marker_file))
  message("Using ", marker_label," marker genes from:", marker_file)
}

#### data output folder ####
data_dir <- here("processed-data","08_bulk_deconvolution", "02_deconvolution_MuSiC")
if(!dir.exists(data_dir)) dir.create(data_dir)


#### load data ####
## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

rownames(sce) <- rowData(sce)$gene_id

## find common genes
common_genes <- intersect(rowData(sce)$gene_id, rowData(rse_gene)$ensemblID)
length(common_genes)

if(marker_label == "FULL"){
  markers <- common_genes
} else {
  message("Input Markers:")
  markers <- scan(marker_file, what="", sep="\n")
  if(!all(markers %in% common_genes)) warning("Markers missing from common genes: ", paste(setdiff(markers, common_genes), collapse = ", "))
  markers <- intersect(markers, common_genes)
}

message(Sys.time(), " - Prep data with ", length(markers), "genes")

#### Run MuSiC ####
message(Sys.time(), " - MuSiC deconvolution")
est_prop_music <- music_prop(bulk.mtx = assays(rse_gene)$counts,
                             sc.sce = sce,
                             markers = markers,
                             clusters = "cellType_broad_hc",
                             samples = "Sample")

save(est_prop_music, file = here(data_dir,paste0("est_prop_music-", marker_label, ".Rdata")))

# slurmjobs::job_single('02_deconvolution_MuSiC_FULL', create_shell = TRUE, memory = '25G', command = "Rscript 02_deconvolution_MuSiC.R FULL")

# slurmjobs::job_single('02_deconvolution_MuSiC_MeanRatio_top25', create_shell = TRUE, memory = '25G', command = "Rscript 02_deconvolution_MuSiC.R MeanRatio_top25 ../../processed-data/08_bulk_deconvolution/markers_MeanRatio_top25.txt")
# slurmjobs::job_single('02_deconvolution_MuSiC_MeanRatio_MAD3', create_shell = TRUE, memory = '25G', command = "Rscript 02_deconvolution_MuSiC.R MeanRatio_MAD3 ../../processed-data/08_bulk_deconvolution/markers_MeanRatio_MAD3.txt")
# slurmjobs::job_single('02_deconvolution_MuSiC_MeanRatio_over2', create_shell = TRUE, memory = '25G', command = "Rscript 02_deconvolution_MuSiC.R MeanRatio_over2 ../../processed-data/08_bulk_deconvolution/markers_MeanRatio_over2.txt")

# slurmjobs::job_single('02_deconvolution_MuSiC_1vALL_top25', create_shell = TRUE, memory = '25G', command = "Rscript 02_deconvolution_MuSiC.R 1vALL_top25 ../../processed-data/08_bulk_deconvolution/markers_1vALL_top25.txt")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
