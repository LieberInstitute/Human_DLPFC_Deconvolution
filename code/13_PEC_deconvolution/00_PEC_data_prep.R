
# pec_data = list(CMC = c(33792, 456460),
#                 DevBrain = c(33201, 102936),
#                 IsoHuB = c(33339, 29675),
#                 LIBD = c(33783, 52214),
#                 MultiomeBrain = c(33818, 141277),
#                 PTSDBrainomics = c(33877, 198572),
#                 SZBDMulti = c(34361, 603209),
#                 `UCLA-ASD` = c(34180, 448524))
# 
# n_cell <- map_int(pec_data,~pluck(.x, 2))
# 
# sort(n_cell)

library("zellkonverter")
library("SingleCellExperiment")
# library("spatialLIBD")
# library("jaffelab")
library("here")
library("sessioninfo")

## rawdata path 
# "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version5"

filepath <- "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version5/PTSDBrainomics_annotated.h5ad"

stopifnot(file.exists(filepath))

message(Sys.time(), " - Reading data from: ", filepath)
sce <- readH5AD(file = filepath)

message("\nSCE Dimesions:")
dim(sce)

print(colnames(colData(sce)))
message("Cell Types:")
## must be syntactically valid
colData(sce)$cellType <- as.factor(make.names(colData(sce)$subclass))
table(sce$cellType)

# identify annotation/cluster labels
rowData(sce)$gene_name <- rownames(sce) # save gene name as column of rowData
rownames(sce) <- rowData(sce)$featureid # have to make row names of object the ensembl id instead of gene names

## Logcounts
# default “X” contain the log-normalized counts
message(Sys.time(), " revert to counts")

## check for all 0s (just first 100 cols for mem)
stopifnot(any(assays(sce)$X[, 1:100] != 0))

counts(sce) <- assays(sce)$X # change to normalized counts
# counts(sce)[counts(sce) != 0] <- (2^counts(sce)[counts(sce) != 0])-1 # Replace just non-zero values
counts(sce)@x <- 2^(counts(sce)@x) - 1 ## remove log2(counts + 1)

table(sce$cellType)

save(sce, filename = here("processed-data", "13_PEC_deconvolution", "sce_PTSDBrainomics.Rdata"))


# slurmjobs::job_single(name = "00_PEC_data_prep", memory = "5G", cores = 1, create_shell = TRUE, command = "Rscript 03_get_est_prop.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
