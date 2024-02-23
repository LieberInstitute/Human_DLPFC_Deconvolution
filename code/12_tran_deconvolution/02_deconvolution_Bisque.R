
library("SingleCellExperiment")
library("BisqueRNA")
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
data_dir <- here("processed-data","12_tran_deconvolution", "02_deconvolution_Bisque")
if(!dir.exists(data_dir)) dir.create(data_dir)

#### load data ####
## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce data
load(here("processed-data", "12_tran_deconvolution", "sce.dlpfc.tran.Rdata"), verbose = TRUE)
table(sce$cellType_broad)

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

message(Sys.time(), " - Prep data with ", length(markers), " genes")

#### Build Expression sets ####
message(Sys.time(), " - Prep Bisque Data")
exp_set_bulk <- ExpressionSet(assayData = assays(rse_gene)$counts[markers,],
                              phenoData=AnnotatedDataFrame(
                                as.data.frame(colData(rse_gene))[c("SAMPLE_ID")]))

exp_set_sce <- ExpressionSet(assayData = as.matrix(assays(sce)$counts[markers,]),
                             phenoData=AnnotatedDataFrame(
                               as.data.frame(colData(sce)[,c("donor", "cellType_broad")])))

### run Bisque ####
exp_set_sce_temp <- exp_set_sce[markers,]
zero_cell_filter <- colSums(exprs(exp_set_sce_temp)) != 0
message("Exclude ",sum(!zero_cell_filter), " cells")
# Exclude 33 cells
exp_set_sce_temp <- exp_set_sce_temp[,zero_cell_filter]

message(Sys.time(), " - Bisque deconvolution")
est_prop_bisque <- ReferenceBasedDecomposition(bulk.eset = exp_set_bulk[markers,],
                                               sc.eset = exp_set_sce_temp,
                                               cell.types = "cellType_broad",
                                               subject.names = "donor",
                                               use.overlap = FALSE)

save(est_prop_bisque, file = here(data_dir, paste0("Tran_est_prop_bisque-",marker_label,".Rdata")))

# slurmjobs::job_single('01_deconvolution_Bisque_FULL', create_shell = TRUE, memory = '100G', command = "Rscript 01_deconvolution_Bisque.R FULL")

# slurmjobs::job_single('01_deconvolution_Bisque_MeanRatio_top25', create_shell = TRUE, memory = '25G', command = "Rscript 01_deconvolution_Bisque.R MeanRatio_top25 ../../processed-data/08_bulk_deconvolution/markers_MeanRatio_top25.txt")
# slurmjobs::job_single('01_deconvolution_Bisque_MeanRatio_MAD3', create_shell = TRUE, memory = '10G', command = "Rscript 01_deconvolution_Bisque.R MeanRatio_MAD3 ../../processed-data/08_bulk_deconvolution/markers_MeanRatio_MAD3.txt")
# slurmjobs::job_single('01_deconvolution_Bisque_MeanRatio_over2', create_shell = TRUE, memory = '10G', command = "Rscript 01_deconvolution_Bisque.R MeanRatio_over2 ../../processed-data/08_bulk_deconvolution/markers_MeanRatio_over2.txt")

# slurmjobs::job_single('01_deconvolution_Bisque_1vALL_top25', create_shell = TRUE, memory = '25G', command = "Rscript 01_deconvolution_Bisque.R 1vALL_top25 ../../processed-data/08_bulk_deconvolution/markers_1vALL_top25.txt")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


