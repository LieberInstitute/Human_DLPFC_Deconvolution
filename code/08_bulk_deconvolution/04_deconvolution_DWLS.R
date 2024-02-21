
# install.packages("DWLS") ## Use cran version

library("DWLS")
library("SingleCellExperiment")
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
data_dir <- here("processed-data","08_bulk_deconvolution", "04_deoncvolution_DWLS")
if(!dir.exists(data_dir)) dir.create(data_dir)

## for marker data from subsets
data_dir_marker <- here("processed-data","08_bulk_deconvolution", "04_deoncvolution_DWLS", marker_label)
if(!dir.exists(data_dir_marker)) dir.create(data_dir_marker)

#### load data ####
## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

## use ensemblIDs
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## load sce data, rm ambiguous
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)
sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]

## use ensemblIDs
rownames(sce) <- rowData(sce)$gene_id

sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)
table(sce$cellType_broad_hc)
# Oligo Excit Inhib 
# 237   521   242 

## find common genes - subset to marker list ##
common_genes <- intersect(rowData(sce)$gene_id, rowData(rse_gene)$ensemblID)
length(common_genes)
# [1] 17804

if(marker_label == "FULL"){
  markers <- common_genes
} else {
  message("Input Markers:")
  markers <- scan(marker_file, what="", sep="\n")
  if(!all(markers %in% common_genes)) warning("Markers missing from common genes: ", paste(setdiff(markers, common_genes), collapse = ", "))
  markers <- intersect(markers, common_genes)
}

message(Sys.time(), " - Prep data with ", length(markers), " genes")
sce <- sce[markers,]
rse_gene <- rse_gene[markers,]

# Build signature from single-cell data
# This function builds a signature matrix using genes identified by the DEAnalysisMAST() function.
message(Sys.time(), "- convert counts to matrix")
sc_matrix <- as.matrix(counts(sce))

message(Sys.time(), "- buildSignatureMatrix")
Signature <- buildSignatureMatrixMAST(scdata=sc_matrix,
                                      id=sce$cellType_broad_hc,
                                      path=data_dir_marker,
                                      diff.cutoff=0.5,
                                      pval.cutoff=0.01)

head(Signature)
dim(Signature)

#trim signature and bulk data to contain the same differentially expressed genes
tr<-trimData(Signature,assays(rse_gene)$counts[,1])

dim(tr$sig)
length(tr$bulk)

#estimate using dampened weighted least squares
message(Sys.time(), "- DWLS")
solDWLS<-solveDampenedWLS(tr$sig,tr$bulk)

Signature <- Signature[rownames(Signature) %in% rownames(rse_gene),]
rse_gene <- rse_gene[rownames(Signature),]

est_prop_dwls <- apply(assays(rse_gene)$counts,2,solveDampenedWLS,S = Signature)
est_prop_dwls <- t(est_prop_dwls)

save(est_prop_dwls, Signature, file = here(data_dir, paste0("est_prop_dwls-", marker_label, ".Rdata")))

# slurmjobs::job_single('04_deconvolution_DWLS_FULL', create_shell = TRUE, memory = '100G', command = "Rscript 04_deconvolution_DWLS.R FULL")

# slurmjobs::job_single('04_deconvolution_DWLS_MeanRatio_top25', create_shell = TRUE, memory = '25G', command = "Rscript 04_deconvolution_DWLS.R MeanRatio_top25 ../../processed-data/08_bulk_deconvolution/markers_MeanRatio_top25.txt")
# slurmjobs::job_single('04_deconvolution_DWLS_MeanRatio_MAD3', create_shell = TRUE, memory = '10G', command = "Rscript 04_deconvolution_DWLS.R MeanRatio_MAD3 ../../processed-data/08_bulk_deconvolution/markers_MeanRatio_MAD3.txt")
# slurmjobs::job_single('04_deconvolution_DWLS_MeanRatio_over2', create_shell = TRUE, memory = '10G', command = "Rscript 04_deconvolution_DWLS.R MeanRatio_over2 ../../processed-data/08_bulk_deconvolution/markers_MeanRatio_over2.txt")

# slurmjobs::job_single('04_deconvolution_DWLS_1vALL_top25', create_shell = TRUE, memory = '25G', command = "Rscript 04_deconvolution_DWLS.R 1vALL_top25 ../../processed-data/08_bulk_deconvolution/markers_1vALL_top25.txt")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()




