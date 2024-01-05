
# install.packages("DWLS") ## Use cran version

library("DWLS")
library("SingleCellExperiment")
library("here")
library("sessioninfo")

#### DWLS example ####

head(DWLS::dataSC_3)

#### Run on our data ####

output_dir <- here("processed-data","08_bulk_deconvolution", "04_deconvoltion_DWLS")
if(!file.exists(output_dir)) dir.create(output_dir)


#### load data ####
## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)
sce <- sce[,sce$cellType_broad_hc != "Ambigious"]

## use ensemblIDs
rownames(sce) <- rowData(sce)$gene_id

# ## subset for test
# sce <- sce[,sce$cellType_broad_hc %in% c("Oligo", "Excit", "Inhib")]
# sce <- sce[sample(rownames(sce), 2000), sample(colnames(sce), 2000)]


sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)
table(sce$cellType_broad_hc)
# Oligo Excit Inhib 
# 237   521   242 

# Build signature from single-cell data
# This function builds a signature matrix using genes identified by the DEAnalysisMAST() function.
message(Sys.time(), "- convert counts to matrix")
sc_matrix <- as.matrix(counts(sce))

message(Sys.time(), "- buildSignatureMatrix")
Signature <- buildSignatureMatrixMAST(scdata=sc_matrix,
                                      id=sce$cellType_broad_hc,
                                      path=output_dir,
                                      diff.cutoff=0.5,
                                      pval.cutoff=0.01)

head(Signature)
dim(Signature)
## use ensemblIDs
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
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

save(est_prop_dwls, file = here("processed-data","08_bulk_deconvolution","est_prop_dwls.Rdata"))

# slurmjobs::job_single(name = "04_deconvolution_DWLS", memory = "100G", cores = 1, create_shell = TRUE, command = "Rscript 04_deconvolution_DWLS.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()




