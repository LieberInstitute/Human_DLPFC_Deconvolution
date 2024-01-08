
library("SingleCellExperiment")
library("here")
library("sessioninfo")

# Tab-delimited tabular input format (.txt) with no double quotations and no missing entries.
data_dir <- here("processed-data" , "08_bulk_deconvolution", "07_deconvolution_CIBERSORTx_prep")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load DLPFC data ####
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

rownames(rse_gene) <- rowData(rse_gene)$ensemblID

rse_counts <- assays(rse_gene)$counts |>
  as.data.frame() |>
  tibble::rownames_to_column("GeneSymbol")

rse_counts[1:5,1:3]
#   GeneSymbol 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291_Br2720_Mid_Cyto
# 1 ENSG00000227232                            42                            72
# 2 ENSG00000278267                             6                             3
# 3 ENSG00000268903                             2                             7
# 4 ENSG00000269981                             1                             1
# 5 ENSG00000279457                           131                           210

write.table(rse_counts, file = here(data_dir, "DLPFC_bulk_counts.txt"), sep = "\t")

## 
## sce data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]

# ## to test
# sce <- sce[,sce$cellType_broad_hc %in% c("Oligo", "Excit", "Inhib")]
# sce <- sce[sample(rownames(sce), 2000), sample(colnames(sce), 2000)]

sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)
table(sce$cellType_broad_hc)

rownames(sce) <- rowData(sce)$gene_id

message(Sys.time(), " - Format sce counts")
sce_counts <- assays(sce)$counts |>
  as.data.frame() |>
  tibble::rownames_to_column("GeneSymbol")

sce_counts[1:5,1:3]

message(Sys.time(), " - Export sce counts to ", data_dir)
write.table(sce_counts, file = here(data_dir, "DLPFC_sc_counts.txt"), sep = "\t")

# slurmjobs::job_single(name = "07_deconvolution_CIBERSORTx_prep", memory = "100G", cores = 1, create_shell = TRUE, command = "Rscript 07_deconvolution_CIBERSORTx_prep.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
