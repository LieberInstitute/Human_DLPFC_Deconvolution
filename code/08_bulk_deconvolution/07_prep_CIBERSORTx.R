
library("SingleCellExperiment")
library("purrr")
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
  tibble::rownames_to_column("GeneSym")

rse_counts[1:5,1:3]
#           GeneSym 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291_Br2720_Mid_Cyto
# 1 ENSG00000227232                            42                            72
# 2 ENSG00000278267                             6                             3
# 3 ENSG00000268903                             2                             7
# 4 ENSG00000269981                             1                             1
# 5 ENSG00000279457                           131                           210

write.table(rse_counts, file = here(data_dir, "DLPFC_bulk_counts.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

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
# colnames(sce) <- colData(sce)$cellType_broad_hc ## colnames as cell type cat

## subset marker genes
marker_files <- list(MeanRatio_top25 = "markers_MeanRatio_top25.txt", 
                     MeanRatio_MAD3 = "markers_MeanRatio_MAD3.txt", 
                     MeanRatio_over2 = "markers_MeanRatio_over2.txt", 
                     `1vALL_top25` = "markers_1vALL_top25.txt")

marker_gene_sets <- map(marker_files, ~scan(here("processed-data", "08_bulk_deconvolution", .x), what="", sep="\n"))
marker_gene_sets <- c(marker_gene_sets, list(FULL = rownames(sce)))

## temp create 2 new inputs
marker_gene_sets <- marker_gene_sets[c("MeanRatio_MAD3", "MeanRatio_over2")]

map_int(marker_gene_sets, length)

walk2(marker_gene_sets, names(marker_gene_sets), function(set, name){
  message(Sys.time(), " - Format sce counts ", name)
  # First, subset to keep just the marker genes. Then keep cells where 5% of
  # the markers have nonzero expression (this prevents errors in SVD due to
  # insufficient column-wise variance)
  sce_sub = sce[set,]
  keep_indices = colSums(assays(sce_sub)$counts > 0) >= 0.05 * length(set)
  message(
    sprintf(
      "Retaining %s%% of cells after filtering low expression",
      round(100 * mean(keep_indices), 1)
    )
  )
  sce_sub = sce_sub[, keep_indices]
  message("Remaining cells of each type:")
  print(table(sce_sub$cellType_broad_hc))

  sce_counts <- assays(sce_sub)$counts |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    as.matrix()

  print(dim(sce_counts))
  colnames(sce_counts) <- c("gene", as.character(sce_sub$cellType_broad_hc))
  # sce_counts[1:5,1:5]

  message(Sys.time(), " - Export")
  write.table(sce_counts, file = here(data_dir, paste0("DLPFC_sc_counts-",name,".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
})


# slurmjobs::job_single(name = "07_prep_CIBERSORTx", memory = "100G", cores = 1, create_shell = TRUE, command = "Rscript 07_prep_CIBERSORTx.R")

## CIBERSORTx Runs
# slurmjobs::job_single(name = "07_deconvolution_CIBERSORTx_FULL", memory = "25G", cores = 1, create_shell = TRUE, command = "TBD")
# slurmjobs::job_single(name = "07_deconvolution_CIBERSORTx_MeanRatio_top25", memory = "25G", cores = 1, create_shell = TRUE, command = "TBD")
# slurmjobs::job_single(name = "07_deconvolution_CIBERSORTx_1vALL_top25", memory = "25G", cores = 1, create_shell = TRUE, command = "TBD")

## Reproducibility information
print("Reproducibility information:")
gc()
Sys.time()
proc.time()
options(width = 120)
session_info()
