
library("SingleCellExperiment")
library("purrr")
library("here")
library("sessioninfo")

# Tab-delimited tabular input format (.txt) with no double quotations and no missing entries.
data_dir <- here("processed-data" , "08_bulk_deconvolution", "07_deconvolution_CIBERSORTx_prep")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)


## sce data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]

sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)
table(sce$cellType_broad_hc)
rownames(sce) <- rowData(sce)$gene_id

## subset marker genes
marker_files <- sprintf("HVG%d0", 1:10)
marker_gene_sets <- map(marker_files, ~scan(here("processed-data", "06_marker_genes","09_HVGs", paste0(.x,".txt")), what="", sep="\n"))

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


# slurmjobs::job_single(name = "07_prep_CIBERSORTx_HVG", memory = "25G", cores = 1, create_shell = TRUE, command = "Rscript 07_prep_CIBERSORTx_HVG.R")

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
