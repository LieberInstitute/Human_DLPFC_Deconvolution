
# install.packages("hspe_0.1.tar.gz")

library("hspe")
library("SingleCellExperiment")
library("jaffelab")
library("tidyverse")
library("here")
library("sessioninfo")

## get args
args = commandArgs(trailingOnly=TRUE)
marker_label <- args[1]

## not using txt list of marker genes, methods needs named list
marker_sets <- c("MeanRatio_top25", "1vALL_top25", "MeanRatio_MAD3", "MeanRatio_over2")
stopifnot(marker_label %in% c(marker_sets, "FULL"))

#### data output folder ####
data_dir <- here("processed-data","08_bulk_deconvolution", "05_deconvolution_hspe")
if(!dir.exists(data_dir)) dir.create(data_dir)

#### load data ####
## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## pseudobulk sce data
sce_pb <- readRDS(here("processed-data", "sce","sce_broad_pseudobulk.rds"))
table(sce_pb$cellType_broad_hc)
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib 
# 18        17        18        18        19        19        19

## use ensemblIDs
rownames(sce_pb) <- rowData(sce_pb)$gene_id

## find common genes
common_genes <- intersect(rowData(sce_pb)$gene_id, rowData(rse_gene)$ensemblID)
message("common genes: ", length(common_genes))
# [1] 13311 

# we can instead explicitly pass a list of markers to hspe specifying the marker genes
# elements of the list correspond one to each cell type in the same order specified either in elements of pure_samples

pure_samples = rafalib::splitit(sce_pb$cellType_broad_hc)
# hspe assumes log2 transformed expressions
mixture_samples = t(assays(rse_gene)$logcounts[common_genes,])
reference_samples = t(assays(sce_pb)$logcounts[common_genes,])

stopifnot(ncol(mixture_samples) == ncol(reference_samples))

est_prop_hspe <- NULL

if(marker_label == "FULL"){
  
  message(Sys.time(), "- hspe")
  ## n_markers...If not specified then top 10% of genes are used. 1vALL ratio method
  est_prop_hspe = hspe(Y = mixture_samples,
                       reference = reference_samples,
                       pure_samples = pure_samples,
                       seed =10524)

} else if(marker_label %in% marker_sets){ ## Run with our markers 
  
  ## load marker gene data
  load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)
  marker_tab <- NULL
  
  if(marker_label == "MeanRatio_top25"){
    
    marker_tab <- marker_stats |> 
      dplyr::filter(gene %in% common_genes, 
                    rank_ratio <= 25) 
    
  } else if(marker_label == "MeanRatio_MAD3"){
    
    marker_tab <- marker_stats |>
      dplyr::filter(gene %in% common_genes, ratio > 1) |>
      dplyr::group_by(cellType.target) |>
      dplyr::filter(ratio > median(ratio) + 3*mad(ratio))
    
  }else if(marker_label == "MeanRatio_over2"){
    
    marker_tab <- marker_stats |> 
      dplyr::filter(gene %in% common_genes, 
                    ratio > 2) 
    
  }else if(marker_label == "1vALL_top25"){
    
    marker_tab <- marker_stats |> 
      dplyr::filter(gene %in% common_genes, 
                    rank_marker <= 25) 
    
  } else {
    stop("non-valid marker label: ", marker_label)
  }
  
  marker_genes <- purrr::map(rafalib::splitit(marker_tab$cellType.target), ~marker_tab$gene[.x])
    
  marker_genes <- marker_genes[names(pure_samples)]
  
  message(Sys.time(), "- hspe w/ markers ", marker_label)
  purrr::map_int(marker_genes, length)
  
  ## notes
  
  ## markers: Marker gene indices...list is a vector of indices (columns of ‘Y’) that will be considered markers of that particular type.
  
  est_prop_hspe = hspe(Y = mixture_samples,
                       reference = reference_samples,
                       pure_samples = pure_samples,
                       markers = marker_genes,
                       seed = 10524)
} else { ## not valid marker condition
  stop()
}

message(Sys.time(), "- Saving")
save(est_prop_hspe, file = here(data_dir,paste0("est_prop_hspe-",marker_label,".Rdata")))

# slurmjobs::job_single('05_deconvolution_hspe_FULL', create_shell = TRUE, memory = '25G', command = "Rscript 05_deconvolution_hspe.R FULL")

# slurmjobs::job_single('05_deconvolution_hspe_MeanRatio_top25', create_shell = TRUE, memory = '25G', command = "Rscript 05_deconvolution_hspe.R MeanRatio_top25")
# slurmjobs::job_single('05_deconvolution_hspe_MeanRatio_MAD3', create_shell = TRUE, memory = '5G', command = "Rscript 05_deconvolution_hspe.R MAD3")
# slurmjobs::job_single('05_deconvolution_hspe_MeanRatio_over2', create_shell = TRUE, memory = '5G', command = "Rscript 05_deconvolution_hspe.R over2")

# slurmjobs::job_single('05_deconvolution_hspe_1vALL_top25', create_shell = TRUE, memory = '25G', command = "Rscript 05_deconvolution_hspe.R 1vALL_top25")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

