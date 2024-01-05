
# install.packages("hspe_0.1.tar.gz")

library("hspe")
library("SingleCellExperiment")
library("jaffelab")
library("here")
library("sessioninfo")

#### hspe example from https://github.com/gjhunt/hspe/blob/main/vign/basic-deconvolution.md ####
# 
# names(shen_orr_ex)
# # [1] "data"       "annotation" "name"
# 
# head(shen_orr_ex$annotation$mixture)
# 
# # RNA-summarized gene expression data
# # The values of the matrix are log2 RMA-summarized gene expressions.
# Y <- shen_orr_ex$data$log
# Y[1:4,1:4]
# #           X1367566_at X1367568_a_at X1367570_at X1367584_at
# # GSM495209    3.396192      7.685769    5.722330    6.628653
# # GSM495210    2.882626      7.759002    6.005583    6.771917
# # GSM495211    3.072980      7.598871    5.741630    6.564820
# # GSM495212    3.168440      7.209959    6.396841    7.040779
# 
# dim(shen_orr_ex$data$log)
# # [1]  42 600
# 
# ## seubset for example
# data = shen_orr_ex$data$log[,c(1:10,201:210,401:410)]
# mixture_proportions = shen_orr_ex$annotation$mixture
# 
# pure_samples = list(Liver=c(1,2,3),Brain=c(4,5,6),Lung=c(7,8,9))
# mixture_samples = data[-(1:9),]
# reference_samples = data[1:9,]
# 
# dim(mixture_samples)
# dim(reference_samples)
# 
# mixture_samples[1:5,1:5]
# reference_samples[1:5:1:5]
# 
# out = hspe(Y=mixture_samples, reference=reference_samples,pure_samples = pure_samples)
# 
# names(out)
# head(out$estimates)
# #           Liver      Brain      Lung
# # GSM495218 0.04487030 0.23356650 0.7215632
# # GSM495219 0.05038660 0.23772982 0.7118836
# # GSM495220 0.05172627 0.23519603 0.7130777
# # GSM495221 0.54595621 0.04624920 0.4077946
# # GSM495222 0.54069326 0.05576177 0.4035450
# # GSM495223 0.54231228 0.05122310 0.4064646

#### Run on our data ####

# output_dir <- here("processed-data","08_bulk_deconvolution", "05_deconvoltion_hspe")
# if(!file.exists(output_dir)) dir.create(output_dir)


#### load data ####
## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce data
# load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)
# sce <- sce[,sce$cellType_broad_hc != "Ambigious"]

## pseudobulk sce data
sce_pb <- readRDS(here("processed-data", "sce","sce_broad_pseudobulk.rds"))
table(sce_pb$cellType_broad_hc)
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib 
# 18        17        18        18        19        19        19

## use ensemblIDs
rownames(sce_pb) <- rowData(sce_pb)$gene_id

common_genes <- intersect(rownames(rse_gene), rownames(sce_pb))
length(common_genes)

# we can instead explicitly pass a list of markers to hspe specifying the marker genes
# elements of the list correspond one to each cell type in the same order specified either in elements of pure_samples

pure_samples = rafalib::splitit(sce_pb$cellType_broad_hc)
# hspe assumes log2 transformed expressions
mixture_samples = t(assays(rse_gene)$logcounts[common_genes,])
reference_samples = t(assays(sce_pb)$logcounts[common_genes,])

ncol(mixture_samples) == ncol(reference_samples)

message(Sys.time(), "- hspe")
est_prop_hspe = hspe(Y = mixture_samples, 
                     reference = reference_samples,
                     pure_samples = pure_samples,
                     seed =10524)

message(Sys.time(), "- Saving")
save(est_prop_hspe, file = here("processed-data","08_bulk_deconvolution","est_prop_hspe.Rdata"))

#### Run with our markers ####

## load marker gene data
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)
# marker_stats

marker_stats |>
  dplyr::filter(gene %in% common_genes,
                rank_ratio <= 25) |>
  dplyr::count(cellType.target)

marker_tab <- marker_stats |> 
  dplyr::filter(gene %in% common_genes, rank_ratio <= 25) 

marker_genes <- purrr::map(rafalib::splitit(marker_tab$cellType.target), ~marker_tab$gene[.x])

marker_genes <- marker_genes[names(pure_samples)]

message(Sys.time(), "- hspe w/ markers")
est_prop_hspe = hspe(Y = mixture_samples, 
                     reference = reference_samples,
                     pure_samples = pure_samples,
                     markers = marker_genes,
                     seed = 10524)

message(Sys.time(), "- Saving")
save(est_prop_hspe, file = here("processed-data","08_bulk_deconvolution","est_prop_hspe_markers.Rdata"))


# slurmjobs::job_single(name = "05_deconvolution_hspe", memory = "25G", cores = 1, create_shell = TRUE, command = "Rscript 05_deconvolution_hspe.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

