
library("SingleCellExperiment")
library("dplyr")
library("here")
library("sessioninfo")

#### load data ####
## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

rownames(sce) <- rowData(sce)$gene_id

## load marker gene data
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)
# marker_stats
dim(marker_stats)

## find common genes
common_genes <- intersect(rowData(sce)$gene_id, rowData(rse_gene)$ensemblID)
length(common_genes)
# [1] 17804


#### Mean Ratio ####
## top 25
markers_mean_ratio_top25 <- marker_stats |> 
  filter(rank_ratio <= 25) |>
  mutate(in_bulk = gene %in% rowData(rse_gene)$ensemblID)

markers_mean_ratio_top25 |>
  filter(in_bulk) |>
  dplyr::count(cellType.target)

# cellType.target     n
# <fct>           <int>
# 1 Astro              24
# 2 EndoMural          24
# 3 Micro              17
# 4 Oligo              25
# 5 OPC                18
# 6 Excit              23
# 7 Inhib              20

write.csv(markers_mean_ratio_top25, file = here("processed-data","08_bulk_deconvolution", "markers_mean_ratio_top25.csv"), row.names = FALSE)

markers_mean_ratio_top25 <- marker_stats |> 
  dplyr::filter(gene %in% common_genes, rank_ratio <= 25) |>
  dplyr::pull(gene)

cat(markers_mean_ratio_top25, sep = "\n", file = here("processed-data","08_bulk_deconvolution", "markers_MeanRatio_top25.txt"))
save(markers_mean_ratio_top25, file = here("processed-data","08_bulk_deconvolution", "markers_mean_ratio_top25.Rdata"))

## 3 MADs

markers_mean_ratio_MAD3 <- marker_stats |>
  filter(gene %in% common_genes, ratio > 1) |>
  group_by(cellType.target) |>
  filter(ratio > median(ratio) + 3*mad(ratio)) |>
  dplyr::pull(gene)

# cellType.target     n
# <fct>           <int>
# 1 Astro              35
# 2 EndoMural          49
# 3 Micro              39
# 4 Oligo              50
# 5 OPC                67
# 6 Excit             157
# 7 Inhib             123
length(markers_mean_ratio_MAD3) #[1] 520
cat(markers_mean_ratio_MAD3, sep = "\n", file = here("processed-data","08_bulk_deconvolution", "markers_MeanRatio_MAD3.txt"))


## MeanRatio > 2
markers_mean_ratio_2 <- marker_stats |>
  filter(gene %in% common_genes, ratio > 2) |>
  # count(cellType.target) |>
  dplyr::pull(gene)

# cellType.target     n
# <fct>           <int>
# 1 Astro              58
# 2 EndoMural          74
# 3 Micro              87
# 4 Oligo              81
# 5 OPC                68
# 6 Excit             137
# 7 Inhib              52

length(markers_mean_ratio_2) #[1] 557
cat(markers_mean_ratio_2, sep = "\n", file = here("processed-data","08_bulk_deconvolution", "markers_MeanRatio_over2.txt"))


## 1vALL top25
markers_1vALL_top25 <- marker_stats |> 
  dplyr::filter(gene %in% common_genes, rank_marker <= 25) |>
  dplyr::pull(gene)

cat(markers_1vALL_top25, sep = "\n", file = here("processed-data","08_bulk_deconvolution", "markers_1vALL_top25.txt"))


# markers_mean_ratio_top25 <- scan(here("processed-data","08_bulk_deconvolution", "markers_MeanRatio_top25.txt"), what="", sep="\n")
# markers_1vALL_top25 <- scan(here("processed-data","08_bulk_deconvolution", "markers_1vALL_top25.txt"), what="", sep="\n")

marker_stats |> 
  mutate(marker_1vALL = gene %in% markers_1vALL_top25,
         marker_MeanRatio = gene %in% markers_mean_ratio_top25) |>
  filter(marker_1vALL & marker_MeanRatio) |>
  count(cellType.target, marker_1vALL, marker_MeanRatio)

## overlap
# # A tibble: 7 × 4
# cellType.target marker_1vALL marker_MeanRatio     n
# <fct>           <lgl>        <lgl>            <int>
#   1 Astro           TRUE         TRUE                 9
# 2 EndoMural       TRUE         TRUE                11
# 3 Micro           TRUE         TRUE                11
# 4 Oligo           TRUE         TRUE                12
# 5 OPC             TRUE         TRUE                12
# 6 Excit           TRUE         TRUE                 3
# 7 Inhib           TRUE         TRUE                 8

  
# slurmjobs::job_single('00_pull_markerse', create_shell = TRUE, memory = '25G', command = "Rscript 00_pull_markers.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


