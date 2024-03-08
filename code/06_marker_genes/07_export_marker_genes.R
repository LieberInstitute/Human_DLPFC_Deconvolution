
library("SummarizedExperiment")
library("tidyverse")
library("here")
library("sessioninfo")
library("jaffelab")

#### Load data ####
## marker_stats
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)

rd <- as.data.frame(rowData(rse_gene))

marker_stats2 <- marker_stats |>
  group_by(cellType.target) |>
  mutate(in_bulk = gene %in% rd$ensemblID,
         MeanRatio_top25 = in_bulk & rank_ratio <= 25,
         MeanRatio_over2 = in_bulk & ratio > 2,
         MeanRatio_MAD3 = in_bulk & (ratio > median(ratio) + 3*mad(ratio)),
         `1vALL_top25` = in_bulk & rank_marker <= 25
         )

marker_stats2 |> 
  ungroup()|> 
  summarise(MeanRatio_top25 = sum(MeanRatio_top25),
            MeanRatio_over2 = sum(MeanRatio_over2),
            MeanRatio_MAD3 = sum(MeanRatio_MAD3),
            `1vALL_top25` = sum(`1vALL_top25`),
            )

marker_genes_top25 <- marker_stats2 |>
  group_by(cellType.target) |>
  # filter(in_bulk, rank_ratio <=25)
  filter(MeanRatio_top25) |>
  arrange(-ratio)

marker_genes_top25 |> count(cellType.target)
## top 25 that overlap w/ rse_genes
# cellType.target     n
# <fct>           <int>
# 1 Astro              24
# 2 EndoMural          24
# 3 Micro              17
# 4 Oligo              25
# 5 OPC                18
# 6 Excit              23
# 7 Inhib              20

## for MR top25, whats the lowest MR?
marker_genes_top25 |>
  arrange(ratio) |>
  slice(1) |>
  arrange(ratio) |>
  select(cellType.target, gene, Symbol,rank_ratio, ratio, logFC)

# cellType.target gene            Symbol     rank_ratio ratio logFC
# <fct>           <chr>           <chr>           <int> <dbl> <dbl>
#   1 Inhib           ENSG00000120049 KCNIP2             25  3.39 0.963
# 2 OPC             ENSG00000231871 IPO9-AS1           25  3.75 1.12 
# 3 Astro           ENSG00000143603 KCNN3              25  4.10 1.25 
# 4 Excit           ENSG00000119125 GDA                25  4.45 3.81 
# 5 EndoMural       ENSG00000176641 RNF152             25  5.99 1.99 
# 6 Oligo           ENSG00000286749 AC008571.2         25  7.08 3.93 
# 7 Micro           ENSG00000169313 P2RY12             25 14.9  2.21 

#### export ####

marker_genes_top25_simple <- marker_genes_top25 |>
  select(cellType.target, gene, Symbol)

marker_genes_ensembl <- map(splitit(marker_genes_top25$cellType.target), ~marker_genes_top25$gene[.x])

save(marker_genes_top25, marker_genes_top25_simple, marker_genes_ensembl, file = here("processed-data","06_marker_genes","marker_genes_top25.Rdata"))

#### Write csv ####

write_csv(marker_stats2, here("processed-data", "06_marker_genes", "marker_stats_broad.csv"))
