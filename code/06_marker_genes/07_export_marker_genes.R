
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

marker_stats <- marker_stats |>
  mutate(in_bulk = gene %in% rd$ensemblID)

marker_stats |> count(cellType.target, in_bulk)


marker_genes_top25 <- marker_stats |>
  group_by(cellType.target) |>
  # filter(in_bulk, rank_ratio <=25)
  filter(in_bulk) |>
  arrange(-ratio) |>
  slice(1:25)

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

marker_genes_top25 |>
  arrange(ratio) |>
  slice(1) |>
  arrange(ratio) |>
  select(cellType.target, gene, Symbol,rank_ratio, ratio, logFC)

# cellType.target gene            Symbol     rank_ratio ratio logFC
# <fct>           <chr>           <chr>           <int> <dbl> <dbl>
# 1 OPC             ENSG00000103489 XYLT1              33  2.80  9.97
# 2 Inhib           ENSG00000132170 PPARG              30  3.00  1.79
# 3 Astro           ENSG00000286288 AL109809.5         26  3.90  1.59
# 4 Excit           ENSG00000106089 STX1A              27  4.36  1.97
# 5 EndoMural       ENSG00000137962 ARHGAP29           26  5.90  3.71
# 6 Oligo           ENSG00000286749 AC008571.2         25  7.08  3.93
# 7 Micro           ENSG00000109794 FAM149A            35  9.36  1.07

#### export ####

marker_genes_top25_simple <- marker_genes_top25 |>
  select(cellType.target, gene, Symbol)

marker_genes_ensembl <- map(splitit(marker_genes_top25$cellType.target), ~marker_genes_top25$gene[.x])

save(marker_genes_top25, marker_genes_top25_simple, marker_genes_ensembl, file = here("processed-data","06_marker_genes","marker_genes_top25.Rdata"))

