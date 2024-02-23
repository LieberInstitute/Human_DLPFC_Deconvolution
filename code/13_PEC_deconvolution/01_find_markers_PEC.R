library("DeconvoBuddies")
library("SingleCellExperiment")
library("tidyverse")
library("ggrepel")
library("here")
library("sessioninfo")

#### prep dirs ####
data_dir <- here("processed-data", "13_PEC_deconvolution", "01_find_markers_PEC")
if(!dir.exists(data_dir)) dir.create(data_dir)

plot_dir <- here("processed-data", "13_PEC_deconvolution", "01_find_markers_PEC")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

## load PEC PTSD Brainomics DLPFC data
load(here("processed-data", "13_PEC_deconvolution", "sce_PTSDBrainomics.Rdata"), verbose = TRUE)

## drop cell types 
sce <- sce[, sce$cellType_broad != "drop"]
# sce$cellType_broad <- droplevels(sce$cellType_broad)

table(sce$cellType_broad)
# Astro EndoMural     Excit     Inhib     Micro     Oligo       OPC 
# 17605       985     86129     36642      9419     35866     10498 
message("any NA: ", any(is.na(sce$cellType_broad)))

## load bulk data find common genes
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
common_genes <- intersect(rowData(sce)$featureid, rowData(rse_gene)$ensemblID)
message("Common genes with bulk data: ", length(common_genes))
# [1] 16734
rm(rse_gene)

#### Run Marker Finding ####
rowData(sce)$Symbol <- rowData(sce)$gene_name

message(Sys.time(), " Find Mean Ratio Genes")
mean_ratio <- get_mean_ratio2(sce, cellType_col = "cellType_broad", assay_name = "X",add_symbol = TRUE)

message(Sys.time(), " Find 1vAll genes")
markers_1vALL <- findMarkers_1vAll(sce, cellType_col = "cellType_broad", mod = "~individualID", assay_name = "counts")

marker_stats <- mean_ratio |> 
  left_join(markers_1vALL) |>
  mutate(in_bulk = gene %in% common_genes,
         MeanRatio_top25 = (rank_ratio <= 25 & in_bulk))

message(Sys.time(), " Save marker gene data")
save(marker_stats, file = here(data_dir, "PEC_marker_stats_broad.Rdata"))
write_csv(marker_stats, file = here(data_dir, "PEC_marker_stats_broad.csv"))

## pull marker sets 
marker_stats |> 
  dplyr::filter(MeanRatio_top25) |>
  count(cellType.target)

markers_mean_ratio_top25 <- marker_stats |> 
  dplyr::filter(gene %in% common_genes, rank_ratio <= 25) |>
  dplyr::pull(gene)

cat(markers_mean_ratio_top25, sep = "\n", file = here("processed-data","13_PEC_deconvolution", "PEC_markers_MeanRatio_top25.txt"))

## hockey stick plots 
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)

mean_ratio_v_logFC <- marker_stats |>
  mutate(top25 = rank_ratio <= 25) |>
  ggplot(aes(x = ratio, y = std.logFC, color = top25)) +
  geom_point(alpha = 0.5) +
  ggrepel::geom_text_repel(aes(label = ifelse(rank_ratio ==1 | rank_marker ==1 , Symbol, "")), 
                  size = 2.5, 
                  color = "black") +
  facet_wrap(~cellType.target, scales = "free_x", nrow = 1) +
  # facet_wrap(~cellType.target, nrow = 1) +
  labs(x = "Mean Ratio", y = "1vALL Standard logFC") +
  theme_bw()
# theme(legend.position = "bottom")

ggsave(mean_ratio_v_logFC, filename = here(plot_dir, "PEC_mean_ratio_v_logFC_free.png"), width = 12, height = 2.5)
ggsave(mean_ratio_v_logFC, filename = here(plot_dir, "PEC_mean_ratio_v_logFC_free.pdf"), width = 12, height = 2.5)

# slurmjobs::job_single('01_find_markers_PEC', create_shell = TRUE, memory = '100G', command = "Rscript 01_find_markers_PEC.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
