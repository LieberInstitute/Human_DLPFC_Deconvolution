
library("DeconvoBuddies")
library("SingleCellExperiment")
library("tidyverse")
library("ggrepel")
library("here")
library("sessioninfo")

## prep plot_dir
plot_dir <- here("plots", "06_marker_genes", "04_marker_gene_plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load sn Data + marker stats ####
## snRNA-seq data
load(here("processed-data","sce","sce_DLPFC.Rdata"), verbose = TRUE)

## drop Ambiguous nuc
sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

table(sce$cellType_broad_hc)
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib 
# 3979      2157      1601     10894      1940     24809     11067 

cell_type_colors <- metadata(sce)$cell_type_colors_broad[levels(sce$cellType_broad_hc)]

rowData(sce)$Symbol <- rowData(sce)$gene_name
rownames(sce) <- rowData(sce)$gene_id

## marker_stats
load(here("processed-data","06_marker_genes","03_find_markers_broad","marker_stats_broad.Rdata"), verbose = TRUE)

#### Mean ratio vs. log FC (Hockey Stick) Plots ####

mean_ratio_v_logFC <- marker_stats |>
  mutate(top25 = rank_ratio <= 25) |>
  ggplot(aes(x = ratio, y = std.logFC, color = top25)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~cellType.target, scales = "free_x", nrow = 1) +
  # facet_wrap(~cellType.target, nrow = 1) +
  labs(x = "Mean Ratio", y = "1vALL Standard logFC") +
  theme_bw() 
  # theme(legend.position = "bottom")

ggsave(mean_ratio_v_logFC, filename = here(plot_dir, "mean_ratio_v_logFC_free.png"), width = 12, height = 4)
ggsave(mean_ratio_v_logFC, filename = here(plot_dir, "mean_ratio_v_logFC_free.pdf"), width = 12, height = 4)

mean_ratio_v_logFC_detail <- marker_stats |>
  mutate(top25 = rank_ratio <= 25) |>
  ggplot(aes(x = ratio, y = std.logFC, color = top25)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(aes(label = ifelse(top25, Symbol,"")), size = 1.5, color = "black") +
  facet_wrap(~cellType.target, scales = "free_x") +
  # facet_wrap(~cellType.target, nrow = 1) +
  geom_hline(yintercept = 1, linetype = "dashed")+
  geom_vline(xintercept = 1, linetype = "dashed")+
  labs(x = "Mean Ratio", y = "1vALL Standard logFC") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(mean_ratio_v_logFC_detail, filename = here(plot_dir, "mean_ratio_v_logFC_detail.png"), height = 10, width = 8)
ggsave(mean_ratio_v_logFC_detail, filename = here(plot_dir, "mean_ratio_v_logFC_detail.pdf"), height = 10, width = 8)


#### plot expression of top Mean Ratio Markers ####
plot_marker_express_ALL(
  sce,
  stats = marker_stats,
  pdf_fn = here(plot_dir, "marker_expression_mean_ratio_broad.pdf"),
  n_genes = 10, #maybe plot 25
  rank_col = "rank_ratio",
  anno_col = "anno_ratio",
  cellType_col = "cellType_broad_hc",
  color_pal = cell_type_colors,
  plot_points = FALSE
)

marker_stats |> arrange(rank_marker) |> select(gene, cellType.target, ratio, rank_ratio, std.logFC, rank_marker, anno_logFC)

plot_marker_express_ALL(
  sce,
  stats = marker_stats,
  pdf_fn = here(plot_dir, "marker_expression_1vALL_broad.pdf"),
  n_genes = 10, #maybe plot 25
  rank_col = "rank_marker",
  anno_col = "anno_logFC",
  cellType_col = "cellType_broad_hc",
  color_pal = cell_type_colors,
  plot_points = FALSE
)

## top 10 Inhib genes ##

plot_marker_express(sce, ma)
