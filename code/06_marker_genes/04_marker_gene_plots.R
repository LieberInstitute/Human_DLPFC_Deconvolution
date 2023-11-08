library("DeconvoBuddies")
library("SingleCellExperiment")
library("tidyverse")
library("ggrepel")
library("here")
library("sessioninfo")

## prep plot_dir
plot_dir <- here("plots", "06_marker_genes", "04_marker_gene_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load sn Data + marker stats ####
## snRNA-seq data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

## drop Ambiguous nuc
sce <- sce[, sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

table(sce$cellType_broad_hc)
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib
# 3979      2157      1601     10894      1940     24809     11067

load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo
# cell_type_colors_broad

# cell_type_colors <- metadata(sce)$cell_type_colors_broad[levels(sce$cellType_broad_hc)]

rowData(sce)$Symbol <- rowData(sce)$gene_name
rownames(sce) <- rowData(sce)$gene_id

## marker_stats
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)

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
    geom_text_repel(aes(label = ifelse(top25, Symbol, "")), size = 1.5, color = "black") +
    facet_wrap(~cellType.target, scales = "free_x") +
    # facet_wrap(~cellType.target, nrow = 1) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_vline(xintercept = 1, linetype = "dashed") +
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
    n_genes = 10, # maybe plot 25
    rank_col = "rank_ratio",
    anno_col = "anno_ratio",
    cellType_col = "cellType_broad_hc",
    color_pal = cell_type_colors_broad,
    plot_points = FALSE
)

plot_marker_express_ALL(
  sce,
  stats = marker_stats,
  pdf_fn = here(plot_dir, "marker_expression_mean_ratio_broad_top25.pdf"),
  n_genes = 25, # maybe plot 25
  rank_col = "rank_ratio",
  anno_col = "anno_ratio",
  cellType_col = "cellType_broad_hc",
  color_pal = cell_type_colors_broad,
  plot_points = FALSE
)

marker_stats |>
    arrange(rank_marker) |>
    select(gene, cellType.target, ratio, rank_ratio, std.logFC, rank_marker, anno_logFC)

plot_marker_express_ALL(
    sce,
    stats = marker_stats,
    pdf_fn = here(plot_dir, "marker_expression_1vALL_broad.pdf"),
    n_genes = 10, # maybe plot 25
    rank_col = "rank_marker",
    anno_col = "anno_logFC",
    cellType_col = "cellType_broad_hc",
    color_pal = cell_type_colors_broad,
    plot_points = FALSE
)

## top 10 Inhib genes ##

top5_Inhib_mr <- plot_marker_express(sce, 
                    stats = marker_stats, 
                    cell_type = "Inhib",
                    n_genes = 5,
                    cellType_col = "cellType_broad_hc",
                    color_pal = cell_type_colors_broad,
                    ncol = 5)

ggsave(top5_Inhib_mr, filename = here(plot_dir, "top5_Inhib_MR.png"), height = 3)

top5_Inhib_1vAll <- plot_marker_express(sce, 
                                     stats = marker_stats, 
                                     cell_type = "Inhib",
                                     n_genes = 5,
                                     rank_col = "rank_marker",
                                     anno_col = "anno_logFC",
                                     cellType_col = "cellType_broad_hc",
                                     color_pal = cell_type_colors_broad,
                                     ncol = 5)

ggsave(top5_Inhib_1vAll, filename = here(plot_dir, "top5_Inhib_1vAll.png"), height = 3)

library(patchwork)

top5_Inhib_1vAll <- top5_Inhib_1vAll + labs(title = "Top 5: 1vALL logFC") + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
top5_Inhib_mr <- top5_Inhib_mr + labs(title = "Top 5: Mean Ratio") 

top5_inhib <- top5_Inhib_1vAll/top5_Inhib_mr
ggsave(top5_inhib, filename = here(plot_dir, "top5_Inhib.png"), height = 4.5, width = 9)
ggsave(top5_inhib, filename = here(plot_dir, "top5_Inhib.pdf"), height = 4.5, width = 9)

## MBP example ##

marker_stats |> filter(Symbol %in% c("MBP", "FOLH1"), cellType.target == "Oligo") |> select(-gene, - cellType.target)
# mean.target cellType   mean ratio rank_ratio Symbol anno_ratio         logFC log.p.value log.FDR std.logFC rank_marker
# <dbl> <fct>     <dbl> <dbl>      <int> <chr>  <chr>              <dbl>       <dbl>   <dbl>     <dbl>       <int>
# 1        1.60 Micro    0.0741 21.6           1 FOLH1  Oligo/Micro: 21.6…  1.37      -8418.  -8410.      1.53          23
# 2        5.45 Micro    2.03    2.68         55 MBP    Oligo/Micro: 2.681 12.6      -11093. -11084.      1.80           3
rownames(sce) <- rowData(sce)$gene_name

sce$Oligo <- ifelse(sce$cellType_broad_hc == "Oligo", "Oligo", "Other")

oligo_example1 <- plot_gene_express(sce, 
                                    genes = c("MBP", "FOLH1"), 
                                    assay_name = "logcounts", 
                                    cat = "Oligo", 
                                    color_pal = cell_type_colors_broad)

ggsave(oligo_example1, filename = here(plot_dir, "Oligo_example1.png"), height = 3, width = 5)

oligo_example2 <- plot_gene_express(sce, 
                  genes = c("MBP", "FOLH1"), 
                  assay_name = "logcounts", 
                  cat = "cellType_broad_hc", 
                  color_pal = cell_type_colors_broad)

ggsave(oligo_example2, filename = here(plot_dir, "Oligo_example2.png"), height = 3, width = 5)


ggsave(oligo_example1/oligo_example2, filename = here(plot_dir, "Oligo_example.png"), height = 6, width = 5)

