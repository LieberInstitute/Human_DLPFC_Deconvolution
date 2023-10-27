
library("SingleCellExperiment")
library("ComplexHeatmap")
library("tidyverse")
library("here")
library("sessioninfo")
library("circlize")

## prep plot_dir
plot_dir <- here("plots", "06_marker_genes", "06_marker_gene_heatmap")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


## load colors
load(here("processed-data", "00_data_prep", "bulk_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)

## load pseudo-bulked data
sce_pb <- readRDS(here("processed-data", "sce","sce_broad_pseudobulk.rds"))
sce_pb
# dim: 15118 128

## marker_stats
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)

rd <- as.data.frame(rowData(rse_gene))

marker_stats <- marker_stats |>
  mutate(in_bulk = gene %in% rd$ensemblID)

# marker_logcounts <- logcounts(sce_pb)[marker_genes_top10$gene, ]
# marker_logcounts[1:5,1:5]
# rownames(marker_logcounts) <- rowData(sce_pb[marker_genes_top10$gene, ])$gene_name
# 
# marker_z_score <- scale(t(marker_logcounts))
# corner(marker_z_score)
# dim(marker_z_score)
# # [1] 128  70

get_marker_z_score <- function(n = 10){
  marker_genes_top <- marker_stats  |>
    group_by(cellType.target) |>
    arrange(-ratio) |>
    slice(1:n)
  
  marker_logcounts <- logcounts(sce_pb)[marker_genes_top$gene, ]
  rownames(marker_logcounts) <- rowData(sce_pb[marker_genes_top$gene, ])$gene_name
  marker_z_score <- scale(t(marker_logcounts))
  ## build annotations
  column_ha <- HeatmapAnnotation(
    cell_type = marker_genes_top$cellType.target,
    col = list(cell_type = cell_type_colors_broad),
    show_legend = FALSE
  )
  return(list(z_score = marker_z_score, marker_col_ha = column_ha))
}

maker_data <- map(c(top5 = 5, top10 = 10, top25 = 25), get_marker_z_score)

## logFC markers
get_FC_z_score <- function(n = 10){
  marker_genes_top <- marker_stats  |>
    group_by(cellType.target) |>
    arrange(-logFC) |>
    slice(1:n)
  
  marker_logcounts <- logcounts(sce_pb)[marker_genes_top$gene, ]
  rownames(marker_logcounts) <- rowData(sce_pb[marker_genes_top$gene, ])$gene_name
  marker_z_score <- scale(t(marker_logcounts))
  ## build annotations
  column_ha <- HeatmapAnnotation(
    cell_type = marker_genes_top$cellType.target,
    col = list(cell_type = cell_type_colors_broad),
    show_legend = FALSE
  )
  return(list(z_score = marker_z_score, marker_col_ha = column_ha))
}

fc_data <- map(c(top5 = 5, top10 = 10, top25 = 25), get_FC_z_score)


## Row annotations
BrNum_colors <- c("#ffb6b4",
                 "#c9a57d",
                 "#ffec9f",
                 "#c1cf83",
                 "#a4e6aa",
                 "#80b6a6",
                 "#61dfca",
                 "#8ffdff",
                 "#fce1ff",
                 "#c5a1bb")

names(BrNum_colors) <- unique(sce_pb$BrNum)

row_ha <- rowAnnotation(
  cell_type = sce_pb$cellType_broad_hc,
  BrNum = sce_pb$BrNum,
  ncells = anno_barplot(sce_pb$ncells),
  col = list(cell_type = cell_type_colors_broad,
             BrNum = BrNum_colors)
)

## just cell type
simple_row_ha <- rowAnnotation(
  cell_type = sce_pb$cellType_broad_hc,
  col = list(cell_type = cell_type_colors_broad)
)

## define color scale for top5

max(rbind(maker_data$top5$z_score, fc_data$top5$z_score))
min(rbind(maker_data$top5$z_score, fc_data$top5$z_score))

col_fun = colorRamp2(c(min(rbind(maker_data$top5$z_score, fc_data$top5$z_score)),
                       0,
                       max(rbind(maker_data$top5$z_score, fc_data$top5$z_score))), 
                     c("blue","white" ,"red"))


## need to fine-tune  params for different  n

## top 5
pdf(here(plot_dir, "markers_heatmap_top5.pdf"), width = 5, height = 5)
Heatmap(maker_data$top5$z_score,
        name = "z score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        right_annotation = simple_row_ha,
        top_annotation = maker_data$top5$marker_col_ha,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE,
        col = col_fun
)
dev.off()

pdf(here(plot_dir, "markers_heatmap_top5_uncluster.pdf"), width = 5, height = 5)
Heatmap(maker_data$top5$z_score,
        name = "z score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        right_annotation = simple_row_ha,
        top_annotation = maker_data$top5$marker_col_ha,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE,
        col = col_fun
)
dev.off()


## logFC version
pdf(here(plot_dir, "markers_heatmap_top5logFC_uncluster.pdf"), width = 5, height = 5)
Heatmap(fc_data$top5$z_score,
        name = "z score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        right_annotation = simple_row_ha,
        top_annotation = fc_data$top5$marker_col_ha,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE,
        col = col_fun
)
dev.off()

pdf(here(plot_dir, "markers_heatmap_top5logFC.pdf"), width = 5, height = 5)
Heatmap(fc_data$top5$z_score,
        name = "z score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        right_annotation = simple_row_ha,
        top_annotation = fc_data$top5$marker_col_ha,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE,
        col = col_fun
)
dev.off()


## top 10 
pdf(here(plot_dir, "markers_heatmap_top10.pdf"), width = 8.5, height = 10.3)
# png(here(plot_dir, "markers_heatmap_layer.png"), height = 1000, width = 100)
Heatmap(maker_data$top10$z_score,
        name = "z score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        right_annotation = row_ha,
        top_annotation = maker_data$top10$marker_col_ha,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE
        # ,
        # col = viridis::viridis(100)
)
dev.off()

pdf(here(plot_dir, "markers_heatmap_top10logFC.pdf"), width = 8.5, height = 10.3)
Heatmap(fc_data$top10$z_score,
        name = "z score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        right_annotation = row_ha,
        top_annotation = fc_data$top10$marker_col_ha,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE
)
dev.off()


## top 25
pdf(here(plot_dir, "markers_heatmap_top25.pdf"), width = 10, height = 11)
# png(here(plot_dir, "markers_heatmap_layer.png"), height = 1000, width = 100)
Heatmap(maker_data$top25$z_score,
        name = "z score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        right_annotation = row_ha,
        top_annotation = maker_data$top25$marker_col_ha,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 4),
        show_row_names = FALSE
        # ,
        # col = viridis::viridis(100)
)
dev.off()


