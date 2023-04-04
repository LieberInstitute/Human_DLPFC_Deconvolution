
library("SummarizedExperiment")
library("tidyverse")
library("ComplexHeatmap")
library("sessioninfo")
library("here")

## prep dirs ##
plot_dir <- here("plots", "06_marker_genes", "01_example_plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

set.seed(3042023)

#### prep example data ####
n_ct <- 4
cell_types <- factor(paste0("cell_type_", 1:n_ct)) ## add 0 pad if over 10

## prep colors
# c("#f2efea", "#fc7753", "#66d7d1", "#403d58", "#dbd56e","#1f271b")
example_colors <- c("#540d6e","#ee4266","#ffd23f","#619B8A", "#1f271b")[1:n_ct]
names(example_colors) <- cell_types

#### Ideal Heatmap ####

# number of marker genes for each cell type
n_gene <- 4
n_total_gene <- n_gene * n_ct
  
marker_matrix <- matrix(rep(cell_types, each = n_total_gene), nrow =  n_total_gene) ==
  matrix(rep(cell_types, each = n_total_gene), nrow =  n_total_gene, byrow = TRUE)

mode(marker_matrix) <- "integer"

gene_names <- paste0("gene_", str_pad(1:n_total_gene, nchar(n_total_gene), "left", pad = "0"))

colnames(marker_matrix) <- cell_types
rownames(marker_matrix) <- gene_names

marker_matrix_var <- marker_matrix + abs(rnorm(n = n_total_gene*n_ct, sd = 0.1))

## annotations
column_ha <- HeatmapAnnotation(
  cell_type = cell_types,
  col = list(cell_type = example_colors),
  show_legend = FALSE
)

row_ha <- rowAnnotation(
  cell_type = rep(cell_types, each = n_gene),
  col = list(cell_type = example_colors),
  show_legend = FALSE
)

## plot heatmap
pdf(here(plot_dir, "ideal_heatmap.pdf"))
Heatmap(marker_matrix,
        name = "Expression",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        top_annotation = column_ha,
        right_annotation = row_ha, 
        rect_gp = gpar(col = "grey25", lwd = 1))
dev.off()


pdf(here(plot_dir, "ideal_heatmap_split.pdf"))
Heatmap(marker_matrix_var,
        name = "Expression",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        top_annotation = column_ha,
        left_annotation = row_ha,
        column_split = cell_types,
        row_split = paste0("markers\n",rep(cell_types, each = n_gene)),
        rect_gp = gpar(col = "grey50", lwd = 1),
        show_row_names = TRUE,
        show_column_names = FALSE
        )
dev.off()

#### Mean Ratio Cartoon ####


# sgejobs::job_single('01_example_plots', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript 01_example_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
