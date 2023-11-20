library("SummarizedExperiment")
library("tidyverse")
library("ComplexHeatmap")
library("sessioninfo")
library("here")

## prep dirs ##
plot_dir <- here("plots", "06_marker_genes", "01_example_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

set.seed(3042023)

#### prep example data ####
n_ct <- 4
cell_types <- factor(paste0("cell_type_", 1:n_ct)) ## add 0 pad if over 10

## prep colors
# c("#f2efea", "#fc7753", "#66d7d1", "#403d58", "#dbd56e","#1f271b")
example_colors <- c("#540d6e", "#ee4266", "#ffd23f", "#619B8A", "#1f271b")[1:n_ct]
names(example_colors) <- cell_types

#### Ideal Heatmap ####

# number of marker genes for each cell type
n_gene <- 4
n_total_gene <- n_gene * n_ct

marker_matrix <- matrix(rep(cell_types, each = n_total_gene), nrow = n_total_gene) ==
    matrix(rep(cell_types, each = n_total_gene), nrow = n_total_gene, byrow = TRUE)

mode(marker_matrix) <- "integer"

gene_names <- paste0("gene_", str_pad(1:n_total_gene, nchar(n_total_gene), "left", pad = "0"))

colnames(marker_matrix) <- cell_types
rownames(marker_matrix) <- gene_names

marker_matrix_var <- marker_matrix + abs(rnorm(n = n_total_gene * n_ct, sd = 0.1))

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
    rect_gp = gpar(col = "grey25", lwd = 1)
)
dev.off()


pdf(here(plot_dir, "ideal_heatmap_split.pdf"), height = 5.5, width = 6)
Heatmap(marker_matrix_var,
    name = "Z Score",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    top_annotation = column_ha,
    right_annotation = row_ha,
    column_split = cell_types,
    row_split = paste0("markers\n", rep(cell_types, each = n_gene)),
    # rect_gp = gpar(col = "grey50", lwd = 1),
    show_row_names = TRUE,
    show_column_names = FALSE
)
dev.off()

## (Neuron) Excit, Oligo, Astro Example

marker_matrix_ct_example <- marker_matrix_var[1:12,1:3]
example_ct <- c("Neuron","Oligo","Astro")
colnames(marker_matrix_ct_example) <- example_ct

example_colors2 <- c(Neuron = "#DB77BA", Oligo = "#E7872B", Astro = "#5AAA46")

column_ha_ct <- HeatmapAnnotation(
  cell_type = example_ct,
  col = list(cell_type = example_colors2),
  show_legend = FALSE
)

row_ha_ct <- rowAnnotation(
  cell_type = rep(example_ct, each = n_gene),
  col = list(cell_type = example_colors2),
  show_legend = FALSE
)

pdf(here(plot_dir, "heatmap_ct_example.pdf"), height = 5, width = 3)
# png(here(plot_dir, "heatmap_ct_example.png"), height = 500, width = 300)
Heatmap(marker_matrix_ct_example,
        name = "Expression",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        # top_annotation = column_ha_ct,
        left_annotation = row_ha_ct,
        # column_split = example_ct,
        # row_split = paste0("markers\n", rep(example_ct, each = n_gene)),
        # rect_gp = gpar(col = "grey50", lwd = 1),
        show_row_names = FALSE,
        show_column_names = TRUE
)
dev.off()

## est prop
set.seed(512)
example_prop <- matrix(data = rep(c(6,3), 4), nrow = 4, byrow = TRUE) - matrix(rnorm(8), nrow = 4)
example_prop[1,1] <- example_prop[1,1] -1
example_prop <- cbind(example_prop, 10 - rowSums(example_prop))/10
rowSums(example_prop)

colnames(example_prop) <- c("Neuron", "Astro","Oligo")
rownames(example_prop) <- paste0("sample_", 1:nrow(example_prop))
example_prop

#          Neuron     Astro      Oligo
# sample_1 0.5291505 0.4182825 0.05256700
# sample_2 0.6108602 0.2917233 0.09741658
# sample_3 0.5631502 0.3516663 0.08518346
# sample_4 0.4060879 0.3460245 0.24788759

example_prop_long <- reshape2::melt(example_prop) |>
  rename(Sample = Var1, cell_type = Var2, prop = value) |>
  mutate(cell_type = factor(cell_type, levels = c("Oligo", "Astro", "Neuron")))

example_prop_plot <- DeconvoBuddies::plot_composition_bar(example_prop_long, 
                                                          sample_col = "Sample",
                                                          x_col = "Sample",
                                                          add_text = FALSE) +
  scale_fill_manual(values = example_colors2) +
  theme_void()

ggsave(example_prop_plot, filename = here(plot_dir, "example_composition.png"), width = 5)

print(example_prop_plot)

#### Mean Ratio Cartoon ####

example_expression <- tibble(
    gene = "specific marker",
    cell_type_1 = abs(rnorm(n = 1000, mean = 5, sd = .5)),
    cell_type_2 = abs(rnorm(n = 1000)),
    cell_type_3 = abs(rnorm(n = 1000)),
    cell_type_4 = abs(rnorm(n = 1000))
) |>
    bind_rows(tibble(
        gene = "unspecific marker",
        cell_type_1 = abs(rnorm(n = 1000, mean = 5, sd = .5)),
        cell_type_2 = abs(rnorm(n = 1000, mean = 4, sd = .5)),
        cell_type_3 = abs(rnorm(n = 1000)),
        cell_type_4 = abs(rnorm(n = 1000))
    )) |>
    pivot_longer(!gene, names_to = "cell_type", values_to = "expression")


# g1_expression <- tibble(gene = "specific marker",
#                         cell_type_1 = floor(rnorm(n = 1000, mean = 10^5, sd = 1000)),
#                         cell_type_2 = floor(rnorm(n = 1000, sd = 1000))
# ) |>
#   pivot_longer(!gene, names_to = "cell_type", values_to = "counts") |>
#   mutate(counts = ifelse(counts < 0, 0, counts),
#          logcounts = log10(counts + 1))

example_violin <- example_expression |>
    ggplot(aes(x = cell_type, y = expression, fill = cell_type)) +
    geom_violin(scale = "width") +
    stat_summary(
        fun = mean,
        geom = "crossbar",
        width = 0.3
    ) +
    facet_wrap(~gene) +
    scale_fill_manual(values = example_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(example_violin, filename = here(plot_dir, "example_expression_violin.png"), height = 5)

#### Check fave markers for example distributions ####


# sgejobs::job_single('01_example_plots', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript 01_example_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
