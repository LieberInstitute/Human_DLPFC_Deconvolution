
library("splatter")
library("scater")
library("DeconvoBuddies")
library("SingleCellExperiment")
library("tidyverse")
library("sessioninfo")
library("here")


## prep dirs ##
plot_dir <- here("plots",  "06_marker_genes", "02_splat_example")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


## load our DLPFC data
# load(here("processed-data","sce","sce_DLPFC.Rdata"), verbose = TRUE)

# Estimate parameters from DLPFC data

# Error in UseMethod("splatEstimate") : 
#   no applicable method for 'splatEstimate' applied to an object of class "c('DelayedMatrix', 'DelayedArray', 'DelayedUnaryIsoOp', 'DelayedUnaryOp', 'DelayedOp', 'Array', 'RectangularData')"
# params <- splatEstimate(sce)
# 

## Try new params
set.seed(40423)
params <- newSplatParams()

sim <- splatSimulate(params,
                     nGenes = 1000,
                     batchCells = 1000, 
                     group.prob = c(0.3, 0.3, 0.2, 0.2),
                     method = "groups")

sim$cell_type <- factor(gsub("Group","cell_type_", sim$Group), levels = paste0("cell_type_", 1:4))
table(sim$cell_type)
# cell_type_1 cell_type_2 cell_type_3 cell_type_4 
# 301         314         184         201 

sim_colors <- c(cell_type_1 = "#540d6e",
                cell_type_2 = "#ee4266",
                cell_type_3 = "#ffd23f",
                cell_type_4 = "#619B8A")

# Use scater to calculate logcounts
sim <- logNormCounts(sim)
# Plot PCA
sim <- runPCA(sim)

pdf(here(plot_dir, "splat_PCA.pdf"))
plotPCA(sim, colour_by = "Group")
dev.off()

# Access the counts
counts(sim)[1:5, 1:5]

## Vis tSNE
sim <- runTSNE(sim, dimred="PCA")

pdf(here(plot_dir, "splat_TSNE.pdf"))
plotReducedDim(sim, dimred="TSNE", colour_by = "cell_type") +
  ggplot2::scale_color_manual(values = sim_colors)
dev.off()

#### Run Marker Finding ####

marker_mean_ratio <- get_mean_ratio2(sim, cellType_col = "cell_type", add_symbol = FALSE)
# marker_1vAll <- findMarkers_1vAll(sim, cellType_col = "cell_type", mod = NULL, add_symbol = FALSE)
## need to make findMarkers_1vAll work with no mod

## temp check 
fm <- scran::findMarkers(sim,
                         groups = sim$cell_type == "cell_type_1",
                         assay.type = "logcounts", test.type = "t",
                         std.lfc = TRUE,
                         direction = "up", pval.type = "all", full.stats = T
)

marker_1vAll <- fm[[2]]$stats.FALSE %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  tibble::add_column(cellType.target = "cell_type_1")

marker_stats <- marker_mean_ratio |> 
  left_join(marker_1vAll) |>
  mutate(Symbol = gene) ## TODO fix in plot_marker_express

hockey_plot <- marker_stats |>
  filter(cellType.target == "cell_type_1") |>
  ggplot(aes(x = ratio, y = logFC)) +
  geom_point()
  
ggsave(hockey_plot, filename = here(plot_dir, "hockey_stick_plot.png"))

plot_marker_express_ALL(sim, marker_stats, 
                        pdf_fn = here(plot_dir, "mean_ratio_markers.pdf"),
                        cellType_col = "cell_type", 
                        color_pal = sim_colors)


expres_plot <- plot_gene_express(sim, genes = c("Gene951","Gene190"), cat = "cell_type")
ggsave(expres_plot, filename = here(plot_dir, "markers_1vALL_ct1.png"))
