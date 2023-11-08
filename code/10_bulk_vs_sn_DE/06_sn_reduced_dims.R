library("SummarizedExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")
library("scater")

plot_dir <- here("plots", "10_bulk_vs_sn_DE", "06_sn_reduced_dims")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)

my_theme <- theme_bw() +
  theme(
    text = element_text(size = 15)
  )

## load our DLPFC data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

## drop Ambiguous nuc
sce <- sce[, sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)
table(sce$cellType_broad_hc)

TSNE_cellTypes_broad <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = cellType_broad_hc)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_color_manual(values = cell_type_colors_broad, name = "Cell\nType") +
  coord_equal() +
  labs(x = "tSNE Dimension 1", y = "tSNE Dimension 2") +
  my_theme +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 3)) +
  theme(legend.position = "bottom", legend.text=element_text(size=10))

ggsave(TSNE_cellTypes_broad, filename = here(plot_dir, "TSNE_cellTypes_broad.png"), width = 4, height = 5)

