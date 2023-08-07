library("tidyverse")
# library("scales")
# library("ggrepel")
# library("patchwork")
library("here")
library("sessioninfo")
# library("DeconvoBuddies")

#### Dir Set-up ####
plot_dir <- here("plots", "03_HALO", "02_nuc_size_QC")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "03_HALO", "02_nuc_size_QC")
if (!dir.exists(data_dir)) dir.create(data_dir)

#### Plot Set-up ####
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)

#### Load Data ####
load(here("processed-data", "03_HALO", "halo_all.Rdata"), verbose = TRUE)

#### Density plots for Nucleus_Area by cell type ####

nuc_area_density <- halo_all |>
  ggplot(aes(x= Nucleus_Area, color = cell_type)) +
  geom_density() +
  theme_bw() +
  scale_color_manual(values = cell_type_colors_halo)

ggsave(nuc_area_density, filename = here(plot_dir, "Nucleus_Area_density.png"), width = 10)


nuc_area_density_sample <- halo_all |>
  ggplot(aes(x= Nucleus_Area, color = Sample)) +
  geom_density() +
  theme_bw() +
  # scale_color_manual(values = cell_type_colors_halo) +
  facet_wrap(~cell_type, ncol = 1) +
  guides(color=guide_legend(ncol =1))

ggsave(nuc_area_density_sample, filename = here(plot_dir, "Nucleus_Area_density_Sample.png"), width = 12, height = 10)

## violin 
nuc_area_violin <- halo_all |>
  ggplot(aes(x= cell_type, y= Nucleus_Area, fill = cell_type)) +
  geom_violin(draw_quantiles = c(.25,.50,.75)) +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme(legend.position = "None") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "blue") + ##minimum size
  geom_text(label = "Global Min: 5 µm2", y = 1, x = "Astro", color = "blue") +
  geom_hline(yintercept = pi * 5^2, linetype = "dashed", color = "red") + ## max size
  geom_text(label = "Reasonable Area: 78 µm2\n(radius = 5µm)", y = 85, x = "Astro", color = "red") +
  labs(y = "Nucleus Area µm2", x = "Cell Type")

ggsave(nuc_area_violin, filename = here(plot_dir, "Nucleus_Area_violin.png"))

nuc_area_violin_sample <- halo_all |>
  ggplot(aes(x= Sample, y= Nucleus_Area, fill = Sample)) +
  geom_violin(draw_quantiles = c(.25,.50,.75)) +
  theme_bw() +
  facet_wrap(~cell_type, ncol = 1) +
  theme(legend.position = "None", 
        axis.text.x = element_text(angle = 90)) 

ggsave(nuc_area_violin_sample, filename = here(plot_dir, "Nucleus_Area_violin_Sample.png"))



