library("tidyverse")
library("SingleCellExperiment")
library("here")
library("sessioninfo")
library("here")
library("patchwork")

#### Set-up ####
plot_dir <- here("plots", "03_HALO", "12_HALO_cell_size")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "03_HALO", "12_HALO_cell_size")
if (!dir.exists(data_dir)) dir.create(data_dir)

## colors
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo

#### load HALO data ####

load(here("processed-data", "03_HALO", "halo_all.Rdata"), verbose = TRUE)

halo_all <- halo_all |>
  select(Sample, cell_type, AKT3_Copies, Nucleus_Area, large_nuc) |>
  mutate(cell_type = factor(cell_type, levels = names(cell_type_colors_halo))) |>
  filter(cell_type != "Other")

#### stats not filtering for large nuc ####

halo_all |>
  group_by(cell_type) |>
  summarise(median_area = median(Nucleus_Area),
            mean_area = mean(Nucleus_Area),
            median_copies = median(AKT3_Copies),
            mean_copies = mean(AKT3_Copies),
            median_area_akt3 = median(Nucleus_Area * AKT3_Copies))

# cell_type median_area mean_area median_copies mean_copies
# <fct>           <dbl>     <dbl>         <dbl>       <dbl>
# 1 Astro            30.8      33.6             3        4.81
# 2 EndoMural        33.1      36.3             3        3.91
# 3 Micro            28.3      29.3             2        3.01
# 4 OligoOPC         29.2      31.4             3        3.18
# 5 Excit            36.5      39.7            11       12.3 
# 6 Inhib            37.8      40.1             8        8.47

(halo_cell_size <- halo_all |>
  filter(!large_nuc) |>
  group_by(cell_type) |>
  summarise(median_area = median(Nucleus_Area),
            mean_area = mean(Nucleus_Area),
            median_copies = median(AKT3_Copies),
            mean_copies = mean(AKT3_Copies),
            median_area_akt3 = median(Nucleus_Area * AKT3_Copies)))

# cell_type median_area mean_area median_copies mean_copies median_area_akt3
# <fct>           <dbl>     <dbl>         <dbl>       <dbl>            <dbl>
# 1 Astro            30.2      31.4             3        4.29             76.5
# 2 EndoMural        32.4      34.0             3        3.77            102. 
# 3 Micro            28.2      28.7             2        2.97             61.4
# 4 OligoOPC         29.0      30.4             3        3.11             68.0
# 5 Excit            35.0      36.2            10       11.5             352. 
# 6 Inhib            36.8      37.4             7        8.11            260. 

## export for MuSiC
MuSiC_cell_size_input <- halo_cell_size |>
  select(cell_type, 
         nuc_area = median_area,
         akt3 = median_copies,
         nuc_area_akt3 =  median_area_akt3) 

write_csv(MuSiC_cell_size_input, file = here(data_dir, "MuSiC_cell_size.csv"))

## unfiltered plots
nuc_area_boxplot_unfiltered <- halo_all |>
  ggplot(aes(x = cell_type, y = Nucleus_Area, fill = cell_type)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw() +
  geom_hline(yintercept = pi * 5^2, color= "red", linetype = "dashed")

ggsave(nuc_area_boxplot_unfiltered, filename = here(plot_dir, "nuc_area_boxplot_unfiltered.png"))

akt3_boxplot_unfiltered <- halo_all |>
  ggplot(aes(x = cell_type, y = AKT3_Copies, fill = cell_type, color = large_nuc)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw()

ggsave(akt3_boxplot_unfiltered, filename = here(plot_dir, "akt3_boxplot_unfiltered.png"))

## filter plots
nuc_area_boxplot<- halo_all |>
  filter(!large_nuc) |>
  ggplot(aes(x = cell_type, y = Nucleus_Area, fill = cell_type)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw()

ggsave(nuc_area_boxplot, filename = here(plot_dir, "nuc_area_boxplot.png"))

akt3_boxplot <- halo_all |>
  filter(!large_nuc) |>
  ggplot(aes(x = cell_type, y = AKT3_Copies, fill = cell_type)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw()

ggsave(akt3_boxplot, filename = here(plot_dir, "akt3_boxplot.png"))

areaXakt3_boxplot <- halo_all |>
  filter(!large_nuc) |>
  ggplot(aes(x = cell_type, y = AKT3_Copies * Nucleus_Area, fill = cell_type)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors_halo) +
  theme_bw()

ggsave(areaXakt3_boxplot, filename = here(plot_dir, "areaXakt3_boxplot.png"))

cell_size_boxplots <- nuc_area_boxplot + theme(legend.position = "None",axis.text.x = element_text(angle = 45, hjust=1)) + 
  akt3_boxplot + theme(legend.position = "None",axis.text.x = element_text(angle = 45, hjust=1)) +
  areaXakt3_boxplot + theme(legend.position = "None",axis.text.x = element_text(angle = 45, hjust=1))

ggsave(cell_size_boxplots, filename = here(plot_dir, "cell_size_boxplots.png"), width = 10, height = 5)

#### scatter plot ####
area_akt3_scatter <- halo_all |>
  # filter(!large_nuc) |>
  ggplot(aes(x = Nucleus_Area, y = AKT3_Copies, color = cell_type)) +
  geom_point(alpha = .5) +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_wrap(~cell_type) +
  theme_bw() +
  geom_vline(xintercept = pi * 5^2, color= "red", linetype = "dashed")

ggsave(area_akt3_scatter, filename = here(plot_dir, "area_akt3_scatter.png"))


