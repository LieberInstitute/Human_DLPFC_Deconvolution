
library("tidyverse")
library("scales")
library("patchwork")
library("here")
library("sessioninfo")

#### Dir Set-up ####
plot_dir <- here("plots", "03_HALO", "11_HALO_hex_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

#### Load Data ####
load(here("processed-data", "03_HALO", "halo_all.Rdata"), verbose = TRUE)

## filter out large nuclei & low confidence sections
halo_all <- halo_all |> 
  filter(!large_nuc,
         Confidence %in% c("OK","High"))

halo_all |> 
  group_by(Combo, Sample, Confidence) |> 
  dplyr::count() |>
  print(n = 25)

# Combo  Sample      Confidence     n
# <chr>  <chr>       <fct>      <int>
#   1 Circle Br2720_post High       36173
# 2 Circle Br2743_ant  OK         57674
# 3 Circle Br3942_ant  High       49580
# 4 Circle Br3942_mid  High       33383
# 5 Circle Br3942_post High       50920
# 6 Circle Br6423_ant  OK         45362
# 7 Circle Br6423_post High       40035
# 8 Circle Br6471_mid  High       47550
# 9 Circle Br8325_mid  OK         32425
# 10 Circle Br8492_mid  High       44308
# 11 Circle Br8492_post High       53008
# 12 Circle Br8667_mid  High       38433
# 13 Star   Br2720_post High       37553
# 14 Star   Br2743_ant  High       49866
# 15 Star   Br3942_mid  OK         34003
# 16 Star   Br6423_ant  OK         46796
# 17 Star   Br6423_post High       41289
# 18 Star   Br6432_post OK         28703
# 19 Star   Br6471_mid  OK         53299
# 20 Star   Br6522_mid  OK         36592
# 21 Star   Br6522_post OK         28093
# 22 Star   Br8325_ant  OK         31145
# 23 Star   Br8325_mid  OK         36050
# 24 Star   Br8492_post High       53709
# 25 Star   Br8667_mid  OK         39648

## example: Br8667_mid 

colnames(halo_all)

# halo_all <- halo_all 

coord_adj <- halo_all |>
  group_by(SAMPLE_ID) |>
  summarize(zeroX = min(XMax),
             zeroY = min(YMax))

halo_all_adj <- halo_all |>
  select(Combo, SAMPLE_ID, Sample, BrNum, Object_Id, XMin, XMax, YMin, YMax, cell_type, AKT3_Copies, Nucleus_Area)|>
  left_join(coord_adj) |>
  mutate(XMax = XMax - zeroX,
         YMax = YMax - zeroY,
         Sample = fct_rev(Sample))

blank_axis <- theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)

#### Star plots ####
star_akt3 <- halo_all_adj |>
  filter(Combo == "Star") |>
  mutate(title = "AKT3 Copies") |>
  ggplot() +
  stat_summary_hex(aes(x = XMax, y = YMax, z = `AKT3_Copies`),
                   fun = mean, bins = 100
  )  +
  scale_fill_continuous(type = "viridis", name = "Mean Copies") + 
  coord_equal() +
  theme_bw() +
  facet_grid(Sample~title) +
  theme(legend.position = "bottom") +
  blank_axis

ggsave(star_akt3, filename = here(plot_dir, "HALO_star_akt3.png"), height = 12, width = 2.5)
ggsave(star_akt3, filename = here(plot_dir, "HALO_star_akt3.pdf"), height = 12, width = 2.5)

# strip.background = element_blank(),
# strip.text.x = element_blank()

star_area <- halo_all_adj |>
  filter(Combo == "Star") |>
  mutate(title = "Nuclear Area") |>
  ggplot() +
  stat_summary_hex(aes(x = XMax, y = YMax, z = `Nucleus_Area`),
                   fun = mean, bins = 100
  )  +
  scale_fill_continuous(type = "viridis", name = "Mean Area") + 
  coord_equal() +
  theme_bw() +
  facet_grid(Sample~title) +
  theme(legend.position = "bottom") +
  blank_axis

ggsave(star_area, filename = here(plot_dir, "HALO_star_area.png"), height = 12, width = 2.5)
ggsave(star_area, filename = here(plot_dir, "HALO_star_area.pdf"), height = 12, width = 2.5)

star_cells <- halo_all_adj |>
  filter(Combo == "Star",
         cell_type != "Other") |>
  ggplot() +
  geom_hex(aes(x = XMax, y = YMax), bins = 100)  +
  scale_fill_continuous(type = "viridis", name = "Number of Cells") + ## top value of 100 for visualization 
  coord_equal() +
  theme_bw() +
  facet_grid(Sample~cell_type) +
  theme(legend.position = "bottom") +
  blank_axis

ggsave(star_cells, filename = here(plot_dir, "HALO_star_cells.png"), height = 12, width = 3.5)
ggsave(star_cells, filename = here(plot_dir, "HALO_star_cells.pdf"), height = 12, width = 3.5)

# ggsave(star_akt3 + star_cells, filename = here(plot_dir, "HALO_star.png"), height = 12, width = 6)

#### Circle plots ####
circle_akt3 <- halo_all_adj |>
  filter(Combo == "Circle") |>
  mutate(title = "AKT3 Copies") |>
  ggplot() +
  stat_summary_hex(aes(x = XMax, y = YMax, z = `AKT3_Copies`),
                   fun = mean, bins = 100
  )  +
  scale_fill_continuous(type = "viridis", name = "Mean Copies") + 
  coord_equal() +
  theme_bw() +
  facet_grid(Sample~title) +
  theme(legend.position = "bottom") +
  blank_axis

ggsave(circle_akt3, filename = here(plot_dir, "HALO_circle_akt3.png"), height = 12, width = 2.5)
ggsave(circle_akt3, filename = here(plot_dir, "HALO_circle_akt3.pdf"), height = 12, width = 2.5)

# strip.background = element_blank(),
# strip.text.x = element_blank()

circle_area <- halo_all_adj |>
  filter(Combo == "Circle") |>
  mutate(title = "Nuclear Area") |>
  ggplot() +
  stat_summary_hex(aes(x = XMax, y = YMax, z = `Nucleus_Area`),
                   fun = mean, bins = 100
  )  +
  scale_fill_continuous(type = "viridis", name = "Mean Area") + 
  coord_equal() +
  theme_bw() +
  facet_grid(Sample~title) +
  theme(legend.position = "bottom") +
  blank_axis

ggsave(circle_area, filename = here(plot_dir, "HALO_circle_area.png"), height = 12, width = 2.5)
ggsave(circle_area, filename = here(plot_dir, "HALO_circle_area.pdf"), height = 12, width = 2.5)

circle_cells <- halo_all_adj |>
  filter(Combo == "Circle",
         cell_type != "Other") |>
  ggplot() +
  geom_hex(aes(x = XMax, y = YMax), bins = 100)  +
  scale_fill_continuous(type = "viridis", name = "Number of Cells") + ## top value of 100 for visualization 
  coord_equal() +
  theme_bw() +
  facet_grid(Sample~cell_type) +
  theme(legend.position = "bottom") +
  blank_axis

ggsave(circle_cells, filename = here(plot_dir, "HALO_circle_cells.png"), height = 12, width = 3.5)
ggsave(circle_cells, filename = here(plot_dir, "HALO_circle_cells.pdf"), height = 12, width = 3.5)


