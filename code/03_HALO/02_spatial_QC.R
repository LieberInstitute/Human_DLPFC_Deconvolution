library("tidyverse") 
library("scales")
library("ggrepel")
library("patchwork")
library("here")
library("sessioninfo")
# library("DeconvoBuddies")

#### Plot Set-up ####
plot_dir <- here("plots", "03_HALO",  "02_spatial_QC")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)

halo_theme <- theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  labs(x = "HALO X Coord", y = "HALO Y Coord") 

blank_axis <- theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)

#### Load Data ###
metadata <- read.csv(here("processed-data","03_HALO","HALO_metadata.csv"))
prop_all <- read.csv(here("processed-data","03_HALO","HALO_cell_type_proportions.csv"))
load(here("processed-data","03_HALO","HALO_Data.Rdata"), verbose = TRUE)

# common_cols <- intersect(colnames(halo_star), colnames(halo_circle))
# 
# halo_all <- rbind(halo_star[,common_cols], halo_circle[,common_cols]) |> left_join(metadata)

dim(halo_all)
# [1] 1936064      34


## get important coordinates for each sample
slide_position <- halo_all |> 
  group_by(Sample_Combo) |>
  summarize(Samp_XMax = max(XMax),
            Samp_XMin = min(XMin),
            Samp_YMax = max(YMax),
            Samp_YMin = min(YMin),
            median_X = median(XMax),
            median_Y = median(YMax)
  )|> 
  left_join(metadata)

star_slide_guide <- slide_position |>
  filter(Combo == "Star") |>
  ggplot(aes(x = median_X, y = median_Y)) +
  geom_rect(aes(xmin=Samp_XMin, xmax=Samp_XMax,
                ymin = Samp_YMin, ymax = Samp_YMax,
                fill = subslide), alpha = 0.2, color = NA) +
  geom_text_repel(aes(label = Sample, color = subslide), size = 2.5)+
  halo_theme+ 
  labs(x = "HALO X Coord", y = "HALO Y Coord", title = "Star Combo") 

ggsave(star_slide_guide + facet_grid(i~Slide), filename = here(plot_dir, "star_slide_guide.png"))
ggsave(star_slide_guide + facet_grid(Round~Slide), filename = here(plot_dir, "star_slide_round_guide.png"))


circle_slide_guide <- slide_position_circle |> 
  filter(Combo == "Circle") |>
  ggplot(aes(x = median_X, y = median_Y)) +
  geom_rect(aes(xmin=Samp_XMin, xmax=Samp_XMax,
                ymin = Samp_YMin, ymax = Samp_YMax,
                fill = subslide), alpha = 0.2, color = NA) +
  geom_text_repel(aes(label = Sample, color = subslide), size = 2.5)+
  halo_theme+ 
  labs(x = "HALO X Coord", y = "HALO Y Coord", title = "Circle Combo") 

ggsave(circle_slide_guide + facet_grid(i~Round), filename = here(plot_dir, "circle_slide_guide.png"))
ggsave(circle_slide_guide + facet_grid(Round~Slide), filename = here(plot_dir, "circle_slide_round_guide.png"))


## Plot Slides
star_slides <- halo_star |> 
  ggplot() +
  geom_rect(aes(xmin=XMin, xmax=XMax,
                ymin = YMin, ymax = YMax, 
                fill = cell_type)) +
  facet_grid(i~Slide) +
  geom_text(data = slide_position_star, aes(x = median_X, y = Samp_YMax + 100, label = Sample), size = 2)+
  scale_fill_manual(values = cell_type_colors_halo) +
  halo_theme+ 
  labs(x = "HALO X Coord", y = "HALO Y Coord", title = "Star Combo") 

ggsave(star_slides, filename = here(plot_dir, "star_slides.png"))


circle_slides <- halo_circle |> 
  ggplot() +
  geom_rect(aes(xmin=XMin, xmax=XMax,
                ymin = YMin, ymax = YMax, 
                fill = cell_type)) +
  geom_text(data = slide_position_circle, aes(x = median_X, y = Samp_YMax + 100, label = Sample), size = 2)+
  facet_grid(i~Slide) +
  scale_fill_manual(values = cell_type_colors_halo) +
  halo_theme+ 
  labs(x = "HALO X Coord", y = "HALO Y Coord", title = "Circle Combo") 

ggsave(circle_slides, filename = here(plot_dir, "circle_slides.png"))

ggsave(star_slides + theme(legend.position = "None") + circle_slides, filename = here(plot_dir, "ALL_slides.png"), width = 12)



star_samples <- halo_star |> 
  ggplot() +
  geom_rect(aes(xmin=XMin, xmax=XMax,
                ymin = YMin, ymax = YMax, 
                fill = cell_type)) +
  facet_wrap(~Sample, scales = 'free_y') +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_halo) 
ggsave(star_samples, filename = here(plot_dir, "star_samples.png"), height = 10, width = 18)


circle_samples <- halo_circle |> 
  ggplot() +
  geom_rect(aes(xmin=XMin, xmax=XMax,
                ymin = YMin, ymax = YMax, 
                fill = cell_type)) +
  facet_wrap(~Sample, scales = 'free_y') +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_halo)

ggsave(circle_samples, filename = here(plot_dir, "circle_samples.png"), height = 10, width = 16)

## plot ct only
star_samples_ct <- halo_star |> 
  filter(cell_type != 'Other') |>
  ggplot() +
  geom_rect(aes(xmin=XMin, xmax=XMax,
                ymin = YMin, ymax = YMax, 
                fill = cell_type)) +
  facet_wrap(~Sample, scales = 'free_y') +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_halo) 

ggsave(star_samples_ct, filename = here(plot_dir, "star_samples_ct.png"), height = 10, width = 18)


circle_samples_ct <- halo_circle |> 
  filter(cell_type != 'Other')|> 
  ggplot() +
  geom_rect(aes(xmin=XMin, xmax=XMax,
                ymin = YMin, ymax = YMax, 
                fill = cell_type)) +
  facet_wrap(~Sample, scales = 'free_y') +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_halo)

ggsave(circle_samples_ct, filename = here(plot_dir, "circle_samples_ct.png"), height = 10, width = 16)


## ALL Samples - need to normalize data
all_samples_ct <- halo_all |> 
  # filter(cell_type != 'Other')|> 
  ggplot() +
  geom_rect(aes(xmin=XMin, xmax=XMax,
                ymin = YMin, ymax = YMax, 
                fill = Combo)) +
  # facet_grid(Combo~Sample, scales = 'free') +
  facet_wrap(~Sample, scales = 'free') +
  # scale_fill_manual(values = cell_type_colors_halo) +
  theme_void() +
  # blank_axis + 
  theme(strip.background = element_blank(),
        strip.text = element_blank())

ggsave(all_samples_ct, filename = here(plot_dir, "ALL_samples_ct.png"), height = 12, width = 16)



#### Hex Plots ####

star_nuc_area_hex <- ggplot(halo_star) +
  stat_summary_hex(aes(x = XMax, y = YMax, z = `Nucleus.Area..µm..`),
                   fun = mean, bins = 100
  ) +
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~Sample, scales = 'free_y') +
  theme_bw() +
  labs(title = "Star Combo", subtitle = "Mean Nucleus Area µm")

ggsave(star_nuc_area_hex, filename = here(plot_dir, "star_hex_nuc_area.png"), height = 10, width = 18)

star_puncta_hex <- 
  # halo_all |> 
  # filter(Combo == "Star") |> 
  ggplot(halo_star) +
  stat_summary_hex(aes(x = XMax, y = YMax, z = `AKT3..Opal.620..Copies`),
                   fun = mean, bins = 100
  ) +
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~Sample_Combo, scales = 'free_y') +
  theme_bw() +
  labs(title = "Star Combo", subtitle = "Mean Number Puncta")

ggsave(star_puncta_hex, filename = here(plot_dir, "star_hex_puncta.png"), height = 10, width = 18)



circle_nuc_area_hex <- ggplot(halo_circle) +
  stat_summary_hex(aes(x = XMax, y = YMax, z = `Nucleus.Area..µm..`),
                   fun = mean, bins = 100
  ) +
  scale_fill_continuous(type = "viridis") +
  facet_wrap(~Sample, scales = 'free_y') +
  theme_bw() +
  labs(title = "Circle Combo", subtitle = "Mean Nucleus Area µm")

ggsave(circle_nuc_area_hex, filename = here(plot_dir, "circle_hex_nuc_area.png"), height = 10, width = 18)

#### print each with grid ####
halo_all$Sample_Combo <- factor(halo_all$Sample_Combo)

walk(levels(halo_all$Sample_Combo), function(s) {
  message(s)
  halo_s <- halo_all %>% filter(Sample_Combo == s)
  
  grid_hex <- ggplot(halo_s) +
    stat_summary_hex(aes(x = XMax, y = YMax, z = `Nucleus.Area..µm..`),
                     fun = mean, bins = 100
    ) +
    scale_fill_continuous(type = "viridis") +
    coord_equal() +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = NA),
      panel.ontop = TRUE,
      panel.grid.major.x = element_line(color = "grey60"),
      panel.grid.major.y = element_line(color = "grey60")
    ) +
    scale_x_continuous(minor_breaks = seq(0, 40000, 1000)) +
    scale_y_continuous(minor_breaks = seq(90000, 0, -1000), trans = "reverse") +
    labs(title = s)
  
  s <- gsub("/", "", s)
  
  ggsave(grid_hex, filename = here(plot_dir, "QC_hex", paste0("grid_hex_", s, ".png")))
})

