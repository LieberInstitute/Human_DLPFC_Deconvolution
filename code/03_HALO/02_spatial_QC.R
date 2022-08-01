library("tidyverse") 
library("scales")
library("here")
library("sessioninfo")
# library("DeconvoBuddies")

#### Plot Set-up ####
plot_dir <- here("plots", "03_HALO",  "02_spatial_QC")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)

#### Load Data ###
metadata <- read.csv(here("processed-data","03_HALO","HALO_metadata.csv"))
prop_all <- read.csv(here("processed-data","03_HALO","HALO_cell_type_proportions.csv"))
load(here("processed-data","03_HALO","HALO_Data.Rdata"), verbose = TRUE)

## Add info
halo_star <- halo_star |> left_join(metadata)

halo_circle <- halo_circle |> left_join(metadata)

## Plot Slides
star_slides <- halo_star |> 
  ggplot() +
  geom_rect(aes(xmin=XMin, xmax=XMax,
                ymin = YMin, ymax = YMax, 
                fill = cell_type)) +
  facet_grid(i~Slide) +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_halo)

ggsave(star_slides, filename = here(plot_dir, "star_slides.png"))


circle_slides <- halo_circle |> 
  ggplot() +
  geom_rect(aes(xmin=XMin, xmax=XMax,
                ymin = YMin, ymax = YMax, 
                fill = cell_type)) +
  facet_grid(i~Slide) +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_halo)

ggsave(circle_slides, filename = here(plot_dir, "circle_slides.png"))



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

## Check out XY locations by Sample
slide_position_star <- halo_star |> 
  group_by(Sample_Combo) |>
  summarize(Samp_XMax = max(XMax),
            Samp_XMin = min(XMin),
            Samp_YMax = max(YMax),
            Samp_YMin = min(XMin),
            median_X = median(XMax),
            median_Y = median(YMax)
            )|> 
  left_join(metadata)

library(ggrepel)

star_slide_guide <- slide_position_star |> 
  ggplot(aes(x = median_X, y = median_Y, color = subslide)) +
  # geom_rect(aes(xmin=Samp_XMin, xmax=Samp_XMax,
  #               ymin = Samp_YMin, ymax = Samp_YMax, 
  #               fill = subslide)) +
  geom_point()+
  geom_text_repel(aes(label = Sample))+
  facet_grid(i~Slide) +
  theme_bw()
         
ggsave(star_slide_guide, filename = here(plot_dir, "star_slide_guide.png"))

star_slide_guide2 <- slide_position_star |> 
  ggplot(aes(x = median_X, y = median_Y, color = subslide)) +
  geom_rect(aes(xmin=Samp_XMin, xmax=Samp_XMax,
                ymin = Samp_YMin, ymax = Samp_YMax,
                fill = subslide), alpha = 0.2) +
  geom_point()+
  facet_wrap(~Sample, scales = 'free_y') +
  theme_bw()

ggsave(star_slide_guide2, filename = here(plot_dir, "star_slide_guide2.png"))

star_samples_ct_outline <-  ggplot() +
  geom_rect(data = halo_star,
            aes(xmin=XMin, xmax=XMax,
                ymin = YMin, ymax = YMax, 
                fill = cell_type)) +
  geom_rect(data = slide_position_star, 
            aes(xmin=Samp_XMin, xmax=Samp_XMax,
                ymin = Samp_YMin, ymax = Samp_YMax,
                fill = subslide), fill = NA, color = "red", linetype = "dashed") +
  facet_wrap(~Sample, scales = 'free_y') +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_halo) 

ggsave(star_samples_ct_outline, filename = here(plot_dir, "star_samples_ct_outline.png"), height = 10, width = 16)

