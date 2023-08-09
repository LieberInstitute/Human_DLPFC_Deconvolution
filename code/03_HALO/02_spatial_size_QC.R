library("tidyverse")
library("scales")
library("ggrepel")
library("patchwork")
library("broom")
library("here")
library("sessioninfo")
# library("DeconvoBuddies")

#### Dir Set-up ####
plot_dir <- here("plots", "03_HALO", "02_spatial_QC")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "03_HALO", "02_spatial_QC")
if (!dir.exists(data_dir)) dir.create(data_dir)

#### Plot Set-up ####
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo

halo_theme <- theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

blank_axis <- theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
)

#### Load Data ####
metadata <- read.csv(here("processed-data", "03_HALO", "HALO_metadata.csv"))
# prop_all <- read.csv(here("processed-data","03_HALO","HALO_cell_type_proportions.csv"))
load(here("processed-data", "03_HALO", "halo_all.Rdata"), verbose = TRUE)

dim(halo_all)
# [1] 1685807      40

## get important coordinates for each sample
slide_position <- halo_all |>
    group_by(Sample, SAMPLE_ID, Combo, subslide, i, Slide) |>
    summarize(
        Samp_XMax = max(XMax),
        Samp_XMin = min(XMin),
        Samp_YMax = max(YMax),
        Samp_YMin = min(YMin),
        median_X = median(XMax),
        median_Y = median(YMax)
    )

star_slide_guide <- slide_position |>
    filter(Combo == "Star") |>
    ggplot(aes(x = median_X, y = median_Y)) +
    geom_rect(aes(
        xmin = Samp_XMin, xmax = Samp_XMax,
        ymin = Samp_YMin, ymax = Samp_YMax,
        fill = subslide
    ), alpha = 0.2, color = NA) +
    geom_text_repel(aes(label = Sample, color = subslide), size = 2.5) +
    halo_theme +
    labs(x = "HALO X Coord", y = "HALO Y Coord", title = "Star Combo")

ggsave(star_slide_guide + facet_grid(i ~ Slide), filename = here(plot_dir, "slide_guide_star.png"))
# ggsave(star_slide_guide + facet_grid(Round~Slide), filename = here(plot_dir, "star_slide_round_guide.png"))


circle_slide_guide <- slide_position |>
    filter(Combo == "Circle") |>
    ggplot(aes(x = median_X, y = median_Y)) +
    geom_rect(aes(
        xmin = Samp_XMin, xmax = Samp_XMax,
        ymin = Samp_YMin, ymax = Samp_YMax,
        fill = subslide
    ), alpha = 0.2, color = NA) +
    geom_text_repel(aes(label = Sample, color = subslide), size = 2.5) +
    halo_theme +
    labs(x = "HALO X Coord", y = "HALO Y Coord", title = "Circle Combo")

ggsave(circle_slide_guide + facet_grid(i ~ Slide), filename = here(plot_dir, "slide_guide_circle.png"))
# ggsave(circle_slide_guide + facet_grid(Round~Slide), filename = here(plot_dir, "circle_slide_round_guide.png"))


## Plot Slides
all_slides <- halo_all |>
    ggplot() +
    geom_rect(aes(
        xmin = XMin, xmax = XMax,
        ymin = YMin, ymax = YMax,
        fill = cell_type
    )) +
    facet_grid(i ~ Slide + Combo) +
    geom_text(data = slide_position, aes(x = median_X, y = Samp_YMax + 100, label = SAMPLE_ID), size = 2) +
    scale_fill_manual(values = cell_type_colors_halo) +
    halo_theme +
    labs(x = "HALO X Coord", y = "HALO Y Coord", title = "Star Combo")

ggsave(all_slides, filename = here(plot_dir, "slides_all.png"), width = 12)


star_samples <- halo_all |>
    filter(Combo == "Star") |>
    ggplot() +
    geom_rect(aes(
        xmin = XMin, xmax = XMax,
        ymin = YMin, ymax = YMax,
        fill = cell_type
    )) +
    facet_wrap(~Sample, scales = "free_y") +
    theme_bw() +
    scale_fill_manual(values = cell_type_colors_halo)

ggsave(star_samples, filename = here(plot_dir, "nuc_samples_Star.png"), height = 10, width = 18)
ggsave(star_samples, filename = here(plot_dir, "nuc_samples_Star.pdf"), height = 10, width = 18)


circle_samples <- halo_all |>
    filter(Combo == "Circle") |>
    ggplot() +
    geom_rect(aes(
        xmin = XMin, xmax = XMax,
        ymin = YMin, ymax = YMax,
        fill = cell_type
    )) +
    facet_wrap(~Sample, scales = "free_y") +
    theme_bw() +
    scale_fill_manual(values = cell_type_colors_halo)

ggsave(circle_samples, filename = here(plot_dir, "nuc_samples_Circle.png"), height = 10, width = 16)
ggsave(circle_samples, filename = here(plot_dir, "nuc_samples_Circle.pdf"), height = 10, width = 16)

## plot ct only
star_samples_ct <- halo_all |>
    filter(Combo == "Star", cell_type != "Other") |>
    ggplot() +
    geom_rect(aes(
        xmin = XMin, xmax = XMax,
        ymin = YMin, ymax = YMax,
        fill = cell_type
    )) +
    facet_wrap(~Sample, scales = "free_y") +
    theme_bw() +
    scale_fill_manual(values = cell_type_colors_halo)

ggsave(star_samples_ct, filename = here(plot_dir, "nuc_ct_samples_star.png"), height = 10, width = 18)
ggsave(star_samples_ct, filename = here(plot_dir, "nuc_ct_samples_star.pdf"), height = 10, width = 18)


circle_samples_ct <- halo_all |>
    filter(Combo == "Circle", cell_type != "Other") |>
    ggplot() +
    geom_rect(aes(
        xmin = XMin, xmax = XMax,
        ymin = YMin, ymax = YMax,
        fill = cell_type
    )) +
    facet_wrap(~Sample, scales = "free_y") +
    theme_bw() +
    scale_fill_manual(values = cell_type_colors_halo)

ggsave(circle_samples_ct, filename = here(plot_dir, "nuc_ct_samples_circle.png"), height = 10, width = 16)
ggsave(circle_samples_ct, filename = here(plot_dir, "nuc_ct_samples_circle.pdf"), height = 10, width = 16)


## ALL Samples
all_samples_ct <- halo_all |>
    # filter(cell_type != 'Other')|>
    ggplot() +
    geom_rect(aes(
        xmin = XMin, xmax = XMax,
        ymin = YMin, ymax = YMax,
        fill = Combo
    )) +
    # facet_grid(Combo~Sample, scales = 'free') +
    facet_wrap(~Sample, scales = "free") +
    # scale_fill_manual(values = cell_type_colors_halo) +
    theme_void() +
    # blank_axis +
    theme(
        strip.background = element_blank(),
        strip.text = element_blank()
    )

ggsave(all_samples_ct, filename = here(plot_dir, "nuc_samples_all.pdf"), height = 12, width = 16)
ggsave(all_samples_ct, filename = here(plot_dir, "nuc_samples_all.png"), height = 12, width = 16)

#### nuc plots by sample ####

## annotate nuclei bigger than biologically reasonable 
halo_all <- halo_all |> mutate(large_nuc = Nucleus_Area > pi * 5^2) ## larger 5µm radius

plot_dir_sample <- here("plots", "03_HALO", "02_spatial_QC", "Sample_Nuc_plots")
if (!dir.exists(plot_dir_sample)) dir.create(plot_dir_sample)

walk(unique(halo_all$Sample), function(s){
  message(s)
  nuc_plot <- halo_all |>
    dplyr::filter(Sample == s)|>
    ggplot() +
    geom_rect(aes(
      xmin = XMin, xmax = XMax,
      ymin = YMin, ymax = YMax,
      fill = large_nuc
    )) +
    facet_wrap(~Combo) +
    scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "black")) +
    coord_equal() +
    theme_bw()
  
  ggsave(nuc_plot, filename = here(plot_dir_sample, paste0("nuc_",s,".png")), width = 10)
})


#### Hex Plots ####

star_nuc_area_hex <- ggplot(halo_star) +
    stat_summary_hex(aes(x = XMax, y = YMax, z = `Nucleus.Area..µm..`),
        fun = mean, bins = 100
    ) +
    scale_fill_continuous(type = "viridis") +
    facet_wrap(~Sample, scales = "free_y") +
    theme_bw() +
    labs(title = "Star Combo", subtitle = "Mean Nucleus Area µm")

ggsave(star_nuc_area_hex, filename = here(plot_dir, "star_hex_nuc_area.png"), height = 10, width = 18)

#### print each with grid ####
plot_dir_hex <- here("plots", "03_HALO", "02_spatial_QC", "QC_hex")
if (!dir.exists(plot_dir_hex)) dir.create(plot_dir_hex)

halo_all$SAMPLE_ID <- factor(halo_all$SAMPLE_ID)

walk(unique(halo_all$Sample), function(s) {
    message(s)
    halo_s <- halo_all %>% filter(Sample == s)

    area_hex <- ggplot(halo_s) +
        stat_summary_hex(aes(x = XMax, y = YMax, z = `Nucleus_Area`),
            fun = mean, bins = 100
        ) +
        scale_fill_continuous(type = "viridis", limits = c(5,100), "Mean Nuc Area\n(max 100)") + ## top value of 100 for visualization 
        coord_equal() +
        theme_bw() +
        # theme(
        #     panel.background = element_rect(fill = NA),
        #     panel.ontop = TRUE,
        #     panel.grid.major.x = element_line(color = "grey60"),
        #     panel.grid.major.y = element_line(color = "grey60")
        # ) +
        # scale_x_continuous(minor_breaks = seq(0, 40000, 1000)) +
        # scale_y_continuous(minor_breaks = seq(90000, 0, -1000), trans = "reverse") +
        facet_wrap(~Combo) +
        labs(title = s)

    # s <- gsub("/", "", s)

    ggsave(area_hex, filename = here(plot_dir_hex,  paste0("hex_MeanNucArea_", s, ".png")), width = 10)
})

## plot mean number of puncta
walk(unique(halo_all$Sample), function(s) {
  message(s)
  halo_s <- halo_all %>% filter(Sample == s)
  
  puncta_hex <- ggplot(halo_s) +
    stat_summary_hex(aes(x = XMax, y = YMax, z = AKT3_Copies),
                     fun = mean, bins = 100
    ) +
    scale_fill_continuous(type = "viridis", name = "Mean AKT3 Copies") + 
    coord_equal() +
    theme_bw() +
    facet_wrap(~Combo) +
    labs(title = s)
  
  # s <- gsub("/", "", s)
  
  ggsave(puncta_hex, filename = here(plot_dir_hex,  paste0("hex_MeanPuncta_", s, ".png")), width = 10)
})


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


## What samples have the most large nuclei?
halo_all |> count(large_nuc)
# large_nuc       n
# 1 FALSE     1631474
# 2 TRUE        54333


filter_n_cells <- halo_all |>
  count(Sample, Combo)|>
  pivot_wider(names_from = "Combo", values_from = "n") |>
  mutate(filter = "All Nuc") |>
  bind_rows(halo_all |>
              filter(!large_nuc) |>
              count(Sample, Combo)|>
              pivot_wider(names_from = "Combo", values_from = "n") |>
              mutate(filter = "NA < 78")) |>
  mutate(error = abs(Star - Circle)/Circle) 

combo_n_cells_scatter_filter <- filter_n_cells |>
  ggplot(aes(Circle, Star)) +
  geom_line(aes(group = Sample), color = "grey50") +
  geom_point(aes(color = filter)) +
  geom_text_repel(aes(label = ifelse(filter == "All Nuc",Sample,"")), size = 2) +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw() 

ggsave(combo_n_cells_scatter_filter, filename = here(plot_dir, "halo_combo_n_cells_scatter_filter.png"))


#### Puncta vs. Nuc Area ####

puncta_v_size <- halo_all %>%
  ungroup() %>%
  group_by(cell_type) %>%
  do(fitPuncta = tidy(lm(AKT3_Copies ~ Nucleus_Area, data = .))) %>%
  unnest(fitPuncta) %>%
  filter(term == "Nucleus_Area") %>%
  mutate(
    p.bonf = p.adjust(p.value, "bonf"),
    p.bonf.sig = p.bonf < 0.05,
    p.bonf.cat = cut(p.bonf,
                     breaks = c(1, 0.05, 0.01, 0.005, 0),
                     labels = c("<= 0.005", "<= 0.01", "<= 0.05", "> 0.05"),
                     include.lowest = TRUE
    ),
    p.fdr = p.adjust(p.value, "fdr"),
    log.p.bonf = -log10(p.bonf)
  )

puncta_v_size_anno <- puncta_v_size %>%
  select(cell_type, estimate, std.error) %>%
  mutate(anno = paste(
    "beta ==", round(estimate, 3)
    # ,"\nse ==", formatC(std.error, format = "e", digits = 1)
  ))


## plot
puncta_NucArea_scater <- halo_all |>
  ggplot(aes(Nucleus_Area, AKT3_Copies, color = cell_type)) +
  geom_point(aes(color = cell_type), size = 0.2, alpha = 0.2) +
  geom_smooth(method = "lm", color = "black") +
  geom_text(data = puncta_v_size_anno, aes(label = anno), parse = TRUE, x = 45, y = 60, size = 3) +
  facet_wrap(~cell_type, nrow = 2) +
  scale_color_manual(values = cell_type_colors_halo) +
  labs(
    x = bquote("Nucleus Area µm"^2),
    y = "Number of Puncta"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  )

ggsave(puncta_NucArea_scater, filename = here(plot_dir, "puncta_NucArea_scater.png"))


puncta_NucArea_hex <- halo_all |>
  filter(cell_type != "Other") |>
  ggplot(aes(Nucleus_Area, AKT3_Copies)) +
  geom_hex() +
  scale_fill_continuous(type = "viridis", trans = "log") +
  geom_smooth(method = "lm", color = "black") +
  # geom_text(data = puncta_v_size_anno, aes(label = anno), parse = TRUE, x = 45, y = 60, size = 3) +
  facet_wrap(~cell_type, nrow = 2) +
  # scale_color_manual(values = cell_type_colors_halo) +
  labs(
    x = bquote("Nucleus Area µm"^2),
    y = "Number of Puncta"
  ) +
  theme_bw() 

ggsave(puncta_NucArea_hex, filename = here(plot_dir, "puncta_NucArea_hex.png"))

puncta_NucArea_hex_filter <- halo_all |>
  filter(cell_type != "Other", !large_nuc) |>
  ggplot(aes(Nucleus_Area, AKT3_Copies)) +
  geom_hex() +
  scale_fill_continuous(type = "viridis", trans = "log") +
  geom_smooth(method = "lm", color = "black") +
  # geom_text(data = puncta_v_size_anno, aes(label = anno), parse = TRUE, x = 45, y = 60, size = 3) +
  facet_wrap(~cell_type, nrow = 2) +
  # scale_color_manual(values = cell_type_colors_halo) +
  labs(
    x = bquote("Nucleus Area µm"^2),
    y = "Number of Puncta"
  ) +
  theme_bw() 

ggsave(puncta_NucArea_hex_filter, filename = here(plot_dir, "puncta_NucArea_hex_filter.png"))



# sgejobs::job_single('07_TREG_boxplots', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 07_TREG_boxplots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

