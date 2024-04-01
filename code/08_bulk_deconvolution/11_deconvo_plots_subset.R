
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")

## prep dirs ##
plot_dir <- here("plots", "08_bulk_deconvolution", "11_deconvo_plots_subset")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors & shapes
load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo
# cell_type_colors_broad
load(here("processed-data","00_data_prep","method_colors.Rdata"), verbose = TRUE)
# method_colors
load(here("processed-data","00_data_prep","library_combo_shapes.Rdata"), verbose = TRUE)
# library_combo_shapes
# library_combo_shapes2

#### load data ####
load(here("processed-data", "08_bulk_deconvolution", "10_get_est_prop_subset","prop_long_subset.Rdata"), verbose = TRUE)
# prop_long_subset
# prop_long_subset_summary_opc
# prop_long_subset_summary

## boxplot across subsets and samples 
all_subsample_boxplot <- prop_long_subset |>
  ggplot(aes(x = SAMPLE_ID, y = prop, color = method)) +
  geom_boxplot() +
  facet_wrap(~cell_type, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
    scale_color_manual(values = method_colors)

ggsave(all_subsample_boxplot, filename = here(plot_dir, "all_subsample_boxplot.png"), height = 10, width = 12)

#### what about median output ? ####
(cor_check <- prop_long_subset_summary |>
    filter(!is.na(RNAscope_prop)) |>
    group_by(method) |>
    summarize(cor = cor(RNAscope_prop, median_prop),
              rmse = Metrics::rmse(RNAscope_prop, median_prop))  |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
    arrange(cor))

# method     cor  rmse cor_anno                
# <chr>    <dbl> <dbl> <chr>                   
# 1 Bisque -0.0185 0.135 "cor:-0.018\nrmse:0.135"
# 2 hspe    0.465  0.137 "cor:0.465\nrmse:0.137" 

(cor_check_library <- prop_long_subset_summary |>
    filter(!is.na(RNAscope_prop)) |>
    group_by(method, library_type, rna_extract, library_combo) |>
    summarize(cor = cor(RNAscope_prop, median_prop),
              rmse = Metrics::rmse(RNAscope_prop, median_prop))  |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
    arrange(cor))

## facet prop
prop_bar_SAMPLE_facet <- prop_long_subset_summary_opc |> 
  mutate(Sample = gsub("_","\n", Sample)) |>
  ggplot(aes(x = library_combo, y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_grid(method~Sample) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Mean Estimated Proportion (1k reps)", x = "Library Type & RNA Extraction Prep", fill = "Cell Type", title = "Subset Experiment") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet_subset.png"), width = 12, height = 5)
ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet_subset.pdf"), width = 12, height = 5)
## scater plot
est_prop_v_RNAscope_scatter <- prop_long_subset_summary |>
  filter(!is.na(RNAscope_prop)) |>
  ggplot() +
  geom_point(aes(x = RNAscope_prop, y = mean_prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check, 
            aes(label = cor_anno,x = .5, y = .6),
            vjust = "inward", hjust = "inward") +
  facet_wrap(~method, nrow = 1) +
  scale_color_manual(values = cell_type_colors_halo) +
  scale_shape_manual(values = library_combo_shapes2) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Mean Estimated Proportion (1k reps)")

ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_subset.png"), height = 4)
ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_subset.pdf"))

#### compare to full data set ####

load(here("processed-data", "08_bulk_deconvolution", "03_get_est_prop","prop_long.Rdata"), verbose = TRUE)
prop_long_opc <- prop_long_opc |> 
  filter(marker == "MeanRatio_top25") |>
  select(SAMPLE_ID, method, cell_type, prop)

prop_long_subset_summary_opc <- prop_long_subset_summary_opc|>
  left_join(prop_long_opc)

(cor_check_full <- prop_long_subset_summary_opc |>
    group_by(method) |>
    summarize(cor = cor(prop, median_prop),
              rmse = Metrics::rmse(prop, median_prop))  |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
    arrange(cor))

# method   cor  rmse cor_anno               
# <chr>  <dbl> <dbl> <chr>                  
# 1 Bisque 0.430 0.135 "cor:0.430\nrmse:0.135"
# 2 hspe   0.927 0.137 "cor:0.927\nrmse:0.137"

est_mean_prop_v_prop_scatter <- prop_long_subset_summary_opc |>
  ggplot() +
  geom_point(aes(x = prop, y = mean_prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check_full, 
            aes(label = cor_anno,x = .25, y = .75),
            vjust = "inward") +
  facet_wrap(~method, nrow = 1) +
  scale_color_manual(values = cell_type_colors_broad) +
  scale_shape_manual(values = library_combo_shapes2) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "Full Data Estimated Proportion", y = "Mean Estimated Proportion (1k reps)")

ggsave(est_mean_prop_v_prop_scatter, filename = here(plot_dir, "est_mean_prop_v_prop_scatter_subset.png"), height = 4)
ggsave(est_mean_prop_v_prop_scatter, filename = here(plot_dir, "est_mean_prop_v_prop_scatter_subset.pdf"))

