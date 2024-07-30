
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")

## prep dirs ##
plot_dir <- here("plots", "13_PEC_deconvolution", "11_deconvo_plots_donor_subset")
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
load(here("processed-data", "13_PEC_deconvolution", "10_get_est_prop_donor_subset","prop_long_donor_subset.Rdata"), verbose = TRUE)
# prop_long_donor_subset
# prop_long_donor_subset_summary_opc
# prop_long_donor_subset_summary

## boxplot across subsets and samples 
all_subsample_boxplot <- prop_long_donor_subset |>
  ggplot(aes(x = SAMPLE_ID, y = prop, color = method)) +
  geom_boxplot() +
  facet_grid(num_donors~cell_type, scales = "free") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
    scale_color_manual(values = method_colors)

ggsave(all_subsample_boxplot, filename = here(plot_dir, "donor_subset_all_subsample_boxplot.png"), height = 10, width = 12)

#### cor vs. mean prop ####
(cor_check <- prop_long_donor_subset_summary |>
    filter(!is.na(RNAscope_prop)) |>
    group_by(method, num_donors) |>
    summarize(cor = cor(RNAscope_prop, mean_prop),
              rmse = Metrics::rmse(RNAscope_prop, median_prop))  |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
    arrange(-cor))


(cor_check_library <- prop_long_donor_subset_summary |>
    filter(!is.na(RNAscope_prop)) |>
    group_by(method,num_donors, library_type, rna_extract, library_combo) |>
    summarize(cor = cor(RNAscope_prop, mean_prop),
              rmse = Metrics::rmse(RNAscope_prop, median_prop))  |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
    arrange(-cor))

prop_long_donor_subset |> head()
  

## facet prop
# prop_bar_SAMPLE_facet <- prop_long_donor_subset_summary_opc |> 
#   mutate(Sample = gsub("_","\n", Sample)) |>
#   ggplot(aes(x = library_combo, y = mean_prop, fill = cell_type)) +
#   geom_bar(stat = "identity") +
#   facet_grid(method~Sample) +
#   scale_fill_manual(values = cell_type_colors_broad) +
#   labs(y = "Mean Estimated Proportion (1k reps)", x = "Library Type & RNA Extraction Prep", fill = "Cell Type", title = "Subset Experiment") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         legend.position = "bottom")
# 
# ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet_subset.png"), width = 12, height = 5)
# ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet_subset.pdf"), width = 12, height = 5)

## scater plot
est_prop_v_RNAscope_scatter <- prop_long_donor_subset_summary |>
  filter(!is.na(RNAscope_prop)) |>
  ggplot() +
  geom_point(aes(x = RNAscope_prop, y = mean_prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check, 
            aes(label = cor_anno,x = .5, y = .7),
            vjust = "inward", hjust = "inward") +
  # facet_grid(num_donors~method) +
  facet_grid(method~num_donors) +
  scale_color_manual(values = cell_type_colors_halo) +
  scale_shape_manual(values = library_combo_shapes2) +
  geom_abline() +
  # coord_equal() +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs( x = "RNAscope Proportion", y = "Mean Estimated Proportion (100 reps)")

ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_donor_subset.png"), width = 12)
ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_donor_subset.pdf"), width = 12)

#### Cor over num dononrs ####

cor_rmse_line_n_donors <- cor_check|>
  ggplot(aes(x = num_donors, y = cor, color = method, group = method)) +
  geom_point(aes(size = rmse), alpha = .7) +
  geom_line() +
  scale_size(range = c(1,8)) +
  theme_bw() +
  labs(x = "Library Type + RNA Extraction") +
  scale_color_manual(values = method_colors)
 
ggsave(cor_rmse_line_n_donors, filename = here(plot_dir, "cor_rmse_line_n_donors.png"), height = 5)


cor_rmse_line_n_donors_library <- cor_check_library|>
  ggplot(aes(x = num_donors, y = cor, color = library_combo, group = library_combo)) +
  geom_point(aes(size = rmse), alpha = .7) +
  geom_line() +
  scale_size(range = c(1,8)) +
  facet_wrap(~method) +
  theme_bw() +
  labs(x = "Library Type + RNA Extraction") +
  theme(legend.position = "bottom")

ggsave(cor_rmse_line_n_donors_library, filename = here(plot_dir, "cor_rmse_line_n_donors_library.png"), width = 12)

#### cor distribution vs. num dononrs ####

cor_subset_boxplot <- cor_check_subset |>
  ggplot(aes(x= as.factor(num_donors), y = cor, color = method)) +
  geom_boxplot() +
  facet_wrap(~method) +
  scale_color_manual(values = method_colors) +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "n donors")

ggsave(cor_subset_boxplot, filename = here(plot_dir, "donor_subset_boxplot_cor.png"), width = 12, height = 4)

rmse_subset_boxplot <- cor_check_subset |>
  ggplot(aes(x= as.factor(num_donors), y = rmse, color = method)) +
  geom_boxplot() +
  facet_wrap(~method) +
  scale_color_manual(values = method_colors) +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "n donors")

ggsave(rmse_subset_boxplot, filename = here(plot_dir, "donor_subset_boxplot_rmse.png"), width = 12, height = 4)



