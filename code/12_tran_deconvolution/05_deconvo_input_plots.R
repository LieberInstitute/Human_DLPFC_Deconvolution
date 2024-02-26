
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")
library("viridis")
library("GGally")

plot_dir <- here("plots", "12_tran_deconvolution", "05_deconvo_input_plots")
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

#### load data - combine input datasets ####
load(here("processed-data", "08_bulk_deconvolution", "03_get_est_prop","prop_long.Rdata"), verbose = TRUE)
prop_long_paired <- prop_long |>
  filter(method %in% c("Bisque", "hspe"),
         marker == "MeanRatio_top25") |>
  mutate(input = "paired")
  
load(here("processed-data", "12_tran_deconvolution", "04_get_est_prop","Tran_prop_long.Rdata"), verbose = TRUE)
prop_long_tran <- prop_long

load(here("processed-data", "13_PEC_deconvolution", "04_get_est_prop","PEC_prop_long.Rdata"), verbose = TRUE)
prop_long <- bind_rows(prop_long, prop_long_tran) |>
  bind_rows(prop_long_paired) |>
  mutate(factor(input, levels = c("paired", "Tran", "PEC")))

prop_long |> filter(is.na(RNAscope_prop)) |> count(input, method)
# input method     n
# <chr> <chr>  <int>
# 1 PEC   Bisque   434
# 2 PEC   hspe     434
# 3 Tran  Bisque   434
# 4 Tran  hspe     434
# 5 paired Bisque   434
# 6 paired hspe     434

#### correlation ####
(cor_check <- prop_long |>
   filter(!is.na(RNAscope_prop)) |>
   group_by(input, method) |>
   summarize(cor = cor(RNAscope_prop, prop),
             rmse = Metrics::rmse(RNAscope_prop, prop))  |>
   mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
   arrange(cor))

# input method   cor  rmse cor_anno               
# <chr> <chr>  <dbl> <dbl> <chr>                  
# 1 Tran   Bisque 0.102 0.196 "cor:0.102\nrmse:0.196"
# 2 PEC    hspe   0.296 0.120 "cor:0.296\nrmse:0.120"
# 3 Tran   hspe   0.503 0.143 "cor:0.503\nrmse:0.143"
# 4 paired Bisque 0.508 0.148 "cor:0.508\nrmse:0.148"
# 5 paired hspe   0.513 0.151 "cor:0.513\nrmse:0.151"
# 6 PEC    Bisque 0.534 0.138 "cor:0.534\nrmse:0.138"

(cor_check_library <- prop_long |>
    filter(!is.na(RNAscope_prop)) |>
    group_by(method, input, library_type, rna_extract, library_combo) |>
    summarize(cor = cor(RNAscope_prop, prop),
              rmse = Metrics::rmse(RNAscope_prop, prop)) |>
    arrange(-cor) |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))
)
# method input  library_type rna_extract library_combo        cor  rmse cor_anno               
# <chr>  <chr>  <chr>        <chr>       <chr>              <dbl> <dbl> <chr>                  
#   1 Bisque PEC    polyA        Cyto        polyA_Cyto         0.675 0.109 "cor:0.675\nrmse:0.109"
# 2 Bisque paired polyA        Cyto        polyA_Cyto         0.672 0.131 "cor:0.672\nrmse:0.131"
# 3 hspe   Tran   polyA        Cyto        polyA_Cyto         0.607 0.135 "cor:0.607\nrmse:0.135"
# 4 hspe   paired polyA        Cyto        polyA_Cyto         0.598 0.144 "cor:0.598\nrmse:0.144"
# 5 Bisque paired polyA        Nuc         polyA_Nuc          0.560 0.136 "cor:0.560\nrmse:0.136"
# 6 Bisque PEC    polyA        Nuc         polyA_Nuc          0.552 0.125 "cor:0.552\nrmse:0.125"
# 7 Bisque PEC    RiboZeroGold Cyto        RiboZeroGold_Cyto  0.524 0.141 "cor:0.524\nrmse:0.141"
# 8 hspe   paired RiboZeroGold Total       RiboZeroGold_Total 0.522 0.152 "cor:0.522\nrmse:0.152"
# 9 Bisque PEC    RiboZeroGold Nuc         RiboZeroGold_Nuc   0.519 0.163 "cor:0.519\nrmse:0.163"
# 10 Bisque PEC    RiboZeroGold Total       RiboZeroGold_Total 0.511 0.145 "cor:0.511\nrmse:0.145"

## cor check line/rank plot
cor_rmse_line <- cor_check_library |>
  mutate(method_in = paste(method, input)) |>
  ggplot(aes(x = library_combo, y = cor, color = method, group = method_in)) +
  geom_point(aes(size = rmse), alpha = .7) +
  geom_line(aes(linetype = input)) +
  scale_size(range = c(1,8)) +
  theme_bw() +
  scale_linetype_manual(values=c(paired = "solid", Tran = "dotted", PEC = "longdash")) +
  labs(x = "Library Type + RNA Extraction") +
  scale_color_manual(values = method_colors)

ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_input.png"), width = 10, height = 4)
ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_input.pdf"), width = 10, height = 4)

cor_rmse_scater_method <- cor_check_library |>
  ggplot(aes(x= cor, y = rmse, color = input, shape = library_combo)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~method, nrow = 1) +
  scale_shape_manual(values = library_combo_shapes2) 

ggsave(cor_rmse_scater_method, filename = here(plot_dir, "cor_rmse_scater_method_input.png"), width = 7, height = 4)
ggsave(cor_rmse_scater_method, filename = here(plot_dir, "cor_rmse_scater_method_input.pdf"), width = 7, height = 4)

## bar plot
prop_bar_SAMPLE_facet <- prop_long |> 
  mutate(Sample = gsub("_","\n", Sample)) |>
  ggplot(aes(x = library_combo, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_grid(method*input~Sample) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", x = "Library Type & RNA Extraction Prep", fill = "Cell Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet_input.png"), width = 12, height = 9)
ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet_input.pdf"), width = 12, height = 9)

## Scatter plots
## all points
est_prop_v_RNAscope_scatter <- prop_long |>
  filter(!is.na(RNAscope_prop),
         input %in% c("PEC", "Tran")) |>
  ggplot() +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check |> 
              filter(input %in% c("PEC", "Tran")) , 
            aes(label = cor_anno,x = .5, y = .75),
            vjust = "inward", hjust = "inward") +
  scale_shape_manual(values = library_combo_shapes2) +
  scale_color_manual(values = cell_type_colors_broad) +
  facet_grid(input~method) +
  # facet_wrap(input~method, ncol = 1) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_input.png"), width = 5, height = 5)
ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_input.pdf"), width = 5, height = 5)

#### ggpair plots ####
input_datasets <- unique(prop_long$input)

input_gg_prop <- map(unique(prop_long$method), function(.x){
  prop_wide <- prop_long |>
    filter(method == .x) |>
    select(SAMPLE_ID, library_combo, cell_type, input, RNAscope = RNAscope_prop, `snRNA-seq` = snRNA_prop, prop) |>
    pivot_wider(names_from = "input", values_from = "prop")
  
  # return(prop_wide)
  message(.x)

  gg_prop <- ggpairs(prop_wide, columns = c("RNAscope", "snRNA-seq", input_datasets), aes(color = cell_type)) +
    scale_color_manual(values = cell_type_colors_broad) +
    scale_fill_manual(values = cell_type_colors_broad) +
    theme_bw() +
    labs(title = .x)

  ggsave(gg_prop, filename = here(plot_dir, paste0("ggpairs_prop_input_",.x,".png")), height = 10, width = 10)
  return(gg_prop)
}
)

