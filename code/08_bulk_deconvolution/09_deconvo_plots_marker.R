
library("tidyverse")
library("sessioninfo")
library("DeconvoBuddies")
library("here")
library("Metrics")
library("survival")
library("viridis")
library("GGally")
library("patchwork")
library("ggrepel")
library("broom")

#### prep dirs & plot info ####
plot_dir <- here("plots", "08_bulk_deconvolution", "09_deconvo_plots_marker")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "08_bulk_deconvolution", "09_deconvo_plots_marker")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

## load colors
load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo
# cell_type_colors_broad
load(here("processed-data","00_data_prep","method_colors.Rdata"), verbose = TRUE)
# method_colors
load(here("processed-data","00_data_prep","library_combo_shapes.Rdata"), verbose = TRUE)
# library_combo_shapes
# library_combo_shapes2

#### load data ####
load(here("processed-data", "08_bulk_deconvolution", "03_get_est_prop","prop_long.Rdata"), verbose = TRUE)

prop_long |>  count(method, marker)|> count(marker)

prop_long |>  filter(!is.na(RNAscope_prop)) |> count(method, marker)
prop_long |> count(cell_type)

#### get n marker genes ####
marker_list <- list.files(here("processed-data","08_bulk_deconvolution"), pattern = "markers_.*.txt", full.names = TRUE)
names(marker_list) <- gsub("markers_|.txt", "", basename(marker_list))

## HVG files
hvg_files <- here("processed-data", "06_marker_genes", sprintf("09_HVGs/HVG%d0.txt", seq(1,10)))
names(hvg_files) <- sprintf("HVG%d0", seq(1,10))

marker_list <- c(marker_list, hvg_files)

n_markers <- c(map_int(marker_list, ~length(scan(.x, what="", sep="\n"))), FULL = 17804)

n_marker_tb <- tibble(marker = names(n_markers), n_markers =n_markers) |>
  arrange(n_markers)
# marker          n_markers
# <chr>               <dbl>
# 1 1vALL_top25           145
# 2 MeanRatio_top25       151
# 3 MeanRatio_MAD3        520
# 4 MeanRatio_over2       557
# 5 FULL                17804

#### compare to RNAscope ####

concordance <- prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  group_by(method, marker) |> 
  group_map(~survival::concordance(RNAscope_prop~prop, data = .x)$concordance) |>
  unlist()

## Overall correlation

(cor_check <- prop_long |>
    filter(!is.na(RNAscope_prop)) |>
    group_by(method, marker) |>
    summarize(cor = cor(RNAscope_prop, prop), #pearson
              cor_spearman = cor(RNAscope_prop, prop, method = "spearman"),
              rmse = Metrics::rmse(RNAscope_prop, prop),
    )  |>
    add_column(concordance = concordance) |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
    arrange(-cor) |>
    left_join(n_marker_tb)  |>
    mutate(marker = factor(marker, 
                           levels = c("FULL", "1vALL_top25", "MeanRatio_MAD3", "MeanRatio_over2", "MeanRatio_top25",
                                      sprintf("HVG%d0", seq(1,10)))))
)
# method     marker            cor cor_spearman  rmse concordance cor_anno                n_markers
# <chr>      <fct>           <dbl>        <dbl> <dbl>       <dbl> <chr>                       <dbl>
#   1 CIBERSORTx HVG10           0.598        0.535 0.200       0.690 "cor:0.598\nrmse:0.200"       811
# 2 hspe       MeanRatio_over2 0.596        0.518 0.215       0.687 "cor:0.596\nrmse:0.215"       557
# 3 CIBERSORTx HVG30           0.590        0.488 0.191       0.667 "cor:0.590\nrmse:0.191"      2434
# 4 CIBERSORTx HVG20           0.590        0.465 0.187       0.657 "cor:0.590\nrmse:0.187"      1623
# 5 hspe       1vALL_top25     0.586        0.458 0.206       0.665 "cor:0.586\nrmse:0.206"       145
# 6 hspe       MeanRatio_MAD3  0.585        0.402 0.212       0.640 "cor:0.585\nrmse:0.212"       520
# 7 CIBERSORTx HVG40           0.577        0.498 0.167       0.667 "cor:0.577\nrmse:0.167"      3245
# 8 CIBERSORTx HVG80           0.576        0.469 0.184       0.662 "cor:0.576\nrmse:0.184"      6490
# 9 CIBERSORTx HVG100          0.575        0.473 0.191       0.658 "cor:0.575\nrmse:0.191"      8113
# 10 CIBERSORTx HVG50           0.570        0.474 0.178       0.660 "cor:0.570\nrmse:0.178"      4056

## plot spearman's vs. pearson's correlation

cor_check |> arrange(-cor_spearman)

cor_method_scatter <- cor_check |>
  filter(!grepl("HVG", marker)) |>
  ggplot(aes(cor, cor_spearman, color = method, shape = marker)) +
  geom_point() +
  scale_color_manual(values = method_colors) +
  theme_bw() +
  coord_equal() +
  geom_abline()

ggsave(cor_method_scatter, filename = here(plot_dir, "Correlation_method_scatter.png"), height = 5)

cor_method_scatter_HVG <- cor_check |>
  filter(!grepl("HVG", marker)) |>
  ggplot(aes(cor, cor_spearman, color = method)) +
  geom_point() +
  geom_text_repel(aes(label = marker), size = 2) +
  scale_color_manual(values = method_colors) +
  theme_bw() +
  coord_equal() +
  geom_abline()

ggsave(cor_method_scatter_HVG, filename = here(plot_dir, "Correlation_method_scatter_HVG.png"), height = 5)

## concordance vs. correlation
cor_concord_scatter <- cor_check |>
  filter(!grepl("HVG", marker)) |>
  ggplot(aes(cor, concordance, color = method, shape = marker)) +
  geom_point() +
  scale_color_manual(values = method_colors) +
  theme_bw() +
  theme(legend.position = "None")

cor_spear_concord_scatter <- cor_check |>
  filter(!grepl("HVG", marker)) |>
  ggplot(aes(cor_spearman, concordance, color = method, shape = marker)) +
  geom_point() +
  scale_color_manual(values = method_colors) +
  theme_bw()

ggsave(cor_concord_scatter + cor_spear_concord_scatter, filename = here(plot_dir, "Correlation_Concordance_scatter.png"), height = 5, width = 9)

## factor method by overall cor from MeanRatio_top25
method_levels <- cor_check |> filter(marker == "MeanRatio_top25") |> arrange(cor) |> pull(method)

cor_check$method <- factor(cor_check$method, levels = method_levels)
prop_long$method <- factor(prop_long$method, levels = method_levels)

## plot overall cor vs. rmse
cor_rmse_scater_overall <- cor_check |>
  ggplot(aes(x= cor, y = rmse, shape= method, color = marker)) +
  geom_point(size = 2) +
  theme_bw() 

ggsave(cor_rmse_scater_overall, filename = here(plot_dir, "cor_rmse_scater_overall.png"), width = 9, height =7)
ggsave(cor_rmse_scater_overall, filename = here(plot_dir, "cor_rmse_scater.pdf"), width = 10, height = 3.4)

## text version
# cor_check |>
#   mutate(method_short = ifelse(substring(method, 1,1) == "B", 
#                                substring(method, 1,2), 
#                                substring(method, 1,1))) |>
#   ggplot(aes(x= cor, y = rmse, color = marker)) +
#   geom_text(aes(label = method_short)) +
#   theme_bw() 


#### Correlation by library type ####
concordance_library <- prop_long |>
  filter(!is.na(RNAscope_prop)) |>
  group_by(method, marker, library_combo) |> 
  group_map(~survival::concordance(RNAscope_prop~prop, data = .x)$concordance) |>
  unlist()

cor_check_library <- prop_long |>
    filter(!is.na(RNAscope_prop)) |>
    group_by(marker, method, library_combo) |>
    summarize(cor = cor(RNAscope_prop, prop),
              cor_spearman = cor(RNAscope_prop, prop, method = "spearman"),
              rmse = Metrics::rmse(RNAscope_prop, prop)) |>
    add_column(concordance = concordance_library) |>
    arrange(-cor) |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)),
           cor_spear_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor_spearman,3), round(rmse,3)),
           ) |>
    left_join(n_marker_tb) |>
    mutate(marker = factor(marker, 
                           levels = c("FULL","1vALL_top25","MeanRatio_MAD3", "MeanRatio_over2", "MeanRatio_top25",
                                      sprintf("HVG%d0", seq(1,10)))))

cor_check_library |> select(1:6)
# marker          method     library_combo   cor cor_spearman  rmse
# <fct>           <fct>      <chr>         <dbl>        <dbl> <dbl>
# 1 HVG10           CIBERSORTx polyA_Cyto    0.700        0.385 0.233
# 2 MeanRatio_over2 hspe       polyA_Cyto    0.690        0.469 0.215
# 3 HVG20           CIBERSORTx polyA_Cyto    0.684        0.388 0.214
# 4 MeanRatio_over2 CIBERSORTx polyA_Cyto    0.684        0.402 0.263
# 5 MeanRatio_top25 Bisque     polyA_Cyto    0.683        0.568 0.123
# 6 HVG30           CIBERSORTx polyA_Cyto    0.682        0.397 0.221
# 7 1vALL_top25     hspe       polyA_Cyto    0.670        0.375 0.205
# 8 MeanRatio_MAD3  hspe       polyA_Cyto    0.669        0.300 0.213
# 9 HVG10           Bisque     polyA_Cyto    0.665        0.442 0.133
# 10 MeanRatio_MAD3  CIBERSORTx polyA_Cyto    0.655        0.173 0.246

## write correlation values
write_csv(cor_check, file = here(data_dir, "est_prop_correlation.csv"))
write_csv(cor_check_library, file = here(data_dir, "est_prop_correlation_library.csv"))

## cor vs rmse ##
cor_rmse_scater <- cor_check_library |>
  filter(!grepl("HVG", marker)) |>
  ggplot(aes(x= cor, y = rmse, color= method, shape = library_combo)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~marker, nrow = 1) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = library_combo_shapes2) 

ggsave(cor_rmse_scater, filename = here(plot_dir, "cor_rmse_scater.png"), width = 10, height = 3.4)
ggsave(cor_rmse_scater, filename = here(plot_dir, "cor_rmse_scater.pdf"), width = 10, height = 3.4)

cor_rmse_scater_HVG <- cor_check_library |>
  filter(grepl("HVG", marker)) |>
  ggplot(aes(x= cor, y = rmse, color= method, shape = library_combo)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~marker, nrow = 2) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = library_combo_shapes2) 

ggsave(cor_rmse_scater_HVG, filename = here(plot_dir, "cor_rmse_scater_HVG.png"), width = 10, height = 6)

## note several values are identical for BayesPrism - maybe use geom_jitter?
cor_check_library |>
  filter(method == "BayesPrism") |>
  arrange(cor)

## factor by method
cor_rmse_scater_method <- cor_check_library |>
  filter(!grepl("HVG", marker)) |>
  ggplot(aes(x= cor, y = rmse, color = marker, shape = library_combo)) +
  geom_point() +
  # geom_jitter(width = 0.25) +
  theme_bw() +
  facet_wrap(~method, nrow = 1) +
  scale_shape_manual(values = library_combo_shapes2) 
  
ggsave(cor_rmse_scater_method, filename = here(plot_dir, "cor_rmse_scater_method.png"), width = 11, height = 4)
ggsave(cor_rmse_scater_method, filename = here(plot_dir, "cor_rmse_scater_method.pdf"), width = 11, height = 4)

cor_rmse_scater_method_HVG <- cor_check_library |>
  ggplot(aes(x= cor, y = rmse, color = marker, shape = library_combo)) +
  geom_point() +
  # geom_jitter(width = 0.25) +
  theme_bw() +
  facet_wrap(~method, nrow = 1) +
  scale_shape_manual(values = library_combo_shapes2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
ggsave(cor_rmse_scater_method_HVG, filename = here(plot_dir, "cor_rmse_scater_method_HVG.png"), width = 12.5, height = 6)
ggsave(cor_rmse_scater_method_HVG, filename = here(plot_dir, "cor_rmse_scater_method_HVG.pdf"), width = 12.5, height = 6)

## cor check line/rank plot
cor_rmse_line <- cor_check_library |>
  filter(marker %in% c('FULL', "1vALL_top25", "MeanRatio_top25")) |>
  mutate(group = paste(method, marker)) |>
  ggplot(aes(x = library_combo, y = cor, color= method)) +
  geom_point(aes(size = 1/rmse), alpha = .7) +
  geom_line(aes(group = group, linetype = marker)) +
  scale_size(range = c(1,8)) +
  theme_bw() +
  scale_linetype_manual(values=c(FULL = "solid", `1vALL_top25` = "dotted", `MeanRatio_top25` = "longdash")) +
  labs(x = "Library Type + RNA Extraction")

ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_markers.png"), width = 10, height = 5)
ggsave(cor_rmse_line, filename = here(plot_dir, "cor_rmse_line_markers.pdf"), width = 10, height = 5)

cor_rmse_line_facet <- cor_check_library |>
  ggplot(aes(x = library_combo, y = cor, color= method)) +
  geom_point(aes(size = rmse), alpha = .7) +
  geom_line(aes(group = method)) +
  scale_size(range = c(1,8)) +
  theme_bw() +
  scale_color_manual(values = method_colors) +
  labs(x = "Library Type + RNA Extraction") +
  facet_wrap(~marker, ncol = 2)

ggsave(cor_rmse_line_facet, filename = here(plot_dir, "cor_rmse_line_markers_facet.png"), width = 15, height = 15)


cor_rmse_line_facet2 <- cor_check_library |>
  filter(marker != "MeanRatio_top25") |>
  filter(!grepl("HVG", marker)) |>
  ggplot(aes(x = library_combo, y = cor, color= method)) +
  geom_point(aes(size = cor), alpha = .7) +
  geom_line(aes(group = method)) +
  scale_size(range = c(1,8)) +
  theme_bw() +
  scale_color_manual(values = method_colors) +
  # scale_linetype_manual(values=c(FULL = "solid", `1vALL_top25` = "dotted", `MeanRatio_top25` = "longdash")) +
  labs(x = "Library Type + RNA Extraction") +
  facet_wrap(~marker, ncol =1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(cor_rmse_line_facet2, filename = here(plot_dir, "cor_rmse_line_markers_facet2.png"), width = 10, height = 9)
ggsave(cor_rmse_line_facet2, filename = here(plot_dir, "cor_rmse_line_markers_facet2.pdf"), width = 10, height = 9)

rmse_line_facet <- cor_check_library |>
  filter(!grepl("HVG", marker)) |>
  ggplot(aes(x = library_combo, y = rmse, color= method)) +
  geom_point(aes(size = cor), alpha = .7) +
  geom_line(aes(group = method)) +
  scale_size(range = c(1,8)) +
  theme_bw() +
  # scale_linetype_manual(values=c(FULL = "solid", `1vALL_top25` = "dotted", `MeanRatio_top25` = "longdash")) +
  scale_color_manual(values = method_colors) +
  labs(x = "Library Type + RNA Extraction") +
  facet_wrap(~marker) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(rmse_line_facet, filename = here(plot_dir, "rmse_line_markers_facet.png"), width = 10, height = 7)
ggsave(rmse_line_facet, filename = here(plot_dir, "rmse_line_markers_facet.pdf"), width = 10, height = 7)

## spearman line rank
cor_spear_rmse_line_top25 <- cor_check_library |>
  filter(marker == "MeanRatio_top25") |>
  filter(!grepl("HVG", marker)) |>
  mutate(group = paste(method, marker)) |>
  ggplot(aes(x = library_combo, y = cor_spearman, color= method)) +
  geom_point(aes(size = rmse), alpha = .7) +
  geom_line(aes(group = group)) +
  # geom_line(aes(group = group)) +
  scale_size(range = c(1,8)) +
  theme_bw() +
  scale_color_manual(values = method_colors) +
  labs(x = "Library Type + RNA Extraction")

ggsave(cor_spear_rmse_line_top25, filename = here(plot_dir, "cor_spear_rmse_line_top25.png"), width = 10, height = 3.4)
ggsave(cor_spear_rmse_line_top25, filename = here(plot_dir, "cor_spear_rmse_line_top25.pdf"), width = 10, height = 3.4)


#### corelation vs. n markers ####
# stats about terms
(cor_check_tidy <- cor_check |>
  do(tidy(lm(cor~n_markers,.))))
# method     term            estimate   std.error statistic  p.value
# <chr>      <chr>              <dbl>       <dbl>     <dbl>    <dbl>
#   1 BayesPrism (Intercept) -0.0774      0.0291         -2.66  1.96e- 2
# 2 BayesPrism n_markers    0.0000203   0.00000472      4.30  8.65e- 4 ***
# 3 Bisque     (Intercept)  0.508       0.00464       109.    1.16e-20
# 4 Bisque     n_markers    0.000000572 0.000000753     0.760 4.61e- 1
# 5 CIBERSORTx (Intercept)  0.557       0.0111         50.3   2.77e-16
# 6 CIBERSORTx n_markers    0.000000981 0.00000179      0.547 5.94e- 1
# 7 DWLS       (Intercept)  0.150       0.0360          4.17  1.09e- 3
# 8 DWLS       n_markers   -0.0000101   0.00000583     -1.73  1.07e- 1
# 9 MuSiC      (Intercept)  0.0462      0.0726          0.636 5.36e- 1
# 10 MuSiC      n_markers   -0.0000232   0.0000118      -1.97  7.02e- 2
# 11 hspe       (Intercept)  0.345       0.0836          4.12  1.20e- 3
# 12 hspe       n_markers   -0.0000284   0.0000136      -2.09  5.67e- 2



# summary statistics for the entire regression, such as R^2 and the F-statistic.
(cor_check_summary <- cor_check |>
  summarise(glance(lm(cor~n_markers))))
# method     r.squared adj.r.squared  sigma statistic  p.value    df logLik    AIC     BIC deviance df.residual  nobs
# <chr>          <dbl>         <dbl>  <dbl>     <dbl>    <dbl> <dbl>  <dbl>  <dbl>   <dbl>    <dbl>       <int> <int>
# 1 BayesPrism    0.587         0.555  0.0817    18.5   0.000865     1  17.4  -28.7  -26.6    0.0867           13    15 ***
# 2 Bisque        0.0425       -0.0311 0.0130     0.577 0.461        1  44.9  -83.8  -81.7    0.00221          13    15
# 3 CIBERSORTx    0.0225       -0.0527 0.0311     0.299 0.594        1  31.9  -57.7  -55.6    0.0125           13    15
# 4 DWLS          0.187         0.124  0.101      2.99  0.107        1  14.2  -22.4  -20.2    0.133            13    15
# 5 MuSiC         0.230         0.171  0.204      3.89  0.0702       1   3.65  -1.29   0.831  0.540            13    15
# 6 hspe          0.252         0.194  0.235      4.37  0.0567       1   1.53   2.94   5.06   0.716            13    15

anno_y <- cor_check_tidy |> 
  select(method, term, estimate) |>
  mutate(term = gsub("\\(|\\)", "", term)) |>
  pivot_wider(values_from = "estimate", names_from = "term") |>
  mutate(y_anno = (15000*n_markers) + Intercept) |>
  left_join(cor_check_summary |> select(method, adj.r.squared, p.value)) |>
  mutate(p_anno = ifelse(p.value < 0.05, "*", ""),
         r2_anno = paste0("R2=", round(adj.r.squared, 2)), p_anno)

cor_n_marker_scatter <- 
  ggplot(data = cor_check, aes(x = n_markers, y = cor, color= method)) +
  geom_smooth(method = "lm", alpha=0.2, linewidth = 0.5) +
  geom_point(size = 2) +
  scale_color_manual(values = method_colors) +
  theme_bw() +
  geom_text(data = anno_y, aes(x = 15000, y = y_anno + 0.05, label = r2_anno))

ggsave(cor_n_marker_scatter, filename = here(plot_dir, "cor_n_marker_scatter.png"), width = 10)
ggsave(cor_n_marker_scatter, filename = here(plot_dir, "cor_n_marker_scatter.pdf"), width = 10)


rmse_check_tidy <- cor_check |>
  do(tidy(lm(rmse~n_markers,.)))

rmse_check_summary <- cor_check |>
  summarise(glance(lm(rmse~n_markers)))

rmse_anno_y <- rmse_check_tidy |> 
  select(method, term, estimate) |>
  mutate(term = gsub("\\(|\\)", "", term)) |>
  pivot_wider(values_from = "estimate", names_from = "term") |>
  mutate(y_anno = (15000*n_markers) + Intercept) |>
  left_join(rmse_check_summary |> select(method, adj.r.squared)) |>
  mutate(r2_anno = paste0("R2=", round(adj.r.squared, 2))) 


rmse_n_marker_scatter <- 
  ggplot(data = cor_check, aes(x = n_markers, y = rmse, color= method)) +
  geom_smooth(method = "lm", alpha=0.2, linewidth = 0.5) +
  geom_point(size = 2) +
  scale_color_manual(values = method_colors) +
  theme_bw() +
  geom_text(data = rmse_anno_y, aes(x = 15000, y = y_anno + 0.005, label = r2_anno))

ggsave(rmse_n_marker_scatter, filename = here(plot_dir, "rmse_n_marker_scatter.png"), width = 10)
ggsave(rmse_n_marker_scatter, filename = here(plot_dir, "rmse_n_marker_scatter.pdf"), width = 10)


cor_n_marker_scatter_r2 <- ggplot(data = cor_check, aes(x = n_markers, y = cor, color = method)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  ggpubr::stat_regline_equation(label.x = with(cor_check,tapply(n_markers,method,quantile,.6)),
                        label.y = with(cor_check,tapply(cor,method,max)- 0.2),
                        aes(label = ..adj.rr.label..),
                        show.legend = FALSE) +
  scale_color_manual(values = method_colors) +
  theme_bw()

ggsave(cor_n_marker_scatter_r2, filename = here(plot_dir, "cor_n_marker_scatter_r2.png"), width = 10)

cor_n_marker_scatter_log <- 
  ggplot() +
  geom_vline(data = n_marker_tb, aes(xintercept = n_markers), color = "skyblue") +
  geom_point(data = cor_check, aes(x = n_markers, y = cor, color= method)) +
  scale_color_manual(values = method_colors) +
  theme_bw() +
  coord_trans(x = "log10") +
  labs(x = "n marker genes (log10 axis)")

ggsave(cor_n_marker_scatter_log, filename = here(plot_dir, "cor_n_marker_scatter_log.png"), width = 10)
ggsave(cor_n_marker_scatter_log, filename = here(plot_dir, "cor_n_marker_scatter_log.pdf"), width = 10)

cor_check |>
  filter(method == "hspe") |>
  mutate(gene_set = paste0(marker, " (", n_markers, ")"),
         gene_set = fct_reorder(gene_set, n_markers, .na_rm = FALSE)) |>
    arrange(gene_set)

cor_marker_scatter <- cor_check |>
  mutate(gene_set = fct_reorder(paste0(marker, " (", n_markers, ")"), n_markers)) |>
  ggplot(aes(x = gene_set, y = cor, color= method)) +
  geom_point(size = 2) +
  geom_line(aes(group = method)) +
  # geom_smooth() + # not working because marker is a factor - fix
  scale_color_manual(values = method_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(cor_marker_scatter, filename = here(plot_dir, "cor_marker_scatter.png"), width = 10)
ggsave(cor_marker_scatter, filename = here(plot_dir, "cor_marker_scatter.pdf"), width = 10)


rmse_marker_scatter <- cor_check |>
  mutate(gene_set = fct_reorder(paste0(marker, " (", n_markers, ")"), n_markers)) |>
  ggplot(aes(x = gene_set, y = rmse, color= method)) +
  geom_point(size = 2) +
  geom_line(aes(group = method)) +
  # geom_smooth() + # not working because marker is a factor - fix
  scale_color_manual(values = method_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(rmse_marker_scatter, filename = here(plot_dir, "rmse_marker_scatter.png"), width = 10)


cor_n_marker_library_scatter <- cor_check_library |>
  # filter(marker != "FULL") |>
  ggplot(aes(x = n_markers, y = cor, color= method)) +
  geom_point(aes(shape = library_combo)) +
  # geom_point(aes(size = rmse), alpha = .7) +
  # geom_line(aes(group = method)) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = library_combo_shapes2) +
  theme_bw() +
  coord_trans(x = "log10")

ggsave(cor_n_marker_scatter, filename = here(plot_dir, "cor_n_marker_scatter.png"))


#### proportion data ####
prop_bar_SAMPLE_facet <- prop_long_opc |> 
  filter(marker %in% c("1vALL_top25", "FULL")) |>
  mutate(Sample = gsub("_","\n", Sample),
         marker = gsub("_","\n", marker)) |>
  ggplot(aes(x = library_combo, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_grid(method*marker~Sample) +
  scale_fill_manual(values = cell_type_colors_broad) +
  labs(y = "Cell Type Proportion", x = "Library Type & RNA Extraction Prep", fill = "Cell Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")

ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet.png"), width = 12, height = 12)
ggsave(prop_bar_SAMPLE_facet, filename = here(plot_dir, "Bulk_prop_SAMPLE_facet.pdf"), width = 12, height = 9)


## Scatter plots
## all points
est_prop_v_RNAscope_scatter <- prop_long |>
  filter(!is.na(RNAscope_prop),
         marker %in% c("1vALL_top25", "FULL")) |>
  ggplot() +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check |> 
              filter(marker %in% c("1vALL_top25", "FULL")) , 
            aes(label = cor_anno,x = .5, y = 1),
            vjust = "inward", hjust = "inward") +
  scale_shape_manual(values = library_combo_shapes2) +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_grid(marker~method) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter.png"), width = 10, height = 5)
ggsave(est_prop_v_RNAscope_scatter, filename = here(plot_dir, "est_prop_v_RNAscope_scatter.pdf"), width = 10, height = 5)

## Other mean ratio approaches 
est_prop_v_RNAscope_scatter_MR <- prop_long |>
  filter(!is.na(RNAscope_prop),
         marker %in% c("MeanRatio_MAD3", "MeanRatio_over2")) |>
  ggplot() +
  geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check |> 
              filter(marker %in% c("MeanRatio_MAD3", "MeanRatio_over2")) , 
            aes(label = cor_anno,x = .5, y = 1),
            vjust = "inward", hjust = "inward") +
  scale_shape_manual(values = library_combo_shapes2) +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_grid(marker~method) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "RNAscope Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_RNAscope_scatter_MR, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_MR.png"), width = 10, height = 5)
ggsave(est_prop_v_RNAscope_scatter_MR, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_MR.pdf"), width = 10, height = 5)

# est_prop_v_RNAscope_scatter_top25_library <- prop_long |>
#   filter(!is.na(RNAscope_prop)) |>
#   ggplot() +
#   scale_color_manual(values = cell_type_colors_halo) +
#   geom_point(aes(x = RNAscope_prop, y = prop, color = cell_type, shape = library_combo)) +
#   geom_text(data = cor_check_library,
#             aes(label = cor_anno,x = .5, y = 1),
#             vjust = "inward", hjust = "inward") +
#   scale_shape_manual(values = library_combo_shapes2) +
#   scale_color_manual(values = cell_type_colors_halo) +
#   facet_grid(library~method) +
#   geom_abline() +
#   # coord_equal() +
#   theme_bw() +
#   labs( x = "RNAscope Proportion", y = "Estimated Proportion")
# 
# ggsave(est_prop_v_RNAscope_scatter_top25_library, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_top25_library.png"), width = 10, height = 9)
# ggsave(est_prop_v_RNAscope_scatter_top25_library, filename = here(plot_dir, "est_prop_v_RNAscope_scatter_top25_library.pdf"), width = 10, height = 9)

#### ggpair plots ####
sn_prop <- read_csv(here("processed-data", "03_HALO", "08_explore_proportions","snRNA_cell_type_proportions.csv")) |>
  select(Sample, cell_type, prop_sn) |>
  mutate(cell_type = gsub("Endo", "EndoMural", cell_type))

marker_sets <- unique(prop_long$marker)

method_gg_prop <- map(levels(prop_long$method), function(.x){
  prop_wide <- prop_long |>
    filter(method == .x) |>
    select(SAMPLE_ID, library_combo, cell_type, marker, RNAscope = RNAscope_prop, `snRNA-seq` = snRNA_prop, prop) |>
    pivot_wider(names_from = "marker", values_from = "prop")
  
  message(.x)
  
  gg_prop <- ggpairs(prop_wide, columns = c("RNAscope", "snRNA-seq", marker_sets), aes(color = cell_type)) +
    scale_color_manual(values = cell_type_colors_halo) +
    scale_fill_manual(values = cell_type_colors_halo) +
    theme_bw() +
    labs(title = .x)
  
  ggsave(gg_prop, filename = here(plot_dir, paste0("ggpairs_prop_",.x,".png")), height = 10, width = 10)
  return(gg_prop)
  }
)

#### est_prop vs. snRNA props? ####

(cor_check_sn <- prop_long_opc |>
    filter(!is.na(snRNA_prop)) |>
    group_by(method, marker) |>
    summarize(cor = cor(snRNA_prop, prop),
              rmse = Metrics::rmse(snRNA_prop, prop))  |>
    mutate(cor_anno = sprintf("cor:%.3f\nrmse:%.3f", round(cor,3), round(rmse,3)))|>
    arrange(-cor))

# method     marker            cor  rmse cor_anno               
# <chr>      <chr>           <dbl> <dbl> <chr>                  
#   1 Bisque     MeanRatio_top25 0.757 0.118 "cor:0.757\nrmse:0.118"
# 2 hspe       MeanRatio_top25 0.710 0.128 "cor:0.710\nrmse:0.128"
# 3 Bisque     FULL            0.699 0.131 "cor:0.699\nrmse:0.131"
# 4 hspe       MeanRatio_over2 0.673 0.176 "cor:0.673\nrmse:0.176"
# 5 Bisque     MeanRatio_MAD3  0.672 0.138 "cor:0.672\nrmse:0.138"
# 6 Bisque     MeanRatio_over2 0.671 0.138 "cor:0.671\nrmse:0.138"
# 7 Bisque     1vALL_top25     0.668 0.138 "cor:0.668\nrmse:0.138"
# 8 hspe       MeanRatio_MAD3  0.667 0.174 "cor:0.667\nrmse:0.174"
# 9 hspe       1vALL_top25     0.663 0.172 "cor:0.663\nrmse:0.172"
# 10 CIBERSORTx MeanRatio_MAD3  0.580 0.180 "cor:0.580\nrmse:0.180"

## factor method by overall cor

method_levels_sn <- cor_check_sn |> filter(marker == "MeanRatio_top25") |> arrange(cor) |> pull(method)

cor_check_sn$method <- factor(cor_check_sn$method, levels = method_levels_sn)
prop_long_opc$method <- factor(prop_long_opc$method, levels = method_levels_sn)

est_prop_v_sn_scatter <- prop_long_opc |>
  filter(!is.na(snRNA_prop)) |>
  ggplot() +
  geom_point(aes(x = snRNA_prop, y = prop, color = cell_type, shape = library_combo)) +
  geom_text(data = cor_check_sn, 
            aes(label = cor_anno,x = .8, y = 1),
            vjust = "inward", hjust = "inward") +
  facet_grid(marker~method) +
  scale_color_manual(values = cell_type_colors_broad) +
  scale_shape_manual(values = library_combo_shapes2) +
  geom_abline() +
  coord_equal() +
  theme_bw() +
  labs( x = "snRNA-seq Proportion", y = "Estimated Proportion")

ggsave(est_prop_v_sn_scatter, filename = here(plot_dir, "est_prop_v_sn_scatter_markers.png"), width = 10, height = 10)
ggsave(est_prop_v_sn_scatter, filename = here(plot_dir, "est_prop_v_sn_scatter_markers.pdf"), width = 10, height = 10)



# sgejobs::job_single('09_deconvo_plots_marker', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript 09_deconvo_plots_marker.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       macOS Sonoma 14.1
# system   x86_64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/New_York
# date     2024-02-02
# rstudio  2023.12.0+369 Ocean Storm (desktop)
# pandoc   3.1.1 @ /usr/local/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
# beachmat               2.18.0    2023-10-24 [1] Bioconductor
# Biobase                2.62.0    2023-10-24 [1] Bioconductor
# BiocGenerics           0.48.1    2023-11-01 [1] Bioconductor
# BiocNeighbors          1.20.2    2024-01-07 [1] Bioconductor 3.18 (R 4.3.2)
# BiocParallel           1.36.0    2023-10-24 [1] Bioconductor
# BiocSingular           1.18.0    2023-10-24 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# bluster                1.12.0    2023-10-24 [1] Bioconductor
# cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.0)
# cluster                2.1.6     2023-12-01 [1] CRAN (R 4.3.0)
# codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.2)
# colorout             * 1.3-0.1   2024-01-10 [1] Github (jalvesaq/colorout@deda341)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DeconvoBuddies       * 0.99.0    2023-05-12 [1] Github (LieberInstitute/DeconvoBuddies@9ce4a42)
# DelayedArray           0.28.0    2023-10-24 [1] Bioconductor
# DelayedMatrixStats     1.24.0    2023-10-24 [1] Bioconductor
# devtools             * 2.4.5     2022-10-11 [1] CRAN (R 4.3.0)
# digest                 0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
# dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.0)
# dqrng                  0.3.2     2023-11-29 [1] CRAN (R 4.3.0)
# edgeR                  4.0.6     2024-01-08 [1] Bioconductor 3.18 (R 4.3.2)
# ellipsis               0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.0)
# fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.3.0)
# fs                     1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb           1.38.5    2023-12-28 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData       1.2.11    2024-01-09 [1] Bioconductor
# GenomicRanges          1.54.1    2023-10-29 [1] Bioconductor
# GGally               * 2.2.0     2023-11-22 [1] CRAN (R 4.3.0)
# ggplot2              * 3.4.4     2023-10-12 [1] CRAN (R 4.3.0)
# ggstats                0.5.1     2023-11-21 [1] CRAN (R 4.3.0)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.0)
# gridExtra              2.3       2017-09-09 [1] CRAN (R 4.3.0)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
# htmltools              0.5.7     2023-11-03 [1] CRAN (R 4.3.0)
# htmlwidgets            1.6.4     2023-12-06 [1] CRAN (R 4.3.0)
# httpuv                 1.6.13    2023-12-06 [1] CRAN (R 4.3.0)
# igraph                 1.6.0     2023-12-11 [1] CRAN (R 4.3.0)
# IRanges                2.36.0    2023-10-24 [1] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.3.0)
# later                  1.3.2     2023-12-06 [1] CRAN (R 4.3.0)
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.0)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.0)
# limma                  3.58.1    2023-10-31 [1] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [1] CRAN (R 4.3.0)
# lubridate            * 1.9.3     2023-09-27 [1] CRAN (R 4.3.0)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-4     2023-11-30 [1] CRAN (R 4.3.0)
# MatrixGenerics         1.14.0    2023-10-24 [1] Bioconductor
# matrixStats            1.2.0     2023-12-11 [1] CRAN (R 4.3.0)
# memoise                2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
# metapod                1.10.1    2023-12-24 [1] Bioconductor 3.18 (R 4.3.2)
# Metrics                0.1.4     2018-07-09 [1] CRAN (R 4.3.0)
# mime                   0.12      2021-09-28 [1] CRAN (R 4.3.0)
# miniUI                 0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgbuild               1.4.3     2023-12-10 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# pkgload                1.3.3     2023-09-22 [1] CRAN (R 4.3.0)
# plyr                   1.8.9     2023-10-02 [1] CRAN (R 4.3.0)
# profvis                0.3.8     2023-05-02 [1] CRAN (R 4.3.0)
# promises               1.2.1     2023-08-10 [1] CRAN (R 4.3.0)
# purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# rafalib                1.0.0     2015-08-09 [1] CRAN (R 4.3.0)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.0)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.0)
# readr                * 2.1.4     2023-02-10 [1] CRAN (R 4.3.0)
# remotes                2.4.2.1   2023-07-18 [1] CRAN (R 4.3.0)
# rlang                  1.1.3     2024-01-10 [1] CRAN (R 4.3.0)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.0)
# rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# rsvd                   1.0.5     2021-04-16 [1] CRAN (R 4.3.0)
# S4Arrays               1.2.0     2023-10-24 [1] Bioconductor
# S4Vectors              0.40.2    2023-11-23 [1] Bioconductor
# ScaledMatrix           1.10.0    2023-10-24 [1] Bioconductor
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.0)
# scran                  1.30.0    2023-10-24 [1] Bioconductor
# scuttle                1.12.0    2023-10-24 [1] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# shiny                  1.8.0     2023-11-17 [1] CRAN (R 4.3.0)
# SingleCellExperiment   1.24.0    2023-10-24 [1] Bioconductor
# SparseArray            1.2.3     2023-12-25 [1] Bioconductor 3.18 (R 4.3.2)
# sparseMatrixStats      1.14.0    2023-10-24 [1] Bioconductor
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.0)
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.0)
# stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.3.0)
# SummarizedExperiment   1.32.0    2023-10-24 [1] Bioconductor
# tibble               * 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                * 1.3.0     2023-01-24 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.3.0)
# timechange             0.2.0     2023-01-11 [1] CRAN (R 4.3.0)
# tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.3.0)
# urlchecker             1.0.1     2021-11-30 [1] CRAN (R 4.3.0)
# usethis              * 2.2.2     2023-07-06 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.0)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.0)
# viridis              * 0.6.4     2023-07-22 [1] CRAN (R 4.3.0)
# viridisLite          * 0.4.2     2023-05-02 [1] CRAN (R 4.3.0)
# withr                  2.5.2     2023-10-30 [1] CRAN (R 4.3.0)
# xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# zlibbioc               1.48.0    2023-10-24 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library
# 
