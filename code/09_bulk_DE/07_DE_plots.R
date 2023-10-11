
library("SummarizedExperiment")
library("tidyverse")
library("EnhancedVolcano")
library("here")
library("sessioninfo")
library("ggrepel")

#### Set up ####

## dirs
plot_dir <- here("plots", "09_bulk_DE", "07_DE_plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data", "00_data_prep", "bulk_colors.Rdata"), verbose = TRUE)
# library_combo_colors
# library_prep_colors
# library_type_colors

#### load data ####

## Simple model ~library_type + Sample 
load(here("processed-data", "09_bulk_DE", "04_DE_library-type", "DE_library-type_simpleMod.Rdata"), verbose = TRUE)
DE_library_type_simple <- DE_out

DE_library_type_simple$gene$Bulk

table(DE_library_type_simple$gene$Bulk$P.Val < 0.05)
# FALSE  TRUE 
# 2462 19283 
table(DE_library_type_simple$gene$Bulk$P.Val < 0.01)
# FALSE  TRUE 
# 2462 19283 

## pval distribution 
pval_dist <- DE_library_type$gene$Bulk |>
  ggplot(aes(x = P.Value)) +
  geom_density()

ggsave(pval_dist, filename = here(plot_dir, "pval_density.png"))

volcano_libray_type <- EnhancedVolcano(DE_library_type_simple$gene$Bulk,
                                       lab = DE_library_type_simple$gene$Bulk$Symbol,
                                       x = 'logFC',
                                       y = 'adj.P.Val',
                                       title = "RiboZero vs. PolyA - Bulk",
                                       subtitle = "Model: ~library_type + Sample")

ggsave(volcano_libray_type, filename = here(plot_dir, "Bulk_DE_volcano_library-type_simple.png"), width = 10, height = 10)

load(here("processed-data", "09_bulk_DE", "04_DE_library-type", "DE_library-type.Rdata"), verbose = TRUE)

## MA plot

MA_plot_custom <- function(de_out = DE_library_type_simple$gene$Bulk,
                           title = NULL,
                           subtitle = NULL){
  de_out |>
    mutate(DE_class = case_when(logFC > 1 & adj.P.Val < 0.05 ~ "RiboZeroGold",
                                logFC < -1 & adj.P.Val < 0.05 ~ "polyA",
                                TRUE ~"Other")) |>
    ggplot(aes(x = log2(meanExprs), y = logFC, color = DE_class)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_text_repel(aes(label = ifelse(DE_class != "Other", Symbol, NA)), size = 2, show.legend=FALSE) +
    scale_color_manual(values = c(library_type_colors, "Other" = "darkgray")) +
    theme_bw() +
    geom_hline(yintercept = c(1,0,-1), linetype = c("dashed", "solid","dashed")) +
    labs(title = title, subtitle = subtitle)
  
} 


volcano_plot_custom <- function(de_out = DE_library_type_simple$gene$Bulk,
                                title = NULL,
                                subtitle = NULL){
  
  pval_lim <- max(de_out$P.Value[de_out$adj.P.Val < 0.05])
  
  de_out |>
    mutate(DE_class = case_when(logFC > 1 & adj.P.Val < 0.05 ~ "RiboZeroGold",
                                logFC < -1 & adj.P.Val < 0.05 ~ "polyA",
                                TRUE ~"Other")) |>
    ggplot(aes(x = logFC, y = -log10(P.Value), color = DE_class)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_text_repel(aes(label = ifelse(DE_class != "Other", Symbol, NA)), size = 2, show.legend=FALSE) +
    scale_color_manual(values = c(library_type_colors, "Other" = "darkgray")) +
    theme_bw() +
    geom_vline(xintercept = c(1,0,-1), linetype = c("dashed", "solid","dashed")) +
    geom_hline(yintercept = -log10(pval_lim), linetype = "dashed") +
    labs(title = title, subtitle = subtitle)
}

volcano_plot  <- volcano_plot_custom(de_out = DE_library_type_simple$gene$Bulk,
                                     title = "RiboZero vs. PolyA - Bulk",
                                     subtitle = "Model: ~library_type + Sample")

ggsave(volcano_plot, filename = here(plot_dir, "Volcano_test.png"))


## library type model

DE_library_type <- read.csv(here("processed-data", "09_bulk_DE", "04_DE_library-type", "DE_library-type_gene_Bulk.csv"))


max(DE_library_type$P.Value[DE_library_type$adj.P.Val > 0.05])

volcano_plot  <- volcano_plot_custom(de_out = DE_library_type,
                                     title = "RiboZero vs. PolyA - Bulk",
                                     subtitle = "Model: ~library_type + Sample + mitoRate + rRNA_rate + totalAssignedGene")

ggsave(volcano_plot, filename = here(plot_dir, "Bulk_DE_volcano_library-type_bulk.png"))

ma_plot  <- MA_plot_custom(de_out = DE_library_type,
                                     title = "RiboZero vs. PolyA - Bulk",
                                     subtitle = "Model: ~library_type + Sample + mitoRate + rRNA_rate + totalAssignedGene")

ggsave(ma_plot, filename = here(plot_dir, "Bulk_DE_MA_library-type_bulk.png"))


limma_ma_plot <- limma::plotMA(DE_library_type)
