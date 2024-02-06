
library("SummarizedExperiment")
library("tidyverse")
library("EnhancedVolcano")
library("here")
library("sessioninfo")
library("ggrepel")
library("jaffelab")
library("UpSetR")

#### Set up ####

## dirs
plot_dir <- here("plots", "10_bulk_vs_sn_DE", "04_DREAM_plots_sn")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data", "00_data_prep", "bulk_colors.Rdata"), verbose = TRUE)
# library_combo_colors
# library_prep_colors
# library_type_colors
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo
# cell_type_colors_broad
#### load data ####
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)

rd <- as.data.frame(rowData(rse_gene)) |> select(gencodeID, ensemblID, gene_type, Symbol, EntrezID)

## marker gene data
load(here("processed-data", "06_marker_genes", "marker_genes_top25.Rdata"), verbose = TRUE)
marker_genes_top25_simple <- marker_genes_top25_simple |> rename(ensemblID = gene)

## library_type
load(here("processed-data", "10_bulk_vs_sn_DE",  "03_DREAM_sn_v_bulk", "DREAM_data-type.Rdata"), verbose = TRUE)
DREAM_data_type <- DREAM_library_type ## fix object name

head(DREAM_data_type$polyA)

rownames(rd) <- rd$ensemblID
rd <- rd[rownames(DREAM_data_type$polyA),]
identical(rownames(rd), rownames(DREAM_data_type$polyA))

DREAM_data_type_long <- map2_dfr(DREAM_data_type, names(DREAM_data_type), function(DE, name){
  DE <- cbind(DE, rd)
  DE$library_type <- name
  return(DE)
}) |>
  as_tibble() |>
  group_by(library_type) |>
  arrange(adj.P.Val) |>
  mutate(DE_class = case_when(logFC > 1 & adj.P.Val < 0.05 ~ "snRNA-seq",
                                         logFC < -1 & adj.P.Val < 0.05 ~ library_type,
                                         TRUE ~"None"),
         rank_p = row_number()) |>
  group_by(library_type, DE_class) |>
  arrange(-abs(logFC)) |>
  mutate(rank_fc = ifelse(DE_class != "None", row_number(), NA)) |> 
  left_join(marker_genes_top25_simple)
 

DREAM_data_type_long |> 
  group_by(library_type, DE_class) |>
  summarise(n_DE = n()) |>
  group_by(library_type) |>
  mutate(percent = 100*n_DE/sum(n_DE))

# library_type DE_class      n_DE percent
# <chr>        <chr>        <int>   <dbl>
# 1 RiboZeroGold None          8902    50.4
# 2 RiboZeroGold RiboZeroGold  4423    25.0
# 3 RiboZeroGold snRNA-seq     4335    24.5
# 4 polyA        None          6402    36.3
# 5 polyA        polyA         5509    31.2
# 6 polyA        snRNA-seq     5749    32.6

DREAM_data_type_long |> filter(DE_class != "None") |> count(DE_class, cellType.target)

#### pval distribution ####
data_type_pval_histo <- DREAM_data_type_long |>
  ggplot(aes(x = P.Value)) +
  geom_histogram(binwidth = 0.05) +
  facet_wrap(~library_type, ncol = 1) +
  theme_bw()

ggsave(data_type_pval_histo, filename = here(plot_dir, "data_type_pval_histogram.png"))

#### Volcano Plots ####

type_max_pval <- DREAM_data_type_long |>
  group_by(library_type) |>
  filter(adj.P.Val < 0.05) |>
  arrange(-P.Value) |>
  slice(1) |>
  select(library_type, P.Value, adj.P.Val) |>
  mutate(logP = -log10(P.Value))
 
data_type_volcano <- DREAM_data_type_long |>
  ggplot(aes(x = logFC, y = -log10(P.Value), color = DE_class)) +
  geom_point(size = .7, alpha = 0.5) +
  geom_text_repel(aes(label = ifelse(DE_class != "None" & ensemblID != Symbol, Symbol, NA)),
  # geom_text_repel(aes(label = ifelse(rank_p < 100 | rank_fc < 100 & ensemblID != Symbol, Symbol, NA)),
                  size = 2, 
                  show.legend=FALSE,
                  color = "black") +
  scale_color_manual(values = c(library_type_colors, "Other" = "darkgray", "snRNA-seq" = "#2f7ec0")) +
  # scale_color_manual(values = c("Bulk" = "#8D6E53", "snRNA-seq" = "#417B5A", "Other" = "darkgray")) +
  # geom_vline(xintercept = rep(c(1,0,-1), 3), linetype = rep(c("dashed", "solid","dashed"),3)) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  geom_hline(data = type_max_pval, aes(yintercept = logP), linetype = "dashed") +
  facet_wrap(~library_type, nrow = 1) +
  theme_bw() 

ggsave(data_type_volcano +
         theme(legend.position = "none"), 
       filename = here(plot_dir, "data_type_Volcano_small.png"), width = 4 , height = 4)
ggsave(data_type_volcano +
         theme(legend.position = "none"),
       filename = here(plot_dir, "data_type_Volcano_small.pdf"), width = 4 , height = 4)

ggsave(data_type_volcano +
         theme(legend.position = "bottom"), 
       filename = here(plot_dir, "data_type_Volcano.png"), width = 12)

DREAM_data_type_long |>
  count(gene_type) |>
  arrange(library_type, DE_class, -n) |>
  print(n = 100)

#### Upset Plot ####
DREAM_data_type_filter <- DREAM_data_type_long |> 
  filter(DE_class != "None") |>
  mutate(prep_class = paste0(DE_class, "_", library_type))

DREAM_data_type_filter |>
  count(prep_class)

DE_libray_type_geneList <- map(splitit(DREAM_data_type_filter$prep_class), ~DREAM_data_type_filter$ensemblID[.x])
map_int(DE_libray_type_geneList, length)

pdf(here(plot_dir, "data_type_upset.pdf"))
# upset(fromList(DE_libray_type_geneList), order.by = "freq", nsets = 6, keep.order = TRUE)
upset(fromList(DE_libray_type_geneList), order.by = "freq", sets = names(DE_libray_type_geneList), keep.order = TRUE)
dev.off()

## cell type markers
data_type_volcano_ct <- DREAM_data_type_long |>
  filter(!is.na(cellType.target)) |>
  ggplot(aes(x = logFC, y = -log10(P.Value), color = cellType.target)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(aes(label = ifelse(!is.na(cellType.target), Symbol, NA)),
                  size = 2, 
                  show.legend=FALSE,
                  color = "black") +
  scale_color_manual(values = cell_type_colors_broad) +
  # geom_vline(xintercept = rep(c(1,0,-1), 3), linetype = rep(c("dashed", "solid","dashed"),3)) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  geom_hline(data = type_max_pval, aes(yintercept = logP), linetype = "dashed") +
  facet_wrap(~library_type, nrow = 1) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(data_type_volcano_ct, filename = here(plot_dir, "data_type_Volcano_cellType_markers.png"), height = 12, width = 12)



# slurmjobs::job_single(name = "04_DREAM_plots_sn", memory = "5G", cores = 1, create_shell = TRUE, command = "Rscript 04_DREAM_plots_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
