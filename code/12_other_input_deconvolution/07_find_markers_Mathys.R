library("DeconvoBuddies")
library("SingleCellExperiment")
library("tidyverse")
library("ggrepel")
library("here")
library("sessioninfo")

#### prep dirs ####
data_dir <- here("processed-data", "12_other_input_deconvolution", "07_find_markers_Mathys")
if(!dir.exists(data_dir)) dir.create(data_dir)

plot_dir <- here("plots", "12_other_input_deconvolution", "07_find_markers_Mathys")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

## load data
load(here("processed-data", "12_other_input_deconvolution", "sce_Mathys.Rdata"), verbose = TRUE)
table(sce$cellType_broad)
table(is.na(sce$cellType_broad))
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib 
# 3392       288      1920     18235      2627     34976      9196 

## load bulk data & find common genes 
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
common_genes <- intersect(rowData(sce)$gene_id, rowData(rse_gene)$ensemblID)
message("Common genes with bulk data: ", length(common_genes))
rm(rse_gene)

#### Run Marker Finding ####
message(Sys.time(), " Find Mean Ratio Genes")
mean_ratio <- get_mean_ratio2(sce, cellType_col = "cellType_broad", add_symbol = TRUE)

message(Sys.time(), " Find 1vAll genes")
markers_1vALL <- findMarkers_1vAll(sce, cellType_col = "cellType_broad", mod = "~individualID")

marker_stats <- mean_ratio |> 
  left_join(markers_1vALL) |>
  mutate(in_bulk = gene %in% common_genes,
         MeanRatio_top25 = (rank_ratio <= 25 & in_bulk))

message(Sys.time(), " Save marker gene data")
save(marker_stats, file = here(data_dir, "Mathys_marker_stats_broad.Rdata"))
write_csv(marker_stats, file = here(data_dir, "Mathys_marker_stats_broad.csv"))

## pull marker sets 
marker_stats |> 
  dplyr::filter(MeanRatio_top25) |>
  dplyr::count(cellType.target)

markers_mean_ratio_top25 <- marker_stats |> 
  dplyr::filter(MeanRatio_top25) |>
  dplyr::pull(gene)

cat(markers_mean_ratio_top25, sep = "\n", file = here("processed-data","12_other_input_deconvolution", "Mathys_markers_MeanRatio_top25.txt"))

## hockey stick plots 
load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)

mean_ratio_v_logFC <- marker_stats |>
  mutate(top25 = rank_ratio <= 25) |>
  ggplot(aes(x = ratio, y = std.logFC, color = top25)) +
  geom_point(alpha = 0.5) +
  ggrepel::geom_text_repel(aes(label = ifelse(rank_ratio ==1 | rank_marker ==1 , Symbol, "")), 
                  size = 2.5, 
                  color = "black") +
  facet_wrap(~cellType.target, scales = "free_x", nrow = 1) +
  labs(x = "Mean Ratio", y = "1vALL Standard logFC") +
  theme_bw()

ggsave(mean_ratio_v_logFC, filename = here(plot_dir, "Mathys_mean_ratio_v_logFC_free.png"), width = 12, height = 2.5)
ggsave(mean_ratio_v_logFC, filename = here(plot_dir, "Mathys_mean_ratio_v_logFC_free.pdf"), width = 12, height = 2.5)

# slurmjobs::job_single('07_find_markers_Mathys', create_shell = TRUE, memory = '25G', command = "Rscript 07_find_markers_Mathys")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
