
library("SummarizedExperiment")
library("tidyverse")
library("EnhancedVolcano")
library("here")
library("sessioninfo")
library("ggrepel")
library("jaffelab")
# library("UpSetR")
library("ComplexUpset")

#### Set up ####

## dirs
plot_dir <- here("plots", "09_bulk_DE", "10_DREAM_plots")
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
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)
marker_genes_top25_simple <- marker_stats |> filter(rank_ratio <= 25) |> select(ensemblID = gene, cellType.target)
marker_genes_top25_simple |> count(cellType.target)

## library_type
load(here("processed-data", "09_bulk_DE","08_DREAM_library-type", "DREAM_library-type.Rdata"), verbose = TRUE)
head(DREAM_library_type$Bulk)
identical(rownames(rd), rownames(DREAM_library_type$Bulk))

DREAM_library_type_long <- map2_dfr(DREAM_library_type, names(DREAM_library_type), function(DE, name){
  DE <- cbind(DE, rd)
  DE$library_prep <- name
  return(DE)
}) |>
  as_tibble() |>
  group_by(library_prep) |>
  arrange(adj.P.Val) |>
  mutate(DE_class = case_when(logFC > 1 & adj.P.Val < 0.05 ~ paste0("RiboZeroGold_",library_prep),
                                         logFC < -1 & adj.P.Val < 0.05 ~ paste0("polyA_", library_prep),
                                         TRUE ~"None"),
         rank_p = row_number()) |>
  group_by(library_prep, DE_class) |>
  arrange(-abs(logFC)) |>
  mutate(rank_fc = ifelse(DE_class != "None", row_number(), NA),
         library_prep = factor(library_prep, levels = c("Cyto", "Bulk", "Nuc")),
         DE_class = gsub("Bulk", "Total", DE_class),
         rna_extract = factor(gsub("Bulk", "Total", library_prep), levels = c("Cyto", "Total", "Nuc"))
         ) |> 
  left_join(marker_genes_top25_simple)

DREAM_library_type_long |> count()

## LibraryPrep
load(here("processed-data", "09_bulk_DE","09_DREAM_library-prep", "DREAM_library-prep.Rdata"), verbose = TRUE)
names(DREAM_library_prep$polyA)
# [1] "Cyto_Nuc"  "Bulk_Nuc"  "Bulk_Cyto"

identical(rownames(rd), rownames(DREAM_library_prep$polyA$Cyto_Nuc))

DREAM_library_prep_long <- map2_dfr(DREAM_library_prep, names(DREAM_library_prep), function(DE_list, type_name){
  DE <- map2_dfr(DE_list, names(DE_list), function(DE, name){
    DE <- cbind(DE, rd)
    DE$library_prep_pair <- name
    return(DE)
  })
  DE$library_type <- type_name
  return(DE)
}) |>
  as_tibble() |>
  mutate(library_prep_pair = gsub("Bulk", "Total", library_prep_pair)) |>
  separate(library_prep_pair, into = c("down", "up"), remove = FALSE) |>
  group_by(library_prep_pair, library_type) |>
  arrange(adj.P.Val) |>
  mutate(DE_class = case_when(logFC > 1 & adj.P.Val < 0.05 ~ paste0(library_type, "_", up),
                              logFC < -1 & adj.P.Val < 0.05 ~ paste0(library_type, "_", down),
                              TRUE ~"None"),
         rank_p = row_number()) |>
  group_by(library_type, library_prep_pair, DE_class) |>
  arrange(-abs(logFC)) |>
  mutate(rank_fc = ifelse(DE_class != "None", row_number(), NA)) |> 
  left_join(marker_genes_top25_simple)

DREAM_library_prep_long |> count()
DREAM_library_prep_long |> count(down, up)

#### pval distribution ####
library_type_pval_histo <- DREAM_library_type_long |>
  ggplot(aes(x = P.Value)) +
  geom_histogram(binwidth = 0.05) +
  facet_wrap(~library_prep, ncol = 1) +
  theme_bw()

ggsave(library_type_pval_histo, filename = here(plot_dir, "library_type_pval_histogram.png"))

library_prep_pval_histo <- DREAM_library_prep_long |>
  ggplot(aes(x = P.Value)) +
  geom_histogram(binwidth = 0.05) +
  facet_grid(library_type~library_prep_pair) +
  theme_bw()

ggsave(library_prep_pval_histo, filename = here(plot_dir, "library_prep_pval_histogram.png"))

#### DEG counts ####

DREAM_library_type_long |> count(DE_class)

DREAM_library_type_count <- DREAM_library_type_long |> 
  group_by(rna_extract, DE_class) |>
  summarise(n_DE = n()) |>
  group_by(rna_extract) |>
  mutate(percent = 100*n_DE/sum(n_DE))

write_csv(DREAM_library_type_count, here("processed-data", "09_bulk_DE","08_DREAM_library-type", "DREAM_library-type_summary_FDR05.csv"))

# rna_extract DE_class            n_DE percent
# <chr>       <chr>              <int>   <dbl>
# 1 Cyto        None                9687   44.5 
# 2 Cyto        RiboZeroGold_Cyto   7109   32.7 
# 3 Cyto        polyA_Cyto          4949   22.8 
# 4 Nuc         None               11840   54.4 
# 5 Nuc         RiboZeroGold_Nuc    5821   26.8 
# 6 Nuc         polyA_Nuc           4084   18.8 
# 7 Total       None               19744   90.8 
# 8 Total       RiboZeroGold_Total   996    4.58
# 9 Total       polyA_Total         1005    4.62
  
DREAM_library_type_anno <- DREAM_library_type_count |>
  filter(DE_class != "None") |>
  mutate(upreg = grepl("RiboZeroGold", DE_class),
         anno_n = sprintf("\n%1d (%.1f%%)", n_DE, percent),
         anno_d = ifelse(upreg, 
                         "RiboZeroGold →", 
                         "← polyA"),
         anno = paste0(anno_d, anno_n),
         anno_x = ifelse(upreg, 1, -1),
         anno_y = 0) 


## prep 
DREAM_library_prep_count <- DREAM_library_prep_long |> 
  group_by(library_type, library_prep_pair, DE_class, up, down) |>
  summarise(n_DE = n()) |>
  group_by(library_type, library_prep_pair) |>
  mutate(percent = 100*n_DE/sum(n_DE))

write_csv(DREAM_library_prep_count, file = here("processed-data", "09_bulk_DE","09_DREAM_library-prep", "DREAM_library-prep_summary_FDR05.csv"))

# library_type library_prep_pair DE_class           n_DE percent
# <chr>        <chr>             <chr>             <int>   <dbl>
# 1 RiboZeroGold Cyto_Nuc          None               17083  78.6  
# 2 RiboZeroGold Cyto_Nuc          RiboZeroGold_Cyto   1887   8.68 
# 3 RiboZeroGold Cyto_Nuc          RiboZeroGold_Nuc    2775  12.8  
# 4 RiboZeroGold Total_Cyto        None               18800  86.5  
# 5 RiboZeroGold Total_Cyto        RiboZeroGold_Cyto    247   1.14 
# 6 RiboZeroGold Total_Cyto        RiboZeroGold_Total  2698  12.4  
# 7 RiboZeroGold Total_Nuc         None               21182  97.4  
# 8 RiboZeroGold Total_Nuc         RiboZeroGold_Nuc     217   0.998
# 9 RiboZeroGold Total_Nuc         RiboZeroGold_Total   346   1.59 
# 10 polyA        Cyto_Nuc          None               16975  78.1  
# 11 polyA        Cyto_Nuc          polyA_Cyto          2449  11.3  
# 12 polyA        Cyto_Nuc          polyA_Nuc           2321  10.7  
# 13 polyA        Total_Cyto        None               15837  72.8  
# 14 polyA        Total_Cyto        polyA_Cyto          2639  12.1  
# 15 polyA        Total_Cyto        polyA_Total         3269  15.0  
# 16 polyA        Total_Nuc         None               20341  93.5  
# 17 polyA        Total_Nuc         polyA_Nuc            848   3.90 
# 18 polyA        Total_Nuc         polyA_Total          556   2.56


DREAM_library_type_anno <- DREAM_library_prep_count |>
  filter(DE_class != "None") |>
  mutate(upreg = grepl(up, DE_class),
         anno_n = sprintf("\n%1d (%.1f%%)", n_DE, percent),
         anno_d = ifelse(upreg, 
                        paste(up, "→"), 
                        paste("←", down)),
         anno = paste0(anno_d, anno_n),
         anno_x = ifelse(upreg, 1, -1),
         anno_y = 0) 

#### Volcano Plots ####

type_max_pval <- DREAM_library_type_long |>
  ungroup() |>
  group_by(rna_extract) |>
  filter(adj.P.Val < 0.05) |>
  arrange(-P.Value) |>
  slice(1) |>
  select(rna_extract, P.Value, adj.P.Val) |>
  mutate(logP = -log10(P.Value))
 
library_type_volcano <- DREAM_library_type_long |>
  ggplot(aes(x = logFC, y = -log10(P.Value), color = DE_class)) +
  geom_point(size = .7, alpha = 0.5) +
  geom_text_repel(aes(label = ifelse(rank_p < 100 | rank_fc < 100, Symbol, NA)),
                  size = 2, 
                  show.legend=FALSE,
                  color = "black") +
  scale_color_manual(values = c(library_combo_colors2, "Other" = "darkgray")) +
  # geom_vline(xintercept = rep(c(1,0,-1), 3), linetype = rep(c("dashed", "solid","dashed"),3)) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  geom_hline(data = type_max_pval, aes(yintercept = logP), linetype = "dashed") +
  geom_text(data = DREAM_library_type_anno, aes(label = anno, x = 15*anno_x, y = anno_y), color = "black", size = 2) +
  facet_wrap(~rna_extract, nrow = 1) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(library_type_volcano, filename = here(plot_dir, "library_type_Volcano_small.png"), width = 6 , height = 4)
ggsave(library_type_volcano, filename = here(plot_dir, "library_type_Volcano.png"), width = 12)

## library prep
prep_max_pval <- DREAM_library_prep_long |>
  filter(adj.P.Val < 0.05) |>
  arrange(-P.Value) |>
  slice(1) |>
  select(library_type, library_prep_pair, P.Value, adj.P.Val) |>
  mutate(logP = -log10(P.Value))

library_prep_volcano <- DREAM_library_prep_long |>
  ggplot(aes(x = logFC, y = -log10(P.Value), color = DE_class)) +
  geom_point(size = .7, alpha = 0.5) +
  geom_text_repel(aes(label = ifelse(rank_p < 100 | rank_fc < 100, Symbol, NA)),
                  size = 3, 
                  show.legend=FALSE,
                  color = "black") +
  scale_color_manual(values = c(library_combo_colors2, "Other" = "darkgray")) +
  # scale_color_manual(values = c(Bulk = "#688E26", Cyto = "#FBB337", Nuc = "#104F55", Other = "darkgray")) +
  # geom_vline(xintercept = rep(c(1,0,-1), 3), linetype = rep(c("dashed", "solid","dashed"),3)) +
  geom_text(data = DREAM_library_type_anno, aes(label = anno2, x = 5*anno_x, y = anno_y), color = "black", size = 3) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  geom_hline(data = prep_max_pval, aes(yintercept = logP), linetype = "dashed") +
  facet_grid(library_type~library_prep_pair) +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size=14))

ggsave(library_prep_volcano, filename = here(plot_dir, "library_prep_Volcano.png"), height = 8.5, width = 11)


#### Upset Plot ####
DREAM_library_type_filter <- DREAM_library_type_long |> 
  filter(DE_class != "None") 

DE_libray_type_geneList <- map(splitit(DREAM_library_type_filter$DE_class), ~DREAM_library_type_filter$ensemblID[.x])
map_int(DE_libray_type_geneList, length)

pdf(here(plot_dir, "library_type_upset.pdf"), width = 9)
# upset(fromList(DE_libray_type_geneList), order.by = "freq", nsets = 6, keep.order = TRUE)
ups <- UpSetR::upset(fromList(DE_libray_type_geneList), 
      order.by = "freq", 
      sets = names(DE_libray_type_geneList), 
      keep.order = TRUE,
      sets.bar.color = library_combo_colors2[names(DE_libray_type_geneList)]
                                        )
print(ups)
dev.off()

## library_prep

case_when(logFC > 1 & adj.P.Val < 0.05 ~ paste0(library_type, "_", up),
                            logFC < -1 & adj.P.Val < 0.05 ~ paste0(library_type, "_", down),
                            TRUE ~"None")

DREAM_library_prep_filter <- DREAM_library_prep_long |> 
  filter(DE_class != "None") |>
  # mutate(prep_class = paste0(library_type  ,"_", library_prep_pair, "_", DE_class))
  # mutate(prep_class = ifelse(grepl(DE_class, up), paste0(library_type,"_",up,"^_",down), paste0(library_type,"_",up,"_",down,"^")))
  mutate(prep_class = paste0(library_type,"_",gsub("_", "vs",library_prep_pair)),
         prep_class2 = case_when(logFC > 1 & adj.P.Val < 0.05 ~ paste0(prep_class, "^", up),
                                 logFC < -1 & adj.P.Val < 0.05 ~ paste0(prep_class, "^", down),
                                 TRUE ~ NA), 
  )

(prep_class_count <- DREAM_library_prep_filter |>
  count(up, prep_class2))
# library_type library_prep_pair DE_class           up    prep_class2                        n
# <chr>        <chr>             <chr>              <chr> <chr>                          <int>
#   1 RiboZeroGold Cyto_Nuc          RiboZeroGold_Cyto  Nuc   RiboZeroGold_CytovsNuc_Cyto     1887
# 2 RiboZeroGold Cyto_Nuc          RiboZeroGold_Nuc   Nuc   RiboZeroGold_CytovsNuc_Nuc      2775
# 3 RiboZeroGold Total_Cyto        RiboZeroGold_Cyto  Cyto  RiboZeroGold_TotalvsCyto_Cyto    247
# 4 RiboZeroGold Total_Cyto        RiboZeroGold_Total Cyto  RiboZeroGold_TotalvsCyto_Total  2698
# 5 RiboZeroGold Total_Nuc         RiboZeroGold_Nuc   Nuc   RiboZeroGold_TotalvsNuc_Nuc      217
# 6 RiboZeroGold Total_Nuc         RiboZeroGold_Total Nuc   RiboZeroGold_TotalvsNuc_Total    346
# 7 polyA        Cyto_Nuc          polyA_Cyto         Nuc   polyA_CytovsNuc_Cyto            2449
# 8 polyA        Cyto_Nuc          polyA_Nuc          Nuc   polyA_CytovsNuc_Nuc             2321
# 9 polyA        Total_Cyto        polyA_Cyto         Cyto  polyA_TotalvsCyto_Cyto          2639
# 10 polyA        Total_Cyto        polyA_Total        Cyto  polyA_TotalvsCyto_Total         3269
# 11 polyA        Total_Nuc         polyA_Nuc          Nuc   polyA_TotalvsNuc_Nuc             848
# 12 polyA        Total_Nuc         polyA_Total        Nuc   polyA_TotalvsNuc_Total           556

DE_libray_prep_geneList <- map(splitit(DREAM_library_prep_filter$prep_class2), ~DREAM_library_prep_filter$ensemblID[.x])
map_int(DE_libray_prep_geneList, length)
unique(DREAM_library_prep_filter$DE_class)

pdf(here(plot_dir, "library_prep_upset.pdf"), width = 11)
# upset(fromList(DE_libray_type_geneList), order.by = "freq", nsets = 6, keep.order = TRUE)
ups <- UpSetR::upset(fromList(DE_libray_prep_geneList), 
                     order.by = "freq", 
                     sets = names(DE_libray_prep_geneList), 
                     keep.order = TRUE,
                     sets.bar.color = library_combo_colors2[prep_class_count$DE_class])
print(ups)
dev.off()

#### marker genes ####
DREAM_library_type_long <- DREAM_library_type_long |> 
  left_join(marker_genes_top25_simple)

DREAM_library_type_long |> filter(DE_class != "None") |> count(cellType.target)
# library_prep DE_class     cellType.target     n
# <chr>        <chr>        <fct>           <int>
#   1 Bulk         RiboZeroGold Astro               1
# 2 Bulk         RiboZeroGold Micro               1
# 3 Bulk         RiboZeroGold NA                994
# 4 Bulk         polyA        EndoMural           1
# 5 Bulk         polyA        Oligo               4
# 6 Bulk         polyA        Excit               1
# 7 Bulk         polyA        Inhib               1
# 8 Bulk         polyA        NA                998
# 9 Cyto         RiboZeroGold Astro              15
# 10 Cyto         RiboZeroGold EndoMural          14


DREAM_type_marker_count <- DREAM_library_type_long |> 
  count(cellType.target) |>
  filter(!is.na(cellType.target)) |>
  separate(DE_class, into = c("library_type", NA), remove = FALSE, sep = "_") |>
  mutate(library_type = factor(library_type, levels =c("polyA", "None", "RiboZeroGold")))


DREAM_type_marker_count_tile <- DREAM_type_marker_count |> 
  ggplot(aes(x = library_type, y = cellType.target, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), color = "white") +
  facet_wrap(~library_prep) +
  theme_bw()

ggsave(DREAM_type_marker_count_tile, filename = here(plot_dir, "DREAM_type_marker_count_tile.png"))
  
## prep
DREAM_library_prep_long <- DREAM_library_prep_long |> 
  left_join(marker_genes_top25_simple)

DREAM_library_prep_long |> 
  filter(DE_class != "None", !is.na(cellType.target)) |> 
  count(cellType.target) |>
  arrange(-n) |>
  print(n = 35)
  
## plots

library_type_volcano_ct <- DREAM_library_type_long |>
  filter(!is.na(cellType.target)) |>
  ggplot(aes(x = logFC, y = -log10(P.Value), color = cellType.target)) +
  geom_point(size = .7, alpha = 0.5) +
  geom_text_repel(aes(label = ifelse(!is.na(cellType.target), Symbol, NA)),
                  size = 2, 
                  show.legend=FALSE,
                  color = "black") +
  scale_color_manual(values = cell_type_colors_broad) +
  geom_text(data = DREAM_library_type_anno, aes(label = anno_d, x = 4*anno_x, y = anno_y), color = "black", size = 5) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  geom_hline(data = type_max_pval, aes(yintercept = logP), linetype = "dashed") +
  facet_wrap(~rna_extract, nrow = 1) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(library_type_volcano_ct, filename = here(plot_dir, "library_type_Volcano_cellType_markers.png"), width = 12)


library_prep_volcano_ct <- DREAM_library_prep_long |>
  filter(!is.na(cellType.target)) |>
  ggplot(aes(x = logFC, y = -log10(P.Value), color = cellType.target)) +
  geom_point(size = .7, alpha = 0.5) +
  geom_text_repel(aes(label = ifelse(!is.na(cellType.target), Symbol, NA)),
                  size = 2, 
                  show.legend=FALSE,
                  color = "black") +
  scale_color_manual(values = cell_type_colors_broad) +
  geom_text(data = DREAM_library_type_anno, aes(label = anno_d, x = 1.5*anno_x, y = anno_y), color = "black", size = 5) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  geom_hline(data = prep_max_pval, aes(yintercept = logP), linetype = "dashed") +
  facet_grid(library_type~library_prep_pair) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(library_prep_volcano_ct, filename = here(plot_dir, "library_prep_Volcano_cellType_markers.png"), height = 12, width = 12)

