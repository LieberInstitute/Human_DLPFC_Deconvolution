
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
plot_dir <- here("plots", "09_bulk_DE", "11_GO_analysis")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data", "00_data_prep", "bulk_colors.Rdata"), verbose = TRUE)
# library_combo_colors
# library_prep_colors
# library_type_colors

#### load data ####
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)

rd <- as.data.frame(rowData(rse_gene)) |> select(gencodeID, ensemblID, gene_type, Symbol, EntrezID)

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
  mutate(DE_class = case_when(logFC > 1 & adj.P.Val < 0.05 ~ "RiboZeroGold",
                                         logFC < -1 & adj.P.Val < 0.05 ~ "polyA",
                                         TRUE ~"None"),
         rank_p = row_number()) |>
  group_by(library_prep, DE_class) |>
  arrange(-abs(logFC)) |>
  mutate(rank_fc = ifelse(DE_class != "None", row_number(), NA))
 

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
})

DREAM_library_prep_long |> count(library_type, library_prep_pair)

DREAM_library_prep_long <- DREAM_library_prep_long |>
  as_tibble() |>
  separate(library_prep_pair, into = c("down", "up"), remove = FALSE) |>
  group_by(library_prep_pair, library_type) |>
  arrange(adj.P.Val) |>
  mutate(DE_class = case_when(logFC > 1 & adj.P.Val < 0.05 ~ up,
                              logFC < -1 & adj.P.Val < 0.05 ~ down,
                              TRUE ~"None"),
         rank_p = row_number()) |>
  group_by(library_type, library_prep_pair, DE_class) |>
  arrange(-abs(logFC)) |>
  mutate(rank_fc = ifelse(DE_class != "None", row_number(), NA))

DREAM_library_prep_long |> count()

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

#### Volcano Plots ####

type_max_pval <- DREAM_library_type_long |>
  filter(adj.P.Val < 0.05) |>
  arrange(-P.Value) |>
  slice(1) |>
  select(library_prep, P.Value, adj.P.Val) |>
  mutate(logP = -log10(P.Value))
 
library_type_volcano <- DREAM_library_type_long |>
  ggplot(aes(x = logFC, y = -log10(P.Value), color = DE_class)) +
  geom_point(size = .7, alpha = 0.5) +
  geom_text_repel(aes(label = ifelse(rank_p < 100 | rank_fc < 100, Symbol, NA)),
                  size = 2, 
                  show.legend=FALSE,
                  color = "black") +
  scale_color_manual(values = c(library_type_colors, "Other" = "darkgray")) +
  # geom_vline(xintercept = rep(c(1,0,-1), 3), linetype = rep(c("dashed", "solid","dashed"),3)) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  geom_hline(data = type_max_pval, aes(yintercept = logP), linetype = "dashed") +
  facet_wrap(~library_prep, nrow = 1) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(library_type_volcano, filename = here(plot_dir, "library_type_Volcano_small.png"), width = 7 , height = 5)
ggsave(library_type_volcano, filename = here(plot_dir, "library_type_Volcano.png"), width = 12)

DREAM_library_type_long |>
  count(gene_type) |>
  arrange(library_prep, DE_class, -n) |>
  print(n = 100)

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
                  size = 2, 
                  show.legend=FALSE,
                  color = "black") +
  scale_color_manual(values = c(Bulk = "#688E26", Cyto = "#FBB337", Nuc = "#104F55", Other = "darkgray")) +
  # geom_vline(xintercept = rep(c(1,0,-1), 3), linetype = rep(c("dashed", "solid","dashed"),3)) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  geom_hline(data = prep_max_pval, aes(yintercept = logP), linetype = "dashed") +
  facet_grid(library_type~library_prep_pair) +
  theme_bw() +
  theme(legend.position = "bottom")

# ggsave(library_prep_volcano, filename = here(plot_dir, "library_type_Volcano_small.png"), width = 7 , height = 5)
ggsave(library_prep_volcano, filename = here(plot_dir, "library_prep_Volcano.png"), height = 12, width = 12)


#### Upset Plot ####
DREAM_library_type_filter <- DREAM_library_type_long |> 
  filter(DE_class != "None") |>
  mutate(prep_class = paste0(DE_class, "_", library_prep))

DREAM_library_type_filter |>
  count(prep_class)
# # Groups:   library_prep, DE_class [6]
# library_prep DE_class     prep_class            n
# <chr>        <chr>        <chr>             <int>
#   1 Bulk         RiboZeroGold RiboZeroGold_Bulk   996
# 2 Bulk         polyA        polyA_Bulk         1005
# 3 Cyto         RiboZeroGold RiboZeroGold_Cyto  7109
# 4 Cyto         polyA        polyA_Cyto         4949
# 5 Nuc          RiboZeroGold RiboZeroGold_Nuc   5821
# 6 Nuc          polyA        polyA_Nuc          4084

DE_libray_type_geneList <- map(splitit(DREAM_library_type_filter$prep_class), ~DREAM_library_type_filter$ensemblID[.x])
map_int(DE_libray_type_geneList, length)

pdf(here(plot_dir, "library_type_upset.pdf"))
# upset(fromList(DE_libray_type_geneList), order.by = "freq", nsets = 6, keep.order = TRUE)
upset(fromList(DE_libray_type_geneList), order.by = "freq", sets = names(DE_libray_type_geneList), keep.order = TRUE)
dev.off()

## library_prep
DREAM_library_prep_filter <- DREAM_library_prep_long |> 
  filter(DE_class != "None") |>
  # mutate(prep_class = paste0(library_type  ,"_", library_prep_pair, "_", DE_class))
  mutate(prep_class = ifelse(DE_class == up, paste0(library_type,"_",up,"^_",down), paste0(library_type,"_",up,"_",down,"^")))

DREAM_library_prep_filter |>
  count(prep_class)
# library_type library_prep_pair DE_class prep_class                  n
# <chr>        <chr>             <chr>    <chr>                   <int>
#   1 RiboZeroGold Bulk_Cyto         Bulk     RiboZeroGold_Cyto_Bulk^  2698
# 2 RiboZeroGold Bulk_Cyto         Cyto     RiboZeroGold_Cyto^_Bulk   247
# 3 RiboZeroGold Bulk_Nuc          Bulk     RiboZeroGold_Nuc_Bulk^    346
# 4 RiboZeroGold Bulk_Nuc          Nuc      RiboZeroGold_Nuc^_Bulk    217
# 5 RiboZeroGold Cyto_Nuc          Cyto     RiboZeroGold_Nuc_Cyto^   1887
# 6 RiboZeroGold Cyto_Nuc          Nuc      RiboZeroGold_Nuc^_Cyto   2775
# 7 polyA        Bulk_Cyto         Bulk     polyA_Cyto_Bulk^         3269
# 8 polyA        Bulk_Cyto         Cyto     polyA_Cyto^_Bulk         2639
# 9 polyA        Bulk_Nuc          Bulk     polyA_Nuc_Bulk^           556
# 10 polyA        Bulk_Nuc          Nuc      polyA_Nuc^_Bulk           848
# 11 polyA        Cyto_Nuc          Cyto     polyA_Nuc_Cyto^          2449
# 12 polyA        Cyto_Nuc          Nuc      polyA_Nuc^_Cyto          2321

DE_libray_prep_geneList <- map(splitit(DREAM_library_prep_filter$prep_class), ~DREAM_library_prep_filter$ensemblID[.x])
map_int(DE_libray_prep_geneList, length)

pdf(here(plot_dir, "library_prep_upset.pdf"))
# upset(fromList(DE_libray_type_geneList), order.by = "freq", nsets = 6, keep.order = TRUE)
upset(fromList(DE_libray_prep_geneList), order.by = "freq", sets = names(DE_libray_prep_geneList), keep.order = TRUE)
dev.off()

