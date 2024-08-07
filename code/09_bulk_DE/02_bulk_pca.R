library("SummarizedExperiment")
library("tidyverse")
library("sessioninfo")
library("here")
library("jaffelab")
library("recount")
library("viridis")
library("ggrepel")
library("GGally")
library("vsn")

## prep dirs ##
plot_dir <- here("plots", "09_bulk_DE", "02_bulk_pca")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data", "00_data_prep", "bulk_colors.Rdata"), verbose = TRUE)
# library_combo_colors
# library_prep_colors
# library_type_colors
load(here("processed-data", "00_data_prep", "library_combo_shapes.Rdata"), verbose = TRUE)
# library_combo_shapes
# library_combo_shapes2

#### Load Data ####
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
pd <- as.data.frame(colData(rse_gene))

# move to colData
pd <- as.data.frame(colData(rse_gene))

focused_qc_metrics <- c(
    "concordMapRate",
    "mitoRate",
    "numMapped",
    "numReads",
    "overallMapRate",
    "totalAssignedGene",
    "totalMapped"
)

pd_simple <- pd |> 
  select(SAMPLE_ID, Sample, round, batch, BrNum, Position, library_type, library_prep, library_combo, library_combo2, qc_class, sex, age, all_of(focused_qc_metrics)) |>
  mutate(round = as.factor(round))

#### check out dispersion patterns ####
pdf(here(plot_dir, "mean_vs_sd_counts.pdf"))
meanSdPlot(assays(rse_gene)$counts, ranks = FALSE)
dev.off()

gene_mean <- tibble(gene = rowData(rse_gene)$Symbol,
                    mean = rowMeans(assays(rse_gene)$counts),
                    sd = apply(assays(rse_gene)$counts, 1, sd))

gene_mean |> arrange(-mean)

mean_v_sd <- gene_mean |>
  ggplot(aes(mean, sd)) +
  geom_point(alpha = 0.2) +
  geom_abline() +
  geom_smooth(color = "red") +
  geom_text_repel(aes(label = ifelse(mean > 2e5, gene, "")), size = 2)

ggsave(mean_v_sd, filename = here(plot_dir, "my_mean_vs_sd_counts.png"))

#### filter genes by Expression ####
gene_rpkm <- getRPKM(rse_gene, "Length")
rse_gene_filter <- rse_gene[rowMeans(gene_rpkm) > 0.1, ]
gene_rpkm_filter <- gene_rpkm[rowMeans(gene_rpkm) > 0.1, ]
table(droplevels(seqnames(rse_gene_filter)))

# pdf(here(plot_dir, "mean_vs_sd_rpkm.pdf"))
# meanSdPlot(gene_rpkm, ranks = FALSE)
# dev.off()

# ## get expression
geneExprs_filter <- log2(gene_rpkm_filter + 1)
assays(rse_gene_filter)$logcounts <- log2(gene_rpkm_filter + 1) ## check


pdf(here(plot_dir, "mean_vs_sd_logrpkm.pdf"))
meanSdPlot(geneExprs_filter, ranks = FALSE)
dev.off()

## calc PCA and vars
set.seed(071524)

pca <- prcomp(t(geneExprs_filter))
pca_vars <- getPcaVars(pca)
pca_vars

pca$x[1:5,1:5]

save(pca, pca_vars, file = here("processed-data", "02_bulk_pca", "bulk_PCA.Rdata"))

pca_vars_lab <- paste0(
  "PC", seq(along = pca_vars), ": ",
  pca_vars, "% Var Expl"
)

## create table with groups and some QC metrics
pca_tab <- pd_simple |> cbind(pca$x[, 1:5])

#### ggpairs for pca ####
gg_pca <- ggpairs(pca_tab,
    mapping = aes(color = library_combo),
    columns = paste0("PC", 1:5),
    upper = "blank"
) +
    scale_color_manual(values = library_combo_colors) +
    scale_fill_manual(values = library_combo_colors)

ggsave(gg_pca, filename = here(plot_dir, "ggpairs_pca.png"), height = 10, width = 10)

gg_pca_position <- ggpairs(pca_tab,
    mapping = aes(color = Position),
    columns = paste0("PC", 1:5),
    upper = "blank"
)

ggsave(gg_pca_position, filename = here(plot_dir, "ggpairs_position.png"), height = 10, width = 10)

## PC1 vs. PC2
pc1v2 <- pca_tab |>
    ggplot(aes(x = PC1, y = PC2, color = library_combo)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = library_combo_colors) +
    labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[2]]) +
    coord_equal()

ggsave(pc1v2, filename = here(plot_dir, "Bulk_PC1vPC2_library_combo.png"))
ggsave(pc1v2 + theme(legend.position = "None"), 
       filename = here(plot_dir, "Bulk_PC1vPC2_library_combo_small.png"),
       width = 3, height = 4)

ggsave(pc1v2 + theme(legend.position = "None"), 
       filename = here(plot_dir, "Bulk_PC1vPC2_library_combo_small.pdf"),
       width = 3, height = 4)

pc1v2_lab <- pca_tab |>
    ggplot(aes(x = PC1, y = PC2, color = library_combo, shape = qc_class)) +
    geom_point() +
    geom_text_repel(aes(label = Sample), color = "black", size = 2) +
    theme_bw() +
    scale_color_manual(values = library_combo_colors) +
    labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[2]]) +
    coord_equal()

ggsave(pc1v2_lab, filename = here(plot_dir, "Bulk_PC1vPC2_library_combo_label.png"))

pc2v5_lab <- pca_tab |>
    ggplot(aes(x = PC2, y = PC5, color = library_combo, shape = qc_class)) +
    geom_point() +
    geom_text_repel(aes(label = Sample), color = "black", size = 2) +
    theme_bw() +
    scale_color_manual(values = library_combo_colors) +
    labs(x = pca_vars_lab[[2]], y = pca_vars_lab[[5]]) +
    coord_equal()

ggsave(pc2v5_lab, filename = here(plot_dir, "Bulk_PC2vPC5_library_combo_label.png"))


pc_test <- pca_tab |>
    ggplot(aes(x = PC1, y = PC2, color = BrNum, shape = qc_class)) +
    geom_point() +
    theme_bw() +
    # scale_color_manual(values = library_combo_colors) +
    labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[2]])

ggsave(pc_test, filename = here(plot_dir, "Bulk_PC1vPC2_BrNum.png"))

pc_test <- pca_tab |>
    ggplot(aes(x = PC1, y = PC2, color = mitoRate, shape = qc_class)) +
    geom_point() +
    theme_bw() +
    scale_color_viridis() +
    labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[2]])

ggsave(pc_test, filename = here(plot_dir, "Bulk_PC1vPC2_mitoRate.png"))


pc_round <- pca_tab |>
  ggplot(aes(x = PC1, y = PC2, color = round, shape = library_combo2)) +
  geom_point() +
  theme_bw() +
  scale_shape_manual(values = library_combo_shapes2, "Library Combo") +
  labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[2]])

ggsave(pc_round , filename = here(plot_dir, "Bulk_PC1vPC2_round.png"), width = 7)


#### compare directly to PD variables ####

pca_long <- pca$x[, 1:6] |>
    as.data.frame() |>
    rownames_to_column("SAMPLE_ID") |>
    pivot_longer(!SAMPLE_ID, names_to = "PC_name", values_to = "PC_val") |>
    right_join(pd) |>
    left_join(tibble(
        PC_lab = pca_vars_lab,
        PC_name = paste0("PC", seq(along = pca_vars))
    )) |>
   mutate(round = as.factor(round))

pca_long

# Position, library_type, library_prep, qc_class

pdf(here(plot_dir, "PCs_vs_groups.pdf"), width = 10)
walk(c("BrNum", "Sample", "round","sex", "Position", "library_type", "library_prep", "library_combo", "batch", "qc_class"), ~ {
    pca_v_cat <- pca_long |>
        ggplot(aes(y = PC_val, x = .data[[.x]], color = .data[[.x]])) +
        geom_boxplot() +
        # geom_text_repel(aes(label = ifelse(qc_class == "warn", Sample,"")), size = 2, color = "black") +
        facet_wrap(~PC_lab) +
        labs(title = .x) +
        theme_bw()  +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # ggsave(pca_v_con, filename = here(plot_dir, paste0("pca_v_",.x,".png")), width = 10)

    print(pca_v_cat)
})
dev.off()

# mitoRate,totalAssignedGene

pdf(here(plot_dir, "PCs_vs_QC_metrics.pdf"), width = 10)
walk(c("age", focused_qc_metrics), ~ {
    pca_v_con <- pca_long |>
        ggplot(aes(x = PC_val, y = .data[[.x]], shape = qc_class, color = library_combo)) +
        geom_point() +
        geom_text_repel(aes(label = ifelse(qc_class == "warn", Sample, "")), size = 2, color = "black") +
        facet_wrap(~PC_lab) +
        scale_color_manual(values = library_combo_colors) +
        labs(title = .x) +
        theme_bw()

    # ggsave(pca_v_con, filename = here(plot_dir, paste0("pca_v_",.x,".png")), width = 10)

    print(pca_v_con)
})
dev.off()

#### Check Br8325_mid ####

## Br8325_mid is an outlier in PC4
ggpairs(pca_tab,
                  mapping = aes(color = Sample == "Br8325_mid"),
                  columns = paste0("PC", 1:5),
                  upper = "blank"
)

pc1_pc4_Br8325_mid <- pca_tab |>
  ggplot(aes(x = PC1, y = PC4, color = Sample == "Br8325_mid", shape = qc_class)) +
  geom_point() +
  theme_bw() +
  # scale_color_viridis() +
  labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[4]])

ggsave(pc1_pc4_Br8325_mid, filename = here(plot_dir, "Bulk_PC1vPC4_Br8325_mid.png"))

pc4_Sample <- pca_tab |>
  ggplot(aes(x = Sample, y = PC4, color = Sample == "Br8325_mid")) +
  geom_boxplot() +
  theme_bw() +
  # scale_color_viridis() +
  labs(y = pca_vars_lab[[4]]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None")

ggsave(pc4_Sample, filename = here(plot_dir, "Bulk_PC4_boxplot_Br8325_mid.png"))



pca_tab |> filter(Sample == "Br8325_mid")

pd_long <-  pd_simple |>
  select(SAMPLE_ID, Sample, qc_class,  library_combo, age, all_of(focused_qc_metrics)) |>
  pivot_longer(!c(SAMPLE_ID, Sample, library_combo, qc_class))

qc_Sample <- pd_long |>
  ggplot(aes(Sample, value)) +
  geom_boxplot(aes(color = Sample == "Br8325_mid"), outlier.shape = NA)+
  geom_point(aes(shape = qc_class, color = library_combo)) +
  scale_color_manual(values = c(library_combo_colors, `TRUE` = "red")) +
  facet_wrap(~name, scales = "free_y")  +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(qc_Sample, filename = here(plot_dir, "Bulk_qc_Sample.png"), width = 11)

#### Check PC over round ####

pc_boxplot_round <- pca_long |>
  ggplot(aes(x = library_type, y = PC_val, color = round)) +
  geom_boxplot()+
  labs(y="") +
  facet_wrap(~PC_name, scales = "free_y")  +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(pc_boxplot_round, filename = here(plot_dir, "Bulk_PCA_round.png"), width = 11)


#### correlation heatmap ####
pd_long <-  pd_simple |>
  select(SAMPLE_ID, BrNum, library_combo, age, all_of(focused_qc_metrics)) |>
  pivot_longer(!c(SAMPLE_ID, BrNum, library_combo))

pca_longer <- pca$x[, 1:6] |>
  as.data.frame() |>
  rownames_to_column("SAMPLE_ID") |>
  pivot_longer(!SAMPLE_ID, names_to = "PC_name", values_to = "PC_val") |>
  right_join(pd_long) 


pca_cor <- pca_longer |>
  group_by(PC_name, name) |>
  summarise(cor = cor(PC_val, value))

pca_cor |>
  arrange(-cor)

pca_cor |>
  ggplot(aes(x = PC_name, y = name, fill = cor)) +
  geom_tile() +
  geom_text(aes(label = round(cor, 3))) +
  scale_fill_gradient2(
    low = "#0571B0", 
    mid = "white", 
    high = "#CA0020", 
    midpoint = 0
  )

## use vairancePartition
library(variancePartition)

C = canCorPairs(~BrNum + Position + library_type + library_prep  + sex + age + 
                  concordMapRate + mitoRate + numMapped + overallMapRate + totalAssignedGene + totalMapped +
                  PC1 + PC2 + PC3 + PC4 + PC5, 
                pca_tab)

C[!grepl("PC", rownames(C)),grepl("PC", colnames(C))]

# Plot correlation matrix
pdf(here(plot_dir, "PCA_canCorPairs.pdf"))
plotCorrMatrix(C)
dev.off()

# sgejobs::job_single('03_qc_pca', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript 03_qc_pca.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
