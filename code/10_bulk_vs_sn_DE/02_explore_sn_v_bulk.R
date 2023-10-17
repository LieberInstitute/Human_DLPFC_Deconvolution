
library("SummarizedExperiment")
library("tidyverse")
library("sessioninfo")
library("here")
library("jaffelab")
library("recount")
# library("viridis")
library("ggrepel")
library("GGally")
# library("vsn")

## prep dirs ##
plot_dir <- here("plots", "10_bulk_vs_sn_DE", "02_explore_sn_v_bulk")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data", "00_data_prep", "bulk_colors.Rdata"), verbose = TRUE)
# library_combo_colors
# library_prep_colors
# library_type_colors

library_combo_colors <- c(library_combo_colors, `snRNA-seq` = "#417B5A")

#### Load Data ####
## bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
# rse_gene

## pseudobulk sce
load(here("processed-data", "10_bulk_vs_sn_DE", "sce_pb_sample.Rdata"), verbose = TRUE)
# sce_pb_sample

#### create combined colData ####
## rm all NA cols
colData(sce_pb_sample) <- colData(sce_pb_sample)[,!apply(apply(colData(sce_pb_sample), 2, is.na), 2, all)]
colData(sce_pb_sample)$data_type <- "snRNA-seq"
colData(sce_pb_sample)$library_combo <- "snRNA-seq"
colData(sce_pb_sample)$library_type <- "snRNA-seq"
colData(sce_pb_sample)$library_prep <- "snRNA-seq"

colnames(colData(sce_pb_sample))
# [1] "Sample"                 "SAMPLE_ID"              "pos"                    "BrNum"                 
# [5] "round"                  "Position"               "age"                    "sex"                   
# [9] "diagnosis"              "high_mito"              "low_sum"                "low_detected"          
# [13] "discard_auto"           "registration_variable"  "registration_sample_id" "ncells"                
# [17] "data_type" 

colData(rse_gene)$data_type <- "Bulk"

(common_colnames <- intersect(colnames(colData(sce_pb_sample)), colnames(colData(rse_gene))))
# [1] "Sample"        "SAMPLE_ID"     "pos"           "BrNum"         "round"         "Position"      "age"          
# [8] "sex"           "diagnosis"     "library_combo" "library_type"  "library_prep" 

colData(sce_pb_sample) <- colData(sce_pb_sample)[,common_colnames]
colData(rse_gene) <- colData(rse_gene)[,common_colnames]


#### Match rowData ####
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

common_genes <- intersect(rownames(sce_pb_sample), rownames(rse_gene))
length(common_genes)
# [1] 17660

## subset to common genes
rse_gene <- rse_gene[common_genes,]
sce_pb_sample <- sce_pb_sample[common_genes,]

rse_combine <- SummarizedExperiment(colData = rbind(colData(rse_gene), 
                                                    colData(sce_pb_sample)),
                                    rowData=rowData(rse_gene),
                                    assays = list(counts = cbind(assays(rse_gene)$counts, 
                                                                 assays(sce_pb_sample)$counts),
                                                  logcounts = cbind(assays(rse_gene)$logcounts, 
                                                                 assays(sce_pb_sample)$logcounts)
                                                  )
                                    )

# move to colData
pd <- as.data.frame(colData(rse_combine))

#### filter genes by Expression ####
gene_rpkm <- getRPKM(rse_combine, "Length")
rse_gene_filter <- rse_gene[rowMeans(gene_rpkm) > 0.1, ]
gene_rpkm_filter <- gene_rpkm[rowMeans(gene_rpkm) > 0.1, ]
table(droplevels(seqnames(rse_gene_filter)))
# 
# pdf(here(plot_dir, "mean_vs_sd_rpkm.pdf"))
# meanSdPlot(gene_rpkm, ranks = FALSE)
# dev.off()
# 
# ## get expression
geneExprs_filter <- log2(gene_rpkm_filter + 1)
# assays(rse_gene_filter)$logcounts <- log2(gene_rpkm_filter + 1) ## check


# pdf(here(plot_dir, "mean_vs_sd_logrpkm.pdf"))
# meanSdPlot(geneExprs_filter, ranks = FALSE)
# dev.off()

## calc PCA and vars
pca <- prcomp(t(geneExprs_filter))
pca_vars <- getPcaVars(pca)
pca_vars_lab <- paste0(
    "PC", seq(along = pca_vars), ": ",
    pca_vars, "% Var Expl"
)

pca_vars

pca$x[1:5,1:5]

## calc PCA and vars
pca2 <- prcomp(t(assays(rse_combine)$logcounts))
pca_vars2 <- getPcaVars(pca2)
pca_vars_lab2 <- paste0(
  "PC", seq(along = pca_vars2), ": ",
  pca_vars, "% Var Expl"
)


pca2$x[1:5,1:5]

## create table with groups and some QC metrics
pca_tab <- pd |> cbind(pca$x[, 1:5])
pca_tab2 <- pd |> cbind(pca2$x[, 1:5])

#### ggpairs for pca ####
gg_pca <- ggpairs(pca_tab,
    mapping = aes(color = library_combo),
    columns = paste0("PC", 1:5),
    upper = "blank"
) +
    scale_color_manual(values = library_combo_colors) +
    scale_fill_manual(values = library_combo_colors)

ggsave(gg_pca, filename = here(plot_dir, "sn_vs_Bulk_ggpairs_pca.png"), height = 10, width = 10)

gg_pca2 <- ggpairs(pca_tab2,
                  mapping = aes(color = library_combo),
                  columns = paste0("PC", 1:5),
                  upper = "blank"
) +
  scale_color_manual(values = library_combo_colors) +
  scale_fill_manual(values = library_combo_colors)

ggsave(gg_pca2, filename = here(plot_dir, "sn_vs_Bulk_ggpairs_pca2.png"), height = 10, width = 10)


gg_pca_position <- ggpairs(pca_tab,
    mapping = aes(color = Position),
    columns = paste0("PC", 1:5),
    upper = "blank"
)

ggsave(gg_pca_position, filename = here(plot_dir, "sn_vs_Bulk_ggpairs_position.png"), height = 10, width = 10)

## PC1 vs. PC2
pc_test <- pca_tab |>
    ggplot(aes(x = PC1, y = PC2, color = library_combo)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = library_combo_colors) +
    labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[2]]) 
    # coord_equal()

ggsave(pc_test, filename = here(plot_dir, "sn_vs_Bulk_PC1vPC2_library_combo.png"))

pc1v2_lab <- pca_tab |>
    ggplot(aes(x = PC1, y = PC2, color = library_combo)) +
    geom_point() +
    geom_text_repel(aes(label = Sample), color = "black", size = 2) +
    theme_bw() +
    scale_color_manual(values = library_combo_colors) +
    labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[2]]) +
    coord_equal()

ggsave(pc1v2_lab, filename = here(plot_dir, "sn_vs_Bulk_PC1vPC2_library_combo_label.png"))

pc_test <- pca_tab |>
    ggplot(aes(x = PC1, y = PC2, color = BrNum)) +
    geom_point() +
    theme_bw() +
    # scale_color_manual(values = library_combo_colors) +
    labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[2]])

ggsave(pc_test, filename = here(plot_dir, "sn_vs_Bulk_PC1vPC2_BrNum.png"))

pc_test <- pca_tab |>
    ggplot(aes(x = PC1, y = PC2, color = mitoRate, shape = qc_class)) +
    geom_point() +
    theme_bw() +
    scale_color_viridis() +
    labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[2]])

ggsave(pc_test, filename = here(plot_dir, "sn_vs_Bulk_PC1vPC2_mitoRate.png"))


# sgejobs::job_single('03_qc_pca', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript 03_qc_pca.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
