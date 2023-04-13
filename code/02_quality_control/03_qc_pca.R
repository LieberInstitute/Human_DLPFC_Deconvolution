library("SummarizedExperiment")
library("tidyverse")
library("sessioninfo")
library("here")
library("jaffelab")
library("recount")
library("viridis")

## prep dirs ##
plot_dir <- here("plots", "02_quality_control", "03_qc_pca")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data","00_data_prep","bulk_colors.Rdata"), verbose = TRUE)

#### Load Data ####
load(here("processed-data","rse","preQC","rse_gene_preQC.Rdata"), verbose = TRUE)
pd <- as.data.frame(colData(rse_gene))

## Add qc record and drop samples
qc_tb <- read.csv(file = here("processed-data", "02_quality_control","QC_record_DLPFC_bulk.csv"))
qc_tb |> filter(qc_class != "pass")

identical(qc_tb$SAMPLE_ID,colnames(rse_gene))
rse_gene$qc_class <- qc_tb$qc_class 
table(rse_gene$qc_class)
# drop pass warn 
# 2  101   10 

## drop 2 poor QC samples
rse_gene <- rse_gene[,rse_gene$qc_class != "drop"]

#move to colData
# rse_gene$library_combo <- paste0(rse_gene$library_type,"_",rse_gene$library_prep)
# table(rse_gene$library_combo)

pd <- as.data.frame(colData(rse_gene))

# qc_variables <- c("numReads", "numMapped", "numUnmapped", "overallMapRate", "concordMapRate", "totalMapped", "mitoMapped","mitoRate", "rRNA_rate", "totalAssignedGene")
pd_simple <- pd |> select(SAMPLE_ID, Sample, BrNum, Position, library_type, library_prep, library_combo, qc_class, mitoRate,totalAssignedGene)

#### filter genes by Expression ####
gene_rpkm <- getRPKM(rse_gene,"Length")
rse_gene_filter <- rse_gene[rowMeans(gene_rpkm) > 0.1,]
gene_rpkm_filter <- gene_rpkm[rowMeans(gene_rpkm) > 0.1,]
table(droplevels(seqnames(rse_gene_filter)))

## get expression
geneExprs_filter <- log2(gene_rpkm_filter+1)
assays(rse_gene_filter)$logcounts <- log2(gene_rpkm_filter+1) ## check 

pca = prcomp(t(geneExprs_filter))
pca_vars = getPcaVars(pca)
pca_vars_lab = paste0("PC", seq(along=pca_vars), ": ",
                      pca_vars, "% Var Expl")

names(pca_vars)
names(pca)

pca$x[,1:5]

pca_tab <- pd_simple |> cbind(pca$x[,1:5])

# missing bulk points?
pc_test <- pca_tab |>
  ggplot(aes(x = PC1, y = PC2, color = library_combo, shape = qc_class)) +
  geom_point()+
  theme_bw() +
  scale_color_manual(values = library_combo_colors) +
  labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[2]]) 

ggsave(pc_test, filename = here(plot_dir, "Bulk_PC1vPC2_library_combo.png"))


pc_test <- pca_tab |>
  ggplot(aes(x = PC1, y = PC2, color = BrNum, shape = qc_class)) +
  geom_point()+
  theme_bw() +
  # scale_color_manual(values = library_combo_colors) +
  labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[2]]) 

ggsave(pc_test, filename = here(plot_dir, "Bulk_PC1vPC2_BrNum.png"))

pc_test <- pca_tab |>
  ggplot(aes(x = PC1, y = PC2, color = mitoRate, shape = qc_class)) +
  geom_point()+
  theme_bw() +
  scale_color_viridis() +
  labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[2]]) 

ggsave(pc_test, filename = here(plot_dir, "Bulk_PC1vPC2_mitoRate.png"))



