library(SummarizedExperiment)
library(tidyverse)
library(sessioninfo)
library(here)

#### Load Big Data ####
load("/dcl01/ajaffe/data/lab/brain_swap/count_data/RNAseq_Collection_postQC_n5761_12dataset_2021-06_geneRSE.Rdata", verbose = TRUE)
dim(rse_gene)
# [1] 58037  5761

pd <- as.data.frame(colData(rse_gene))

table(pd$Dataset)
# Astellas_DG         BrainSeq_Phase1   BrainSeq_Phase2_DLPFC   BrainSeq_Phase2_HIPPO BrainSeq_Phase3_Caudate 
#         263                     727                     453                     447                     464 
# BrainSeq_Phase4and5                Habenula            Nicotine_NAc          psychENCODE_BP        psychENCODE_Mood 
#                 490                      69                     235                      17                    1091 
# PTSD_BrainOmics                 VA_PTSD 
#             225                    1280 

colnames(pd)
# [1] "SAMPLE_ID"         "RNum"              "RIN"               "Region"            "Dataset"          
# [6] "BrNum"             "Dx"                "Age"               "Sex"               "Race"             
# [11] "numReads"          "numMapped"         "numUnmapped"       "mitoMapped"        "totalMapped"      
# [16] "overallMapRate"    "concordMapRate"    "mitoRate"          "rRNA_rate"         "totalAssignedGene"
# [21] "bamFile"           "rna_preSwap_BrNum"

## Subset and QC data
qc_variables <- c("numReads", "overallMapRate","mitoRate", "rRNA_rate", "totalAssignedGene")

pd_qc_long <- pd %>%
  select(SAMPLE_ID, Dataset, all_of(qc_variables)) %>%
  pivot_longer(!c(SAMPLE_ID, Dataset), names_to = "qc_var") %>%
  mutate(`Data Status` = "old")

qc_boxplots <- ggplot(data = pd_qc_long, aes(x = Dataset, y = value)) +
  geom_boxplot(aes(fill = `Data Status`)) +
  facet_wrap(~qc_var, scales = "free_y", ncol = 2 ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15))

ggsave(qc_boxplots, filename = here("plots", "quality_control", "qc_boxplots.png"), width = 10, height = 12)

  
