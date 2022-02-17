library(SummarizedExperiment)
library(tidyverse)
library(sessioninfo)
library(here)
library(readxl)
library(ggrepel)
library(jaffelab)
library(plotly)

#### Load Big Data ####
load("/dcl01/ajaffe/data/lab/brain_swap/count_data/RNAseq_Collection_postQC_n5761_12dataset_2021-06_geneRSE.Rdata", verbose = TRUE)
dim(rse_gene)
# [1] 58037  5761

bsp2_library_type <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/misc/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx') %>%
  mutate(RNum = paste0("R", RNum)) %>%
  select(RNum, library_type = Protocol)

table(bsp2_library_type$library_type)
# RiboZeroGold  RiboZeroHMR 
# 47          405 


pd_big <- as.data.frame(colData(rse_gene)) %>%
  left_join(bsp2_library_type) %>%
  replace_na(list(library_type = "RiboZeroGold")) %>%
  mutate(library_type = ifelse(Dataset == "BrainSeq_Phase1", "polyA", library_type))

table(pd_big$library_type)
# PolyA RiboZeroGold  RiboZeroHMR 
# 727         4630          404 

table(pd_big$Dataset)
# Astellas_DG         BrainSeq_Phase1   BrainSeq_Phase2_DLPFC   BrainSeq_Phase2_HIPPO BrainSeq_Phase3_Caudate 
#         263                     727                     453                     447                     464 
# BrainSeq_Phase4and5                Habenula            Nicotine_NAc          psychENCODE_BP        psychENCODE_Mood 
#                 490                      69                     235                      17                    1091 
# PTSD_BrainOmics                 VA_PTSD 
#             225                    1280 

colnames(pd_big)
# [1] "SAMPLE_ID"         "RNum"              "RIN"               "Region"            "Dataset"          
# [6] "BrNum"             "Dx"                "Age"               "Sex"               "Race"             
# [11] "numReads"          "numMapped"         "numUnmapped"       "mitoMapped"        "totalMapped"      
# [16] "overallMapRate"    "concordMapRate"    "mitoRate"          "rRNA_rate"         "totalAssignedGene"
# [21] "bamFile"           "rna_preSwap_BrNum" "library_type" 

#### Load New Data ####
load(here("processed-data","01_SPEAQeasy","round1_2021-10-19","rse_gene.Rdata"), verbose = TRUE)
pd_new <- as.data.frame(colData(rse_gene))

## pass fail matrix
pass_fail <- pd_new %>%
  select(SAMPLE_ID, basic_statistics:adapter_content) %>%
  pivot_longer(!SAMPLE_ID, names_to = "Test")

pass_fail %>% group_by(Test) %>%
  count(value)

pass_fail_plot <- ggplot(pass_fail, aes(x = Test, y = SAMPLE_ID)) +
  geom_tile(aes(fill = value)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(pass_fail_plot, filename = here("plots","02_quality_control","pass_fail_tile.png"))

## Subset and QC data
qc_variables <- c("numReads", "numMapped", "numUnmapped", "overallMapRate", "concordMapRate", "totalMapped", "mitoMapped","mitoRate", "rRNA_rate", "totalAssignedGene")

pd_big_qc_long <- pd_big %>%
  select(SAMPLE_ID, Dataset, library_type, all_of(qc_variables)) %>%
  pivot_longer(!c(SAMPLE_ID, Dataset, library_type), names_to = "qc_var") %>%
  mutate(`Data Status` = "old")


pd_new_qc_long <- pd_new %>%
  select(SAMPLE_ID, Dataset, library_type, all_of(qc_variables)) %>%
  pivot_longer(!c(SAMPLE_ID, Dataset, library_type), names_to = "qc_var") %>%
  mutate(`Data Status` = "new")

pd_qc_long <- rbind(pd_new_qc_long, pd_big_qc_long)

qc_boxplots <- ggplot(data = pd_qc_long, aes(x = Dataset, y = value)) +
  geom_boxplot(aes(fill = `library_type`)) +
  facet_wrap(~qc_var, scales = "free_y", ncol = 5 ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 17))

ggsave(qc_boxplots, filename = here("plots", "02_quality_control", "qc_boxplots.png"), width = 24, height = 12)


#### Just New data plotly boxplots###
pd_new_qc_key <- pd_new_qc_long %>%
  mutate(anno = ss(SAMPLE_ID, "-", 3)) %>%
  highlight_key(~anno)

qc_new_boxplots <-  pd_new_qc_key %>%
  ggplot(aes(x = Dataset, y = value)) +
  geom_boxplot(aes(fill = `library_type`), outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(~qc_var, scales = "free_y", ncol = 5 ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 17)) +
  labs(title = "Human DLPFC Deconvolution", subtitle = "Bulk Round 1")

ggsave(qc_new_boxplots, filename = here("plots", "02_quality_control", "qc_new_boxplots.png"), width = 20, height = 12)

qc_new_boxplots_plotly <- ggplotly(qc_new_boxplots, tooltip = c("colour","text"))

htmlwidgets::saveWidget(highlight(qc_new_boxplots_plotly,
                                  on = "plotly_click",
                                  off = "plotly_doubleclick",
                                  selectize = TRUE,
                                  dynamic = TRUE,
                                  persistent = FALSE,),
                        selfcontained = FALSE,
                        file = here("plots", "02_quality_control", "qc_new_boxplots.html"))


pd_new %>% filter(totalMapped < 5e7)
