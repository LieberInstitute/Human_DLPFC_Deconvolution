library("SummarizedExperiment")
library("tidyverse")
library("sessioninfo")
library("here")
library("readxl")
library("ggrepel")
library("jaffelab")
library("scuttle")
library("GGally")
library("patchwork")

## prep dirs ##
plot_dir <- here("plots", "02_quality_control", "01_prelim_bulk_qc_check")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load Big Data ####
load("/dcl01/ajaffe/data/lab/brain_swap/count_data/RNAseq_Collection_postQC_n5761_12dataset_2021-06_geneRSE.Rdata", verbose = TRUE)
dim(rse_gene)
# [1] 58037  5761

bsp2_library_type <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/misc/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx') |>
  mutate(RNum = paste0("R", RNum)) |>
  select(RNum, library_type = Protocol)

table(bsp2_library_type$library_type)
# RiboZeroGold  RiboZeroHMR 
# 47          405 


pd_big <- as.data.frame(colData(rse_gene)) |>
  left_join(bsp2_library_type) |>
  replace_na(list(library_type = "RiboZeroGold")) |>
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

## historic cutoffs

historic_cutoffs <- pd_big |>
  group_by(library_type, Dataset) |>
  summarise(concordMapRate_min = min(concordMapRate),
            mitoMapped_max = max(mitoMapped),
            mitoRate_max = max(mitoRate),
            numMapped_min = min(numMapped),
            numReads_min = min(numReads),
            overallMapRate_min = min(overallMapRate),
            rRNArate_max = max(rRNA_rate),
            totalAssignedGene_min = min(totalAssignedGene),
            totalMapped_min = min(totalMapped)
            )

historic_cutoff_long  <- historic_cutoffs |>
  pivot_longer(!c(library_type, Dataset), names_to = "cutoff")

historic_cutoff_box  <- historic_cutoff_long |>
  ggplot(aes(x = library_type, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_text_repel(aes(label = Dataset), color = "blue", size = 1)+
  facet_wrap(~cutoff, scales = "free_y")

ggsave(historic_cutoff_box, filename = here(plot_dir, "historic_cutoff_boxplot.png"))

## sanity check numbers
historic_cutoff_median <- historic_cutoff_long |>
  filter(library_type == "RiboZeroGold") |>
  group_by(cutoff) |>
  summarise(Median = median(value)) |>
  separate(cutoff, into = c("qc_var", "type"), sep = "_") |>
  mutate(qc_var = gsub("rRNArate", "rRNA_rate", qc_var))

# qc_var            type   Median
# <chr>             <chr>   <dbl>
# 1 concordMapRate    min   5.22e-1
# 2 mitoMapped        max   1.19e+7
# 3 mitoRate          max   9.96e-2
# 4 numMapped         min   4.93e+7
# 5 numReads          min   6.75e+7
# 6 overallMapRate    min   5.47e-1
# 7 rRNA_rate         max   1.16e-5
# 8 totalAssignedGene min   2.94e-1
# 9 totalMapped       min   4.66e+7

#### Load New Data ####
load(here("processed-data","rse","preQC","rse_gene_preQC.Rdata"), verbose = TRUE)
pd_new <- as.data.frame(colData(rse_gene)) |>
  mutate(Dataset = seq_set) ## eval as dataset here - but keep as seq_set in rse_gene

pd_new |>
  count(Dataset, round, library_type)
#         Dataset library_type  n
# 1 2107UNHS-0291        polyA 12
# 2 2107UNHS-0293 RiboZeroGold 12
# 3    AN00000904        polyA 44
# 4    AN00000906 RiboZeroGold 45

## pass fail matrix
pass_fail <- pd_new |>
  select(SAMPLE_ID,library_type, Dataset, basic_statistics:adapter_content) |>
  pivot_longer(!c(SAMPLE_ID,library_type,Dataset), names_to = "Test")

pass_fail |> group_by(Test) |>
  count(value)


pass_fail |> group_by(Dataset, Test) |>
  count(value) |>
  filter(value=="FAIL")

# Dataset       Test                        value     n
# <chr>         <chr>                       <chr> <int>
# 1 2107UNHS-0291 per_sequence_gc_content     FAIL      1
# 2 2107UNHS-0291 sequence_duplication_levels FAIL     12
# 3 2107UNHS-0293 per_base_sequence_content   FAIL      1
# 4 2107UNHS-0293 per_sequence_gc_content     FAIL     10
# 5 2107UNHS-0293 sequence_duplication_levels FAIL      4
# 6 AN00000904    per_sequence_gc_content     FAIL      4
# 7 AN00000904    sequence_duplication_levels FAIL     44
# 8 AN00000906    per_base_sequence_content   FAIL      4
# 9 AN00000906    per_sequence_gc_content     FAIL     41
# 10 AN00000906    sequence_duplication_levels FAIL     18


pass_fail_plot <- ggplot(pass_fail, aes(x = Test, y = SAMPLE_ID)) +
  geom_tile(aes(fill = value)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~library_type+Dataset, scales = "free_y")

ggsave(pass_fail_plot, filename = here(plot_dir,"pass_fail_tile.png"), height = 10, width = 10)

## Subset and QC data
qc_variables <- c("numReads", "numMapped", "numUnmapped", "overallMapRate", "concordMapRate", "totalMapped", "mitoMapped","mitoRate", "rRNA_rate", "totalAssignedGene")

pd_big_qc_long <- pd_big |>
  select(SAMPLE_ID, Dataset, library_type, all_of(qc_variables)) |>
  pivot_longer(!c(SAMPLE_ID, Dataset, library_type), names_to = "qc_var") |>
  mutate(`Data Status` = "old")


pd_new_qc_long <- pd_new |>
  select(SAMPLE_ID, Dataset, library_type, library_prep, all_of(qc_variables)) |>
  pivot_longer(!c(SAMPLE_ID, Dataset, library_type, library_prep), names_to = "qc_var") |>
  mutate(`Data Status` = "new")

## compare other LIBD datasets
pd_qc_long <- rbind(pd_new_qc_long |> select(-library_prep), 
                    pd_big_qc_long)

qc_boxplots_LIBD_exp <- ggplot(data = pd_qc_long, aes(x = Dataset, y = value)) +
  geom_boxplot(aes(fill = `library_type`, color = `Data Status`)) +
  facet_wrap(~qc_var, scales = "free_y", ncol = 5 ) +
  scale_color_manual(values = c(new = "black", old = "grey50"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 17)) +
  geom_hline(data = historic_cutoff_median, aes(yintercept = Median), color = "darkgreen")

ggsave(qc_boxplots_LIBD_exp , filename = here(plot_dir, "qc_boxplots_LIBD_experiment.png"), width = 24, height = 12)


qc_boxplots_LIBD_lt <- ggplot(data = pd_qc_long, aes(x = library_type, y = value)) +
  geom_boxplot(aes(fill = `library_type`, color = `Data Status`)) +
  facet_wrap(~qc_var, scales = "free_y", ncol = 5 ) +
  scale_color_manual(values = c(new = "black", old = "grey50"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 17))+
  geom_hline(data = historic_cutoff_median, aes(yintercept = Median), color = "darkgreen")


ggsave(qc_boxplots_LIBD_lt, filename = here(plot_dir, "qc_boxplots_LIBD_library-type.png"), width = 24, height = 12)


qc_boxplots_LIBD_lt2 <- ggplot(data = pd_qc_long, aes(x = library_type, y = value)) +
  geom_boxplot(aes(fill = `library_type`, color = `Data Status`)) +
  facet_wrap(~qc_var, scales = "free_y", ncol = 5 ) +
  scale_color_manual(values = c(new = "black", old = "grey50"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 17))+
  geom_hline(data = historic_cutoff_median, aes(yintercept = Median), color = "darkgreen")


ggsave(qc_boxplots_LIBD_lt, filename = here(plot_dir, "qc_boxplots_LIBD_library-type.png"), width = 24, height = 12)


#### Combined Data automatic outliers ####
library(scuttle)

focused_qc_metrics <- c("numReads",
                        "numMapped",
                        "overallMapRate", ## high cor with concordMapRate, -numUnmapped
                        "totalAssignedGene",
                        "mitoRate", # cor with mitoMapped
                        "rRNA_rate")

names(focused_qc_metrics) <- focused_qc_metrics
tail <- c(numReads = "lower", 
          numMapped = "lower",
          overallMapRate = "lower",
          totalAssignedGene = "lower",
          mitoRate = "higher",
          rRNA_rate = "higher")

## apply isOutlier
auto_outliers <- map2(focused_qc_metrics, tail, ~isOutlier(rse_gene[[.x]], type = .y, nmads = 4))

## find MAD cutoffs
auto_cutoff <- map_dfr(auto_outliers, ~attr(.x,"thresholds")) |>
  add_column(qc_var = focused_qc_metrics, .before = 1)

map(auto_outlier, table)

auto_outliers_combine <- Reduce("|", auto_outliers)
table(auto_outliers_combine)
# FALSE  TRUE 
# 85    28 

auto_outlier_tb <- tibble(SAMPLE_ID = pd_new$SAMPLE_ID, auto_outlier = auto_outliers_combine)

### boxplot for global isOutlier cutoffs
qc_boxplot_isOutlier <- pd_new_qc_long |>
  left_join(auto_outlier_tb) |>
  filter(qc_var %in% focused_qc_metrics) |>
  ggplot(aes(x = library_type, y = value)) +
      geom_boxplot(aes(fill = library_type), outlier.shape = NA, alpha = 0.5) +
      geom_jitter(aes(color = auto_outlier)) +
      facet_wrap(~qc_var, scales = "free_y", nrow = 2) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "isOutlier cutoff - 4 MADs") +
      scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "grey50")) +
      geom_hline(data = auto_cutoff, aes(yintercept = lower), color = "red", linetype = "dashed") +
      geom_hline(data = auto_cutoff, aes(yintercept = higher), color = "blue", linetype = "dashed")+
  geom_hline(data = mdd_cutoffs, aes(yintercept = cutoff), color = "purple", linetype = "dotted")
#     
ggsave(qc_boxplot_isOutlier, filename = here(plot_dir, paste0("qc_boxplots_isOutlier.png")), width = 12)

## by library prep?

#### New Data Only Boxplots ####
mdd_cutoffs <- tibble(qc_var = c("overallMapRate","totalAssignedGene","numReads","rRNA_rate"),
                      cutoff = c(0.5, 0.3, 10^7.25, 1e-3)
                  )


qc_boxplots_lt <- pd_new_qc_long |>
  group_by(library_type) |>
  group_map(
    ~{qc_boxplot <- ggplot(.x, aes(x = Dataset, y = value)) +
      geom_boxplot(outlier.shape = NULL) +
      geom_jitter(width = .2) +
      facet_wrap(~qc_var, scales = "free_y", nrow = 2) +
      scale_color_manual(values = c(new = "black", old = "grey50"))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = .y) +
      geom_hline(data = mdd_cutoffs, aes(yintercept = cutoff), color = "red", linetype = "dashed")
    
    ggsave(qc_boxplot, filename = here(plot_dir, paste0("qc_boxplots_",.y,".png")), width = 12)
    return(qc_boxplot)
    }
    )

ggsave(qc_boxplots_lt[[1]]/ qc_boxplots_lt[[2]], filename = here(plot_dir, "qc_boxplots.png"), width = 24, height = 12)

## color by library_prep
pd_new_qc_long |>
  group_by(library_type) |>
  group_map(
    ~{qc_boxplot <- ggplot(.x, aes(x = Dataset, y = value)) +
      geom_boxplot(aes(fill = library_prep), outlier.shape = NA, alpha = 0.5) +
      geom_point(aes(color = library_prep), position = position_jitterdodge()) +
      facet_wrap(~qc_var, scales = "free_y", nrow = 2) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = .y) +
      geom_hline(data = mdd_cutoffs, aes(yintercept = cutoff), color = "red", linetype = "dashed")
    
    ggsave(qc_boxplot, filename = here(plot_dir, paste0("qc_boxplots_",.y,".png")), width = 12)
    # return(qc_boxplot)
    }
  )

#### gg pairs ####
pd_new |> 
  group_by(library_type) |>
  group_map(~{
    
    qc_ggpair <- ggpairs(.x, 
                         columns = qc_variables,
                         ggplot2::aes(colour=library_prep, alpha = .6)) + 
      theme_bw() +
      labs(title = .y)
    
    ggsave(qc_ggpair, filename = here(plot_dir, paste0("qc_ggpairs_",.y,"_all.png")), width = 14, height = 14) 
  })


pd_new |>
  filter(numMapped < 5e7) |>
  select(SAMPLE_ID,library_type, numMapped)

pd_new |>
  filter(library_type == "polyA",
         totalAssignedGene < 0.3) |>
  select(SAMPLE_ID, totalAssignedGene)

pd_new |>
  filter(library_type == "RiboZeroGold",
         totalAssignedGene < 0.2) |>
  select(SAMPLE_ID, totalAssignedGene)

#### Hand picked samples 
bad_samples <- c("2107UNHS-0293_Br2720_Mid_Nuc", #low mapped
                 "AN00000904_Br2743_Ant_Cyto", ## low mapped, low totalAssignedGene
                 "AN00000906_Br2743_Ant_Cyto",
                 "AN00000906_Br2743_Ant",
                 "AN00000904_Br2743_Ant",
                 "AN00000906_Br8325_Mid_Nuc",
                 "AN00000906_Br8492_Mid_Nuc"
                 )

pd_new |> 
  group_by(library_type) |>
  mutate(drop = SAMPLE_ID %in% bad_samples) |>
  count(drop)

pd_new |> 
  group_by(library_type) |>
  mutate(drop = SAMPLE_ID %in% bad_samples) |>
  group_map(~{
    
    qc_ggpair <- ggpairs(.x, 
                         columns = qc_variables,
                         ggplot2::aes(colour=drop, alpha = .6)) + 
      theme_bw() +
      labs(title = .y)
    
    ggsave(qc_ggpair, filename = here(plot_dir, paste0("qc_ggpairs_",.y,"_drop.png")), width = 14, height = 14) 
  })




pd_new |> 
  group_by(library_type) |>
  group_map(~{
    
    qc_ggpair <- ggpairs(.x, 
                         columns = focused_qc_metrics,
                         ggplot2::aes(colour=library_prep, alpha = .6)) + 
      theme_bw() +
      labs(title = .y)
    
    ggsave(qc_ggpair, filename = here(plot_dir, paste0("qc_ggpairs_",.y,"_focus.png")), width = 10, height = 10) 
  })


#### ERCC data ####
ercc_boxplot <- pd_new |>
  ggplot(aes(x = Dataset, y = ERCCsumLogErr, fill = library_type)) +
  geom_boxplot(outlier.shape = NULL) +
  geom_jitter(width = .2) +
  facet_wrap(~round, scales = "free")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(subtitle = "No ERCC spike in for Round 1")

ggsave(ercc_boxplot, filename = here(plot_dir, "ERCCsumLogErr_boxplot.png"))

## there was no ERCC spike in for round1 so those values are not 
pd_new <- pd_new |>
  mutate(ERCCsumLogErr = ifelse(round == 1, NA, ERCCsumLogErr))
  

## no RIN?
"RIN" %in% colnames(pd_new)

#### Metrics vs. ERCC ###
ercc_check <- pd_new |>
  select(SAMPLE_ID, ERCCsumLogErr, library_prep, round) |>
  # filter(!is.na(ERCCsumLogErr)) |>
  replace_na(list(ERCCsumLogErr = Inf)) |>
  left_join(pd_new_qc_long |>
              filter(qc_var %in% focused_qc_metrics),
            multiple = "all")
  

ercc_check |>
  group_by(library_type) |>
  group_map(~{
    
    ercc_scatter <- ggplot(.x, aes(x = ERCCsumLogErr, y = value, color = library_prep, shape = factor(round))) +
      geom_point() +
      facet_wrap(~qc_var, scales = "free_y", nrow = 1) +
      labs(title = .y)+
      geom_hline(data = mdd_cutoffs, aes(yintercept = cutoff), color = "red", linetype = "dashed")
    
    ggsave(ercc_scatter, filename = here(plot_dir, paste0("ERCC_scatter_", .y,".png")), width = 14, height = 5)
    
  })
  


# sgejobs::job_single('01_prelim_bulk_qc_check', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript 01_prelim_bulk_qc_check.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

