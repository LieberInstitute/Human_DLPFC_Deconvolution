
library("tidyverse")
library("SummarizedExperiment")
library("here")
library("sessioninfo")

#### Plot Setup ####
plot_dir = here("plots","00_data_prep","02_update_experiment_tile")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

pos_df <- tibble(Position = c("Anterior", "Middle", "Posterior"), 
                 pos = c("ant", "mid", "post"))

## load colors
load(here("processed-data", "00_data_prep", "bulk_colors.Rdata"), verbose = TRUE)

#### snRNA-seq ####

load(here("processed-data", "sce", "sce_pd.Rdata"), verbose = TRUE)

sn_samples <- sce_pd |>
  dplyr::count(Sample, SAMPLE_ID, Position, pos, round, BrNum, age, sex) 

## should be the same as before (QC happened previously)
sn_n_samp <- sn_samples |>
  dplyr::count(Sample, Position, BrNum) |>
  dplyr::mutate(data_type = "snRNA-seq")

sn_demo <- sn_samples |> group_by(BrNum, age, sex) |> summarize(n_snRNA = n())

#### bulk data ####

## post QC bulk data (dropped 3 samples)
load(here("processed-data","rse","rse_gene.Rdata"), verbose = TRUE)

bulk_samples <- colData(rse_gene) |>
  as.data.frame() |>
  select(Sample, SAMPLE_ID, Position, pos, library_type, rna_extract, library_combo2, round, BrNum, age, sex) |>
  mutate(rna_extract = factor(rna_extract, levels = c("Nuc", "Total", "Cyto")))

bulk_n_samp <- bulk_samples |>
  dplyr::count(Sample, Position, BrNum, library_type) |>
  mutate(data_type = paste0("Bulk RNA-seq\n", library_type)) |>
  select(-library_type)

bulk_demo <- bulk_samples |> group_by(BrNum, age, sex) |> summarize(n_bulk = n())

## type x prep
library_prep_type_count <- bulk_samples |>
  dplyr::count(rna_extract, library_type, library_combo2)

bulk_tile <- library_prep_type_count |>
  mutate(library_type = gsub("RiboZeroGold", "RiboZero\nGold", library_type))|>
  ggplot(aes(rna_extract, library_type, fill = library_combo2)) +
  geom_tile() +
  geom_text(aes(label = n), color = "white", fontface = "bold") +
  scale_fill_manual(values = library_combo_colors2) +
  theme_bw() + 
  theme(legend.position='none') +
  labs(x = "RNA Extraction", y = "Library Type") +
  coord_flip()
  
ggsave(bulk_tile, filename = here(plot_dir, "bulk_sample_tile.png"), height = 3, width = 2)
ggsave(bulk_tile, filename = here(plot_dir, "bulk_sample_tile.pdf"), height = 3, width = 2)
 
#### HALO Info ####
halo_info <- read.csv(here("processed-data","03_HALO", "01_import_HALO_data","HALO_metadata.csv"))
head(halo_info)

halo_info |> dplyr::count(Combo, Confidence)

halo_samp <- halo_info |>
  filter(Confidence %in% c("High", "OK"))|> 
  select(SAMPLE_ID, BrNum, Position, Combo) |>
  left_join(pos_df) |>
  mutate(Sample = paste0(BrNum, "_", pos))

halo_n_samp <- halo_samp  |>
  dplyr::count(Sample, Position, BrNum) |>
  mutate(data_type = "RNA Scope")

halo_combo_samp <- halo_samp  |>
  dplyr::group_by(Sample, Position, BrNum) |>
  summarize(combo = paste0(Combo, collapse = "_")) |>
  mutate(n = ifelse(combo == "Star_Circle", "Both", combo))

halo_demo <- halo_samp |> group_by(BrNum) |> summarize(n_RNAScope = n())

#### ALL DLPFC SAMPLES ####
all_dlpfc <- do.call("rbind", list(sn_n_samp, bulk_n_samp, halo_n_samp))

## Do we have 4 matched assays for each sample?
all_dlpfc <- all_dlpfc |> 
  left_join(all_dlpfc |> dplyr::count(Sample)|> mutate(Matched = n ==4) |> select(-n)) 

all_dlpfc

experiment_tile <- all_dlpfc |> 
  mutate(n = as.factor(n)) |>
  ggplot(aes(Position, BrNum, fill = n))+
  geom_tile()+
  geom_text(aes(label = n, color = Matched))+
  scale_color_manual(values =c(`FALSE` = "red", `TRUE` = "black")) +
  facet_wrap(~data_type, nrow = 1) +
  theme_bw()

ggsave(experiment_tile, filename = here(plot_dir, "experiment_tile.png"), width = 10)

ggsave(experiment_tile + 
         theme(axis.text.x=element_text(angle=45, hjust=1)),
       filename = here(plot_dir, "experiment_tile_small.pdf"), height = 5)

#### simplify n plot ####

all_dlpfc |>
  dplyr::group_by(data_type) |>
  summarize(sum = sum(n))

bulk_count <- bulk_samples |>
  dplyr::count(fraction = rna_extract, library_type, fill = library_combo2) |>
  mutate(data_type = "Bulk RNA-seq")
# 
# fraction library_type              fill  n    data_type
# 1     Bulk RiboZeroGold RiboZeroGold_Bulk 19 Bulk RNA-seq
# 2     Bulk        polyA        polyA_Bulk 19 Bulk RNA-seq
# 3     Cyto RiboZeroGold RiboZeroGold_Cyto 19 Bulk RNA-seq
# 4     Cyto        polyA        polyA_Cyto 18 Bulk RNA-seq
# 5      Nuc RiboZeroGold  RiboZeroGold_Nuc 17 Bulk RNA-seq
# 6      Nuc        polyA         polyA_Nuc 18 Bulk RNA-seq

sn_count <- sn_samples |>
  mutate(data_type = "snRNA-seq") |>
  dplyr::count(data_type) |>
  mutate(fill = "snRNA-seq", fraction = "Nuc", library_type = "polyA")

halo_count <- halo_samp |>
  dplyr::count(Combo) |>
  select(fill = Combo, fraction = Combo, n) |>
  mutate(data_type = "RNAscope",library_type = "RNAscope")

simple_count <- bind_rows(bulk_count, sn_count, halo_count)

simple_count |> group_by(data_type) |> summarize(sum = sum(n))

simple_tile <- simple_count |>
  ggplot(aes(fraction, library_type, fill = fill)) +
  geom_tile() +
  geom_text(aes(label = n), color = "white") +
  scale_fill_manual(values = c(library_combo_colors2, `snRNA-seq` = "#2f7ec0", Star = "cyan", Circle = "pink")) +
  theme_bw() + 
  facet_wrap(~data_type, ncol = 1, scales = "free") +
  theme(legend.position='none') 

ggsave(simple_tile, filename = here(plot_dir, "simple_sample_tile.png"), height = 4, width = 4)

#### by tissue block ####
bulk_qc <- read_csv(here("processed-data", "02_quality_control", "preQC_colData.csv")) |>
  mutate(pass_qc = !(auto_drop | drop_pca),
         data_type = "Bulk RNA-seq")|>
  select(BrNum, Sample, data_type, pass_qc, subsample = library_combo2)

bulk_qc |> dplyr::count(pass_qc)

qc_tab <- sn_n_samp |>
  select(BrNum, Sample, data_type) |>
  mutate(pass_qc = TRUE, subsample = data_type)|> 
  rbind(halo_info |>
          mutate(pass_qc = Confidence %in% c("High", "OK"),
                 subsample = Combo,
                 data_type = "RNAScope/IF") |>
          select(BrNum, Sample, data_type, pass_qc, subsample)) |>
  rbind(bulk_qc) |>
  mutate(data_type = factor(data_type, levels= c("RNAScope/IF","snRNA-seq","Bulk RNA-seq")))

length(unique(qc_tab$Sample))
# [1] 22

table(jaffelab::ss(unique(qc_tab$Sample),'_',2))
# ant  mid post 
# 7    9    6 

qc_count_assay <- qc_tab |> 
  group_by(data_type) |> 
  summarise(preQC = n(), 
            postQC = sum(pass_qc)) |>
  mutate(data_type_anno = paste0(data_type,"\n(n=", preQC, ":", postQC, ")")) |>
  mutate(data_type_anno = fct_reorder(data_type_anno, as.integer(data_type)))
# data_type   preQC postQC data_type_anno          
# <fct>       <int>  <int> <fct>                   
# 1 RNAScope/IF    42     25 "RNAScope/IF\n(n=42:25)"
# 2 snRNA-seq      19     19 "snRNA-seq\n(n=19:19)"  
# 3 Bulk RNA      113    110 "Bulk RNA\n(n=113:110)" 

qc_count_sub <- qc_tab |> 
  group_by(data_type, subsample) |> 
  summarise(preQC = n(), 
            postQC = sum(pass_qc)) |>
  mutate(subsample_anno = paste0(subsample," (n=", preQC, ":", postQC, ")"))
#   data_type   subsample         preQC postQC subsample_anno             
#   <fct>       <chr>             <int>  <int> <chr>                      
# 1 RNAScope/IF Circle                21     12 Circle (n=21:12)            
# 2 RNAScope/IF Star                  21     13 Star (n=21:13)              
# 3 snRNA-seq   snRNA-seq             19     19 snRNA-seq (n=19:19)         
# 4 Bulk RNA    RiboZeroGold_Cyto     19     19 RiboZeroGold_Cyto (n=19:19) 
# 5 Bulk RNA    RiboZeroGold_Nuc      19     17 RiboZeroGold_Nuc (n=19:17)  
# 6 Bulk RNA    RiboZeroGold_Total    19     19 RiboZeroGold_Total (n=19:19)
# 7 Bulk RNA    polyA_Cyto            19     18 polyA_Cyto (n=19:18)        
# 8 Bulk RNA    polyA_Nuc             18     18 polyA_Nuc (n=18:18)         
# 9 Bulk RNA    polyA_Total           19     19 polyA_Total (n=19:19)  

qc_tile <- qc_tab |> 
  left_join(qc_count_assay |> select(-preQC, -postQC)) |>
  left_join(qc_count_sub |> select(-preQC, -postQC)) |>
  # ggplot(aes(x = Sample, y = subsample, fill = pass_qc)) +
  ggplot(aes(x = Sample, y = subsample_anno, fill = pass_qc)) +
  geom_tile(color = "black") +
  # facet_wrap(~data_type, ncol = 1, scales = "free_y") + 
  facet_grid(data_type_anno~BrNum, scales="free", switch="y") + 
  scale_y_discrete(position = "right") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y.left = ,
        legend.position = "bottom") +
  labs(x = "Tissue Blocks (n=22)", y = "")

ggsave(qc_tile, filename = here(plot_dir, "qc_tile.png"), height = 5, width = 7.25)
ggsave(qc_tile, filename = here(plot_dir, "qc_tile.pdf"), height = 5, width = 7.25)

#### build demographics table ####
br_qc <- qc_tab |> 
  group_by(data_type, BrNum) |> 
  summarise(n_qc = paste0(sum(pass_qc),'(', n(), ')')) |>
  pivot_wider(names_from = "data_type", values_from = n_qc)

all_demo <- sn_demo |>
  select(-n_snRNA) |>
  mutate(PrimaryDx = "Control") |>
  left_join(br_qc)

write_csv(all_demo, file = here("processed-data", "00_data_prep", "donor_demographics.csv"))
