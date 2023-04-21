library("SummarizedExperiment")
library("here")
library("tidyverse")
library("SingleCellExperiment")
library("sessioninfo")

#### Plot Setup ####
plot_dir = here("plots","00_data_prep","data_standards")
if(!dir.exists(plot_dir)) dir.create(plot_dir)


pos_df <- tibble(Position = c("Anterior", "Middle", "Posterior"), 
                pos = c("ant", "mid", "post"))

#### Review terms in snRNA-seq ####
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

sce_pd <- colData(sce) |> as.data.frame()

## save colData for eazy access
save(sce_pd, file = here("processed-data", "sce", "sce_pd.Rdata"))

sn_samples <- sce_pd |>
  count(Sample, SAMPLE_ID, Position, pos, round, BrNum, age, sex) 

#         Sample SAMPLE_ID  Position  pos  round  BrNum   age sex    n
# 1   Br2720_mid      1c-k    Middle  mid round1 Br2720 48.22   F 3101
# 2  Br2720_post      9c_k Posterior post round3 Br2720 48.22   F 5911
# 3   Br2743_ant      8c_k  Anterior  ant round3 Br2743 61.54   M 2861
# 4   Br2743_mid    round0    Middle  mid round0 Br2743 61.54   M 2723
# 5   Br3942_ant     11c_k  Anterior  ant round3 Br3942 47.53   M 5205
# 6   Br3942_mid     15c_k    Middle  mid round4 Br3942 47.53   M 4282
# 7   Br6423_ant     14c_k  Anterior  ant round4 Br6423 51.73   M 3898
# 8  Br6423_post     12c_k Posterior post round4 Br6423 51.73   M 4067
# 9   Br6432_ant      2c-k  Anterior  ant round1 Br6432 48.88   M 3059
# 10  Br6471_ant      3c-k  Anterior  ant round1 Br6471 55.46   M 3212
# 11  Br6471_mid      6c_k    Middle  mid round2 Br6471 55.46   M 4724
# 12  Br6522_mid      4c_k    Middle  mid round2 Br6522 33.39   M 4004
# 13 Br6522_post      5c_k Posterior post round2 Br6522 33.39   M 4275
# 14  Br8325_ant     16c_k  Anterior  ant round5 Br8325 57.62   F 4707
# 15  Br8325_mid     17c_k    Middle  mid round5 Br8325 57.62   F 4020
# 16  Br8492_mid     10c_k    Middle  mid round3 Br8492 53.40   F 4997
# 17 Br8492_post      7c_k Posterior post round2 Br8492 53.40   F 2661
# 18  Br8667_ant     18c_k  Anterior  ant round5 Br8667 37.33   F 5774
# 19  Br8667_mid     13c_k    Middle  mid round4 Br8667 37.33   F 4123

## need to fix:
# round -> Round (?)

sn_n_samp <- sn_samples |>
  count(Sample, Position, BrNum) |>
  mutate(data_type = "snRNA-seq")

#### Bulk info ####
bulk_samples <- read.csv(here("processed-data","01_SPEAQeasy","data_info.csv"))
head(bulk_samples)[,1:7]
#                       SAMPLE_ID       Dataset  BrNum location library_prep library_type round
# 1 2107UNHS-0291_Br2720_Mid_Cyto 2107UNHS-0291 Br2720      Mid         Cyto        polyA     1
# 2  2107UNHS-0291_Br2720_Mid_Nuc 2107UNHS-0291 Br2720      Mid          Nuc        polyA     1
# 3      2107UNHS-0291_Br2720_Mid 2107UNHS-0291 Br2720      Mid         Bulk        polyA     1
# 4 2107UNHS-0291_Br6432_Ant_Cyto 2107UNHS-0291 Br6432      Ant         Cyto        polyA     1
# 5  2107UNHS-0291_Br6432_Ant_Nuc 2107UNHS-0291 Br6432      Ant          Nuc        polyA     1
# 6      2107UNHS-0291_Br6432_Ant 2107UNHS-0291 Br6432      Ant         Bulk        polyA     1

bulk_samples |> count(location)

## Add Sample(?) BrXXXX_mid to match other data - how to handle Library prep?
# round -> Round (?)
# location -> Position & full name

## temp fix 
bulk_samples <- bulk_samples |>
  mutate(pos = tolower(location), 
         Sample = paste0(BrNum, "_", pos)) |>
  left_join(pos_df)

bulk_samples |> count(library_prep, library_type)

bulk_n_samp <- bulk_samples |>
  count(Sample, Position, BrNum, library_type) |>
  mutate(data_type = paste0("Bulk RNA-seq\n", library_type)) |>
  select(-library_type)


#### HALO Info ####
halo_info <- read.csv(here("processed-data","03_HALO","HALO_metadata.csv"))
head(halo_info)
#         SAMPLE_ID Round Section  BrNum Pos  Combo  Position     Sample Slide subslide i subslide2 n_nuc
# 1 R1_2720M_Circle    R1   2720M Br2720   M Circle    Middle Br2720_mid     5        A 1        A1 28239
# 2 R1_6432A_Circle    R1   6432A Br6432   A Circle  Anterior Br6432_ant     5        D 1        D1 61378
# 3 R1_6432M_Circle    R1   6432M Br6432   M Circle    Middle Br6432_mid     2        A 1        A1 48252
# 4 R1_6432P_Circle    R1   6432P Br6432   P Circle Posterior Br6432_pos     5        B 1        B1 47951
# 5 R1_6471A_Circle    R1   6471A Br6471   A Circle  Anterior Br6471_ant     5        C 1        C1 43790
# 6 R2_6471M_Circle    R2   6471M Br6471   M Circle    Middle Br6471_mid     4        C 1        C1 59392


halo_samp <- halo_info |> 
  select(SAMPLE_ID, BrNum, Position, Combo) |>
  left_join(pos_df) |>
  mutate(Sample = paste0(BrNum, "_", pos))

halo_n_samp <- halo_samp |>
  count(Sample, Position, BrNum) |>
  mutate(data_type = "RNA Scope")

#### Use these colnames ####
## SAMPLE_ID - unique sample identifier
## Sample - BrXXXX_mid - identifies tissue ... hmm maybe call tissue?
## BrNum - BrXXXX
## Position - Anterior, Middle, Posterior
## pos - short lowercase Position?

#### ALL DLPFC SAMPLES ####

all_dlpfc <- do.call("rbind", list(sn_n_samp, bulk_n_samp, halo_n_samp))

## Do we have 4 matched assays for each sample?
all_dlpfc <- all_dlpfc |> 
  left_join(all_dlpfc |> count(Sample)|> mutate(Matched = n ==4) |> select(-n)) 

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
       filename = here(plot_dir, "experiment_tile_small.png"), height = 5)

# sgejobs::job_single('data_standards', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript data_standards.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
