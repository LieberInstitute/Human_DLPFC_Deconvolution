library("SummarizedExperiment")
library("here")
library("tidyverse")

#### Plot Setup ####
plot_dir = here("plots","00_data_prep","data_standards")
if(!dir.exists(plot_dir)) dir.create(plot_dir)


pos_df <- tibble(Position = c("Anterior", "Middle", "Posterior"), 
                pos = c("ant", "mid", "post"))

#### Review terms in snRNA-seq ####
sn_samples <- read.csv("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/03_build_sce/sample_info.csv")
head(sn_samples)
#        Sample file_id    region subject  round    n n_high_mito n_low_sum n_low_detect n_discard_auto
# 1  Br2720_mid    1c-k    middle  Br2720 round1 3743         544       110          306            642
# 2 Br2720_post    9c_k posterior  Br2720 round3 5950          26         0           13             39
# 3  Br2743_ant    8c_k  anterior  Br2743 round3 3198         337         0            0            337
# 4  Br2743_mid  round0    middle  Br2743 round0 3426         586        76          412            703
# 5  Br3942_ant   11c_k  anterior  Br3942 round3 6269        1064         0            0           1064
# 6  Br3942_mid   15c_k    middle  Br3942 round4 4309          27         0            0             27

## need to fix:
# round -> Round (?)
# subject -> BrNum
# region -> Position & capitalize first letter

## temp fix 
sn_samples <- sn_samples |>
  rename(Position = region, BrNum = subject) |>
  mutate(Position = str_to_title(Position))

sn_n_samp <- sn_samples |>
  count(Sample, Position, BrNum) |>
  mutate(data_type = "snRNA-seq")

#### Bulk info ####
bulk_samples <- read.csv(here("processed-data","01_SPEAQeasy","data_info.csv"))
head(bulk_samples)
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
  mutate(data_type = paste("Bulk RNA-seq", library_type)) |>
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
