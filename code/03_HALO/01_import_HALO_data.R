
library("tidyverse")
library("scales")
library("here")
library("sessioninfo")
library("DeconvoBuddies")
library("readxl")

#### Set-up ####
# plot_dir <- here("plots", "03_HALO", "01_import_HALO_data")
# if (!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "03_HALO", "01_import_HALO_data")
if (!dir.exists(data_dir)) dir.create(data_dir)

load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo

#### Metadata ####
kelsey_notes <- read_excel(here("raw-data", "HALO", "Annotation_refinement_20230922", "Deconvolution_HALO_Analysis_Refined_Annotations.xlsx"),
                           sheet = 2) |>
  mutate(Glare = grepl("GLARE", `Comments_Issues`))

colnames(kelsey_notes)

# kelsey_notes |> filter(grepl("2723", Section))

## Br2723 -> Br2743 Fix may result in some mismatches with server files
kelsey_notes <- kelsey_notes |>
    mutate(Section = gsub("2723", "2743", Section)) |>
    mutate(typo = !Vectorize(grepl)(Section, `Path to Final csv File`))

kelsey_notes |> count(typo)
kelsey_notes |> count(`Confidence in overall quality`)

## Refine details about slides
slide_tab <- kelsey_notes |>
    select(Section, Slide, Combo = Combination, Glare, Excluded, Confidence = `Confidence in overall quality`) |>
    separate(Slide, into = c("Slide", "subslide"), sep = ", ") |>
    group_by(Slide, subslide, Combo) |>
    mutate(
        i = row_number(),
        subslide2 = factor(paste0(subslide, i)),
        Slide = factor(Slide)
    ) |>
    ungroup()

slide_tab |> count(Slide, Combo)
slide_tab |> count(Combo, Slide, subslide2)

slide_tab |> count(Excluded, Confidence)
# Excluded Confidence     n
# <chr>    <chr>      <int>
# 1 Maybe    Low            1
# 2 No       High          13
# 3 No       Low            7
# 4 No       OK            13
# 5 Yes      Excluded       8

## Position Names
position_codes <- data.frame(
    Pos = c("A", "M", "P"),
    pos = c("ant", "mid", "post"),
    Position = factor(c("Anterior", "Middle", "Posterior"))
)


metadata <- kelsey_notes |>
    separate(Section, into = c("BrNum", "Pos"), sep = "(?<=[0-9])(?=[A-Z])", remove = FALSE) |>
    left_join(position_codes) |>
    mutate(
        BrNum = paste0("Br", BrNum),
        Sample = paste0(BrNum, "_", pos),
        combo = toupper(Combination),
        SAMPLE_ID = paste0(BrNum, Pos, "_", combo)
    ) |>
    as_tibble() |>
    select(SAMPLE_ID, Sample, BrNum, Position, Pos, pos, Section, Round, Combo = Combination) |>
    left_join(slide_tab) |>
  mutate(Confidence = factor(Confidence, levels = c("High", "OK", "Low", "Excluded")))

# SAMPLE_ID    Sample      BrNum  Position  Pos   pos   Section Round Combo Slide subslide Glare     i subslide2
# 1 Br6432A_STAR Br6432_ant  Br6432 Anterior  A     ant   6432A       1 Star  6     D        FALSE     1 D1
# 2 Br2720M_STAR Br2720_mid  Br2720 Middle    M     mid   2720M       1 Star  6     A        FALSE     1 A1
# 3 Br6432P_STAR Br6432_post Br6432 Posterior P     post  6432P       1 Star  6     B        FALSE     1 B1
# 4 Br6471A_STAR Br6471_ant  Br6471 Anterior  A     ant   6471A       1 Star  6     C        FALSE     1 C1
# 5 Br6432M_STAR Br6432_mid  Br6432 Middle    M     mid   6432M       1 Star  2     B        TRUE      1 B1
# 6 Br8325A_STAR Br8325_ant  Br8325 Anterior  A     ant   8325A       5 Star  7     A        FALSE     1 A1
# 7 Br8325M_STAR Br8325_mid  Br8325 Middle    M     mid   8325M       5 Star  7     B        FALSE     1 B1
# 8 Br8667A_STAR Br8667_ant  Br8667 Anterior  A     ant   8667A       5 Star  7     C        FALSE     1 C1
# 9 Br6522M_STAR Br6522_mid  Br6522 Middle    M     mid   6522M       2 Star  5     A        FALSE     1 A1
# 10 Br6522P_STAR Br6522_post Br6522 Posterior P     post  6522P       2 Star  5     B        FALSE     1 B1

metadata |> count(Round)
metadata |> count(Combo)
# Combo     n
# <chr>      <int>
# 1 Circle        21
# 2 Star          21
metadata |> count(BrNum)
# # A tibble: 10 × 2
# BrNum      n
# <chr>  <int>
# 1 Br2720     4
# 2 Br2743     2
# 3 Br3942     6
# 4 Br6423     4
# 5 Br6432     6
# 6 Br6471     4
# 7 Br6522     4
# 8 Br8325     4
# 9 Br8492     4
# 10 Br8667     4

metadata |> count(Position)
# Position      n
# <fct>     <int>
# 1 Anterior     14
# 2 Middle       16
# 3 Posterior    12

metadata |> count(Combo, Slide, subslide2)

## Write metadata
write_csv(metadata, file = here(data_dir, "HALO_metadata.csv"))

#### Halo files ####
## file system from Sophia 2022 work
halo_runs <- c("prelim", "Deconvolution_HALO_analysis", "Algorithm_Check_20220914", "Algorithm_Check_20220920", "Algorithm_Check_20221020")
names(halo_runs) <- c("prelim", "Aug10_2022", "Sep14_20222", "Sep20_2022", "Oct20_2022")

halo_files <- map(halo_runs, function(dir) {
    hf <- list.files(path = here("raw-data", "HALO", dir), full.names = TRUE, recursive = TRUE, pattern = ".csv")

    sample <- gsub("^.*(\\d\\d\\d\\d[APM])_.*$", "\\1", basename(hf))
    combo <- gsub("^.*(Star|Circle|STAR|CIRCLE).*$", "\\1", basename(hf))
    names(hf) <- paste0("Br", sample, "_", toupper(combo))
    return(hf)
})

map(halo_files, names)

## maybe delete in files?
halo_files$Sep20_2022$`BrHALO_Allsections_round4_AKT3_Circle_Fused.tif_object_Data.csv_CIRCLE` <- NULL

## Add 2023 Refined annotations 
Sep22_2023 <- list.files(here("raw-data", "HALO", "Annotation_refinement_20230922"), recursive = TRUE, pattern = "object_results.csv")
names(Sep22_2023) <- gsub("tar","TAR", gsub("ircle","IRCLE",map_chr(strsplit(Sep22_2023, split = "/"),2)))
Sep22_2023 <- map(Sep22_2023, ~here("raw-data", "HALO", "Annotation_refinement_20230922", .x))
all(file.exists(unlist(Sep22_2023)))

halo_files$Sep22_2023 <- Sep22_2023

## any mismatches?
map(halo_files, ~ all(names(.x) %in% metadata$SAMPLE_ID))
# map(halo_files, ~.x[!names(.x) %in% metadata$SAMPLE_ID])

map_int(halo_files, length)

## Mismatches resolved
# mismatch_files <- unlist(map(halo_files, ~.x[!names(.x) %in% metadata$SAMPLE_ID]))

#### Rename files ####
## renamed files
# rename_files <- gsub("2723", "2743", mismatch_files)
# rename_files <- gsub("8776M", "8667M", rename_files)
# rename_files <- gsub("8876M", "8667M", rename_files)
#
# write.csv(data.frame(oldname = mismatch_files,
#            newname = rename_files),
#           file = here("processed-data", "03_HALO", "rename_log.csv"))

# file.rename(mismatch_files, rename_files)

# mismatch_files[mismatch_files == rename_files]

halo_file_table <- tibble(
    name = names(unlist(halo_files)),
    JHPCE_path = unlist(halo_files)
) |>
    separate(name, into = c("HALO_run", "SAMPLE_ID"), sep = "\\.") |>
    mutate(
        basename = basename(JHPCE_path),
        rescan = grepl("RESCAN", JHPCE_path)
    )
# |>
#   pivot_wider(names_from = "HALO_run", values_from = JHPCE_path)

## Duplicate data?
halo_file_table |>
    count(SAMPLE_ID, HALO_run) |>
    filter(n > 1)

# SAMPLE_ID    HALO_run     n
# <chr>        <chr>    <int>
# 1 Br6432A_STAR Oct20        3
# 2 Br6432A_STAR Sep20        2

halo_file_table |>
    group_by(SAMPLE_ID, HALO_run) |>
    filter(n() > 1) |>
    select(-JHPCE_path)

# HALO_run SAMPLE_ID    basename
# <chr>    <chr>        <chr>
# 1 Sep20    Br6432A_STAR STAR_6432A_HALO_allsections_AKT3_Star_Fused.tif_object_Data.csv
# 2 Sep20    Br6432A_STAR STAR_6432A_HALO_RESCAN_6432A_Star_AKT3_Fused.tif_object_Data.csv ## Use rescan?

# 3 Oct20    Br6432A_STAR STAR_2720M_HALO_RESCAN_6432A_Star_AKT3_Fused.tif_object_Data.csv ## Two SAMPLE_IDs?
# 4 Oct20    Br6432A_STAR STAR_6432A_HALO_allsections_AKT3_Star_Fused.tif_object_Data.csv
# 5 Oct20    Br6432A_STAR STAR_6432A_HALO_RESCAN_6432A_Star_AKT3_Fused.tif_object_Data.csv ## Use rescan?

#### Read csv files ####

read_halo <- function(fn) {
    halo_data <- data.table::fread(fn)
    halo_data$JHPCE_path <- fn
    return(halo_data)
}

## Read Annotation refinement files from Sep 2023
halo_tables <- map(halo_files$Sep22_2023, read_halo)

#### Cell Type Annotation ####

ct_markers <- tibble(
    cellType = c("Endo", "Astro", "Inhib", "Excit", "Micro", "Oligo"),
    marker = c("CLDN5", "GFAP", "GAD1", "SLC17A7", "TMEM119", "OLIG2"),
    Combo = rep(c("CIRCLE", "STAR"), each = 3)
)

write_csv(ct_markers, here(data_dir, "CellType_markers.csv"))

# cellType marker  Combo
# <chr>    <chr>   <chr>
# 1 Endo     CLDN5   Circle
# 2 Astro    GFAP    Circle
# 3 Inhib    GAD1    Circle
# 4 Excit    SLC17A7 Star
# 5 Micro    TMEM119 Star
# 6 Oligo    OLIG2   Star

## CIRCLE COMBO
halo_circle <- do.call("rbind", halo_tables[grep("CIRCLE", names(halo_tables))])
dim(halo_circle)
# [1] 767393     46

halo_circle |> count(GAD1, GFAP, CLDN5)
# GAD1 GFAP CLDN5      n
# 1:    0    0     0 474934
# 2:    0    0     1  47488
# 3:    0    1     0 167441
# 4:    1    0     0  77530

## cell cell types
halo_circle <- halo_circle |>
    mutate(
        n_marker = GAD1 + GFAP + CLDN5,
        cell_type = case_when(
            n_marker > 1 ~ "Multi",
            GAD1 == 1 ~ "Inhib",
            GFAP == 1 ~ "Astro",
            CLDN5 == 1 ~ "Endo",
            TRUE ~ "Other"
        )
    )

head(halo_circle)

summary(halo_circle$`AKT3 (Opal 570) Copies`)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.000   2.000   4.000   6.017   8.000  64.000


## no GAD1 copies?
head(halo_circle[, 1:5])
head(halo_circle[, grepl("Copies", colnames(halo_circle))])
# GFAP..Alexa.594..Copies AKT3..Opal.570..Copies CLDN5..Alexa.488..Copies


halo_circle |> count(cell_type)
#    cell_type      n
# 1:     Astro 167441
# 2:      Endo  47488
# 3:     Inhib  77530
# 4:     Other 474934

# circle_n_nuc <- halo_circle |>
#   group_by(SAMPLE_ID) |>
#   summarize(n_nuc = n())
#
# circle_prop <- halo_circle |>
#   count(SAMPLE_ID, cell_type) |>
#   right_join(metadata, .) |>
#   mutate(prop = n/n_nuc)
#
# circle_prop_wide <- circle_prop |>
#   mutate(prop = round(prop, 2)) |>
#   select(SAMPLE_ID, cell_type, prop) |>
#   pivot_wider(names_from = "cell_type", values_from = "prop")

#### Star combo ####

## need to fix R5 colnames SLC17A -> SLC17A7
# cn <- colnames(halo_tables$R1_2720M_Star)
cn <- colnames(halo_tables$Br2720P_STAR)

map(halo_tables[grep("STAR", names(halo_tables))], colnames)

map(halo_tables[grep("STAR", names(halo_tables))], ~ setdiff(colnames(.x), cn))

## fix typos
colnames(halo_tables$Br8325A_STAR) <- gsub("SLC17A\\&", "SLC17A7", colnames(halo_tables$Br8325A_STAR))
colnames(halo_tables$Br8325M_STAR) <- gsub("SLC17A\\&", "SLC17A7", colnames(halo_tables$Br8325M_STAR))
colnames(halo_tables$Br8667A_STAR) <- gsub("SLC17A\\&", "SLC17A7", colnames(halo_tables$Br8667A_STAR))


all(map_lgl(halo_tables[grep("STAR", names(halo_tables))], ~ all(colnames(.x) == cn)))

## rbind all star tables
halo_star <- do.call("rbind", halo_tables[grep("STAR", names(halo_tables))])
dim(halo_star)
halo_star |> count(SLC17A7, TMEM119, OLIG2)
#   SLC17A7 TMEM119 OLIG2      n
# 1:       0       0     0 355661
# 2:       0       0     1  68625
# 3:       0       1     0  25930
# 4:       1       0     0 144790

halo_star <- halo_star |>
    mutate(
        n_marker = SLC17A7 + TMEM119 + OLIG2,
        cell_type = case_when(
            n_marker > 1 ~ "Multi",
            SLC17A7 == 1 ~ "Excit",
            TMEM119 == 1 ~ "Micro",
            OLIG2 == 1 ~ "Oligo",
            TRUE ~ "Other"
        )
    )

halo_star |> count(cell_type)
# cell_type      n
# 1     Excit 199546
# 2     Micro  34504
# 3     Oligo  66579
# 4     Other 623519


halo_star |> count(JHPCE_path)
#
# star_prop <- halo_star |>
#   count(SAMPLE_ID, cell_type) |>
#   right_join(metadata, .) |>
#   mutate(prop = n/n_nuc)
#
# star_prop_wide <- star_prop |>
#   mutate(prop = round(prop, 2)) |>
#   select(SAMPLE_ID, cell_type, prop) |>
#   pivot_wider(names_from = "cell_type", values_from = "prop")
#
# star_prop |> filter(cell_type == 'Other', prop == 1) |> select(Sample, SAMPLE_ID, Round)

## save raw out put with more columns from HALO
save(halo_star, halo_circle, file = here(data_dir, "HALO_Data.Rdata"))

#### Build HALO ALL ####

# "AKT3 (Opal 620) Copies"
colnames(halo_star)[grepl("AKT3.*Copies", colnames(halo_star))] <- "AKT3_Copies"
# "AKT3 (Opal 570) Copies"
colnames(halo_circle)[grepl("AKT3.*Copies", colnames(halo_circle))] <- "AKT3_Copies"

# find common columns
(common_cols <- intersect(colnames(halo_star), colnames(halo_circle)))

## make more usable
common_cols_rename <- gsub(" \\(.*\\)", "", common_cols)
common_cols_rename <- gsub("/| ", "_", common_cols_rename)

## subset & rename
halo_star <- halo_star[, ..common_cols]
colnames(halo_star) <- common_cols_rename

halo_circle <- halo_circle[, ..common_cols]
colnames(halo_circle) <- common_cols_rename

## Exclude data from these files
# exclude_data <- c(
#     "STAR_2720M_HALO_RESCAN_6432A_Star_AKT3_Fused.tif_object_Data.csv",
#     "STAR_6432A_HALO_allsections_AKT3_Star_Fused.tif_object_Data.csv"
# )

halo_all <- rbind(halo_star, halo_circle) |>
  left_join(halo_file_table) |>
  left_join(metadata) |>
  relocate(SAMPLE_ID, Sample, BrNum, Position, Object_Id, cell_type, AKT3_Copies, Nucleus_Area, XMin, XMax, YMin, YMax) |>
  # filter(!basename %in% exclude_data) |>
  mutate(large_nuc = Nucleus_Area > pi * 5^2) |>
  as_tibble()


halo_all |>
    group_by(SAMPLE_ID, basename, Confidence) |>
    summarize(n = n())

halo_all |>
    distinct(SAMPLE_ID, basename) |>
    count(SAMPLE_ID) |>
    filter(n > 1)

halo_all |> count(Confidence)
metadata |> count(Confidence)

save(halo_all, file = here("processed-data", "03_HALO", "halo_all.Rdata"))

# sgejobs::job_single('01_import_HALO_data', create_shell = TRUE, memory = '10G', command = "Rscript 01_import_HALO_data.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
