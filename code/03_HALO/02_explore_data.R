library("tidyverse")
library("scales")
library("here")
library("sessioninfo")
library("DeconvoBuddies")

#### Plot Set-up ####
plot_dir <- here("plots", "03_HALO", "01_explore_data")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

load(here("processed-data", "00_data_prep", "cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo

#### Metadata ####
kelsey_notes <- read.csv(here("processed-data", "03_HALO", "Deconvolution_HALO_Analysis_kelsey_googlesheet.csv")) |>
    mutate(Glare = grepl("GLARE", Comments.Issues))

## Br2723 -> Br2743 Fix may result in some mismatches with server files
kelsey_notes <- kelsey_notes |>
    mutate(Section = gsub("2723", "2743", Section)) |>
    mutate(typo = !Vectorize(grepl)(Section, `Path.to..Final..csv.File`))

kelsey_notes |> count(typo)

## Refine details about slides
slide_tab <- kelsey_notes |>
    select(Section, Slide, Combo = Combination, Glare) |>
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
    left_join(slide_tab)

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
# BrNum      n
# <chr>  <int>
# 1 Br2720     4
# 2 Br2723     2 # error - probably Br2743
# 3 Br3942     6
# 4 Br6423     4
# 5 Br6432     6
# 6 Br6471     4
# 7 Br6522     4
# 8 Br8325     4
# 9 Br8492     4
# 10 Br8667    4

metadata |> count(Position)
# Position      n
# <fct>     <int>
# 1 Anterior     14
# 2 Middle       16
# 3 Posterior    12

metadata |> count(Combo, Slide, subslide2)

## Write metadata
write_csv(metadata, file = here("processed-data", "03_HALO", "HALO_metadata.csv"))

#### Halo files ####
halo_runs <- c("prelim", "Deconvolution_HALO_analysis", "Algorithm_Check_20220914", "Algorithm_Check_20220920", "Algorithm_Check_20221020")
names(halo_runs) <- c("prelim", "Aug10", "Sep14", "Sep20", "Oct20")

halo_files <- map(halo_runs, function(dir) {
    hf <- list.files(path = here("raw-data", "HALO", dir), full.names = TRUE, recursive = TRUE, pattern = ".csv")

    sample <- gsub("^.*(\\d\\d\\d\\d[APM])_.*$", "\\1", basename(hf))
    combo <- gsub("^.*(Star|Circle|STAR|CIRCLE).*$", "\\1", basename(hf))
    names(hf) <- paste0("Br", sample, "_", toupper(combo))
    return(hf)
})

map(halo_files, names)

## maybe delete in files?
halo_files$Sep20$`BrHALO_Allsections_round4_AKT3_Circle_Fused.tif_object_Data.csv_CIRCLE` <- NULL

map(halo_files, length)

## any mismatches?
map(halo_files, ~ all(names(.x) %in% metadata$SAMPLE_ID))
# map(halo_files, ~.x[!names(.x) %in% metadata$SAMPLE_ID])

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

## Eventually run all of the HALO runs
halo_tables <- map(halo_files$Oct20, read_halo)

## TODO move when git is fixed
# #### Exploratory n nuc plots ####
# summary(metadata$n_nuc)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# # 22677   40586   46286   46097   53010   63293
# sd(metadata$n_nuc)
# # [1] 9953.509
#
# n_nuc_col <- ggplot(metadata, aes(x = Sample, y = n_nuc, fill = Combo)) +
#   geom_col(position = "dodge") +
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 45, hjust=1))
#
# ggsave(n_nuc_col, filename = here(plot_dir, "n_nuc_col.png"), width = 10)
#
#
# n_nuc_box <- ggplot(metadata, aes(x = Slide, y = n_nuc, fill = Combo)) +
#   geom_boxplot(position = "dodge") +
#   theme_bw()
#
# ggsave(n_nuc_box, filename = here(plot_dir, "n_nuc_box.png"))
#
# ## checkout distribution by combo
# n_nuc_box_combo <- ggplot(metadata, aes(x = Combo, y = n_nuc, fill = Combo)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.2)+
#   theme_bw() +
#   theme(legend.position = "None")
#
# ggsave(n_nuc_box_combo, filename = here(plot_dir, "n_nuc_box_combo.png"))
#
# n_nuc_box_round <- ggplot(metadata, aes(x = Round, y = n_nuc, fill = Combo)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   geom_point(aes(color = Glare), position=position_jitterdodge(jitter.width = 0.2))+
#   theme_bw()
#
# ggsave(n_nuc_box_round, filename = here(plot_dir, "n_nuc_box_round.png"))
#
# n_nuc_box_pos <- ggplot(metadata, aes(x = Position, y = n_nuc, fill = Position)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.2)+
#   facet_wrap(~Combo) +
#   theme_bw() +
#   theme(legend.position = "None")
#
# ggsave(n_nuc_box_pos, filename = here(plot_dir, "n_nuc_box_position.png"))
#
# ## Explore Slide + Splide postion of Samples
# slide_tile_nuc <-  metadata |>
#   ggplot(aes(x = Slide, y = subslide2, fill = n_nuc)) +
#   geom_tile() +
#   geom_text(aes(label = Sample), size = 2.5)+
#   facet_wrap(~Combo) +
#   scale_fill_viridis()+
#   theme_bw()
#
# ggsave(slide_tile_nuc, filename = here(plot_dir, "n_nuc_slide_tile.png"), width = 10)
#
#
# slide_tile_round <-  metadata |>
#   ggplot(aes(x = Slide, y = subslide2, fill = Round)) +
#   geom_tile() +
#   geom_text(aes(label = Sample), size = 2.5)+
#   facet_wrap(~Combo) +
#   theme_bw()
#
# ggsave(slide_tile_round, filename = here(plot_dir, "round_slide_tile.png"), width = 10)

#### Cell Type Annotation ####

ct_markers <- tibble(
    cellType = c("Endo", "Astro", "Inhib", "Excit", "Micro", "Oligo"),
    marker = c("CLDN5", "GFAP", "GAD1", "SLC17A7", "TMEM119", "OLIG2"),
    Combo = rep(c("CIRCLE", "STAR"), each = 3)
)

write_csv(ct_markers, here("processed-data", "03_HALO", "CellType_markers.csv"))

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
# [1] 1011916      46

halo_circle |> count(GAD1, GFAP, CLDN5)
#    GAD1 GFAP CLDN5      n
# 1:    0    0     0 616882
# 2:    0    0     1  59853
# 3:    0    1     0 240751
# 4:    1    0     0  94430

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
# cell_type      n
# 1:     Astro 240751
# 2:      Endo  59853
# 3:     Inhib  94430
# 4:     Other 616882

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
cn <- colnames(halo_tables$Br6432A_STAR)

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
# 1:       0       0     0 473785
# 2:       0       0     1 111398
# 3:       0       1     0  39152
# 4:       1       0     0 157718

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

save(halo_star, halo_circle, file = here("processed-data", "03_HALO", "HALO_Data.Rdata"))

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
exclude_data <- c(
    "STAR_2720M_HALO_RESCAN_6432A_Star_AKT3_Fused.tif_object_Data.csv",
    "STAR_6432A_HALO_allsections_AKT3_Star_Fused.tif_object_Data.csv"
)

halo_all <- rbind(halo_star, halo_circle) |>
    left_join(halo_file_table) |>
    left_join(metadata) |>
    relocate(SAMPLE_ID, Sample, BrNum, Position, Object_Id, cell_type, AKT3_Copies, Nucleus_Area, XMin, XMax, YMin, YMax) |>
    filter(!basename %in% exclude_data)

halo_all |>
    group_by(SAMPLE_ID, basename) |>
    summarize(n = n())

halo_all |>
    distinct(SAMPLE_ID, basename) |>
    count(SAMPLE_ID) |>
    filter(n > 1)

save(halo_all, file = here("processed-data", "03_HALO", "halo_all.Rdata"))

## TODO move when git is fixed, end previous script at exporting HALO_all
#### Composition Plots ####
circle_colors <- cell_type_colors_halo[c("Astro", "Endo", "Inhib", "Other")]
star_colors <- cell_type_colors_halo[c("Excit", "Micro", "Oligo", "Other")]


cell_type_prop <- halo_all |>
    filter(!basename %in% exclude_data) |>
    group_by(basename, SAMPLE_ID, Sample, Combo, rescan, cell_type) |>
    count() |>
    group_by(basename, SAMPLE_ID, Sample, Combo, rescan) |>
    mutate(prop = n / sum(n))

cell_type_prop |> filter(Sample == "Br8667_mid")

circle_prop_bar <- plot_composition_bar(circle_prop, sample_col = "Sample", x_col = "Sample") +
    scale_fill_manual(values = circle_colors) +
    labs(title = "Circle Combo") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(circle_prop_bar, filename = here(plot_dir, "prop_bar_circle.png"), width = 12)

circle_count_bar <- halo_circle |>
    count(SAMPLE_ID, cell_type) |>
    ggplot(aes(x = SAMPLE_ID, y = n, fill = cell_type)) +
    geom_col() +
    scale_fill_manual(values = circle_colors) +
    labs(title = "Circle Combo") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(circle_count_bar, filename = here(plot_dir, "count_bar_circle.png"), width = 12)


star_prop_bar <- plot_composition_bar(star_prop, sample_col = "Sample", x_col = "Sample") +
    scale_fill_manual(values = star_colors) +
    labs(title = "Star Combo") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(star_prop_bar, filename = here(plot_dir, "prop_bar_star.png"), width = 12)

star_count_bar <- halo_star |>
    count(SAMPLE_ID, cell_type) |>
    ggplot(aes(x = SAMPLE_ID, y = n, fill = cell_type)) +
    geom_col() +
    scale_fill_manual(values = star_colors) +
    labs(title = "star Combo") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(star_count_bar, filename = here(plot_dir, "count_bar_star.png"), width = 12)

##
# all_count_bar <- halo_all |>
#   count(Sample, Combo, cell_type) |>
#   ggplot(aes(x = Sample, y = n,fill = cell_type))+
#   geom_col()+
#   scale_fill_manual(values = cell_type_colors_halo) +
#   facet_wrap(~Combo, ncol = 1) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust=1))
#
# ggsave(all_count_bar, filename = here(plot_dir, "count_bar_ALL.png"), width = 12)

all_count_bar <- cell_type_prop |>
    ggplot(aes(x = Sample, y = n, fill = cell_type)) +
    geom_col() +
    scale_fill_manual(values = cell_type_colors_halo) +
    facet_wrap(~Combo, ncol = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Number of Cells")

ggsave(all_count_bar, filename = here(plot_dir, "count_bar_ALL_Oct20.png"), width = 12)
ggsave(all_count_bar, filename = here(plot_dir, "count_bar_ALL_Oct20_small.png"), width = 6, height = 5)


## Combine ##
# prop_all <- rbind(circle_prop, star_prop)

prop_boxplots <- prop_all |>
    ggplot(aes(x = Round, y = prop, color = cell_type)) +
    geom_boxplot() +
    scale_color_manual(values = cell_type_colors_halo) +
    facet_wrap(~cell_type, scales = "free_y")

ggsave(prop_boxplots, filename = here(plot_dir, "prop_boxplots.png"))

prop_boxplots <- prop_all |>
    filter(cell_type != "Other") |>
    ggplot(aes(x = cell_type, y = prop, fill = cell_type)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.2, colour = "black", pch = 21) +
    scale_fill_manual(values = cell_type_colors_halo) +
    # scale_color_manual(values = cell_type_colors_broad) +
    theme_bw() +
    theme(legend.position = "None") +
    labs(title = "RNAscope Cell Type Proportions")

ggsave(prop_boxplots, filename = here(plot_dir, "prop_boxplots.png"))


prop_boxplot_position <- prop_all |>
    ggplot(aes(x = Position, y = prop, fill = cell_type)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.2, colour = "black", pch = 21) +
    scale_fill_manual(values = cell_type_colors_halo) +
    facet_wrap(~cell_type, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "None") +
    labs(title = "RNAscope Cell Type Proportions")

ggsave(prop_boxplot_position, filename = here(plot_dir, "prop_boxplot_position.png"))


prop_other_adj <- prop_all |>
    filter(cell_type != "Other") |>
    group_by(Sample) |>
    summarize(
        cell_type = "Other_est",
        prop = 1 - sum(prop)
    )

prop_all_adj <- prop_all |>
    filter(cell_type != "Other") |>
    select(Sample, cell_type, prop) |>
    rbind(prop_other_adj)

prop_all_adj |>
    group_by(Sample) |>
    summarize(sum(prop))

adj_prop_bar <- plot_composition_bar(prop_all_adj, sample_col = "Sample", x_col = "Sample", min_prop_text = .01) +
    scale_fill_manual(values = cell_type_colors_halo) +
    labs(title = "Estimated Sample Compositions") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(adj_prop_bar, filename = here(plot_dir, "prop_bar_adj.png"), width = 12)


prop_all |>
    filter(prop != 1) |>
    group_by(Combo, cell_type) |>
    summarize(
        min = min(prop),
        mean = mean(prop),
        median = median(prop),
        max = max(prop)
    )
# Combo  cell_type      min   mean median    max
# <chr>  <chr>        <dbl>  <dbl>  <dbl>  <dbl>
# 1 Circle Astro     0.0487   0.177  0.136  0.592
# 2 Circle Endo      0.000291 0.0396 0.0374 0.0877
# 3 Circle Inhib     0.0479   0.103  0.110  0.142
# 4 Circle Other     0.325    0.680  0.705  0.840
# 5 Star   Excit     0.0174   0.246  0.274  0.342
# 6 Star   Micro     0.000915 0.0459 0.0426 0.0948
# 7 Star   Oligo     0.000131 0.0868 0.0308 0.319
# 8 Star   Other     0.490    0.625  0.629  0.741

write_csv(prop_all, file = here("processed-data", "03_HALO", "HALO_cell_type_proportions.csv"))

## cell type correlations
prop_wide <- prop_all |>
    filter(cell_type != "Other") |>
    ungroup() |>
    select(Sample, cell_type, prop) |>
    pivot_wider(names_from = "cell_type", values_from = "prop")

library("GGally")

prop_pairs <- ggpairs(prop_wide, columns = 2:7)
ggsave(prop_pairs, filename = here(plot_dir, "prop_pairs.png"))

n_wide <- prop_all |>
    filter(cell_type != "Other") |>
    ungroup() |>
    select(Sample, cell_type, n) |>
    pivot_wider(names_from = "cell_type", values_from = "n")

n_pairs <- ggpairs(n_wide, columns = 2:7)
ggsave(n_pairs, filename = here(plot_dir, "n_nuc_pairs.png"))
