library("googlesheets4")
info <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1bld6g-7MN2G18b8hwXEkC0dDOhWlzmFTAmpcPYdgXjI/edit?usp=sharing",
  "Sheet1"
)
info <- info[seq_len(42), ] ## Drop the extra rows with info
stopifnot(max(table(with(info, paste(Section, "_", Combination)))) == 1)

jhpce <- c(
  "raw-data/RNAscope_images/raw_images_to_share/HALO_Round3_BO_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round3_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_RESCAN_2720M_B1R1_Circle_AKT3_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section3_Br6471_Ant_Circle.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round2_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section3_Br6471_Ant_Star.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_Allsections_round4_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_Round3_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_B1_R5_Star_Br8325M_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_B1_R5_Circle_Br8325A_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_B1_R5_Circle_Br8667A_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section2_Br6432_Post_Star.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_RESCAN_6432M_B1R1_BO_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section1_Br6432_Ant_Circle.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_B1_R5_Star_Br8667A_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round2_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section1_Br2720_Mid_Star.tif",
  "raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section1_Br6432_Mid_Circle.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_B1_R5_Star_Br8325A_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section1_Br6432_Post_Circle.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_B1_R5_Circle_Br8325M_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_round3_BO_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/raw_images_to_share/HALO_RESCAN_6432A_Star_AKT3_Fused.tif"
)

kelsey <- gsub(".*(\\\\|/)", "", paste0(info$`Path to Fused/Stitched Images`, ".tif"))
m <- match(kelsey, basename(jhpce))
stopifnot(all(!is.na(m)))

info$jhpce_path <- jhpce[m]

## To print the info so I can add it to Google Sheets
cat(paste0(info$jhpce_path, "\n"))

raw-data/RNAscope_images/raw_images_to_share/HALO_RESCAN_6432A_Star_AKT3_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section1_Br2720_Mid_Star.tif
raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section2_Br6432_Post_Star.tif
raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section3_Br6471_Ant_Star.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_RESCAN_6432M_B1R1_BO_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_B1_R5_Star_Br8325A_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_B1_R5_Star_Br8325M_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_B1_R5_Star_Br8667A_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round2_AKT3_Star_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round2_AKT3_Star_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round2_AKT3_Star_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round2_AKT3_Star_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_Round3_AKT3_Star_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_Round3_AKT3_Star_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_Round3_AKT3_Star_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_Round3_AKT3_Star_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_round3_BO_AKT3_Star_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_AKT3_Star_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_AKT3_Star_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_AKT3_Star_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_AKT3_Star_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section1_Br6432_Ant_Circle.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_RESCAN_2720M_B1R1_Circle_AKT3_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section1_Br6432_Post_Circle.tif
raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section3_Br6471_Ant_Circle.tif
raw-data/RNAscope_images/raw_images_to_share/Fused_HALO_section1_Br6432_Mid_Circle.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_B1_R5_Circle_Br8325A_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_B1_R5_Circle_Br8325M_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_B1_R5_Circle_Br8667A_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round2_AKT3_Circle_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round2_AKT3_Circle_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round2_AKT3_Circle_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round2_AKT3_Circle_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round3_AKT3_Circle_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round3_AKT3_Circle_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round3_AKT3_Circle_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_allsections_round3_AKT3_Circle_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_Round3_BO_AKT3_Circle_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_Allsections_round4_AKT3_Circle_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_Allsections_round4_AKT3_Circle_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_Allsections_round4_AKT3_Circle_Fused.tif
raw-data/RNAscope_images/raw_images_to_share/HALO_Allsections_round4_AKT3_Circle_Fused.tif
