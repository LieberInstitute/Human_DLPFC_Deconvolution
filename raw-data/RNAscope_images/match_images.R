library("googlesheets4")
info <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1bld6g-7MN2G18b8hwXEkC0dDOhWlzmFTAmpcPYdgXjI/edit?usp=sharing",
  "Sheet1"
)
info <- info[-43, ] ## Drop the extra row with info
stopifnot(max(table(with(info, paste(Section, "_", Combination)))) == 1)

jhpce <- c(
  "raw-data/RNAscope_images/03102022_RIF_DECON_B1_R5_Circle/03102022_RIF_DECON_B1_R5_Circle/HALO/HALO_B1_R5_Circle_Br8667A_Fused.tif",
  "raw-data/RNAscope_images/32422_RIF_DECON_B2_R4_Circle/32422_RIF_DECON_B2_R4_Circle/HALO/HALO_section4_Br8667M_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/32422_RIF_DECON_B2_R4_Circle/32422_RIF_DECON_B2_R4_Circle/HALO/HALO_section2_Br6423P_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/32422_RIF_DECON_B2_R4_Circle/32422_RIF_DECON_B2_R4_Circle/HALO/HALO_section3_Br3942M_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/32422_RIF_DECON_B2_R4_Circle/32422_RIF_DECON_B2_R4_Circle/HALO/HALO_section1_Br6423A_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_BO_Star/031822_RIF_DECON_B2_R3_BO_Star/REDO_HALO/HALO_round3_BO_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_BO_Star/031822_RIF_DECON_B2_R3_BO_Star/HALO/HALO_section2_Br3942A_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/03102022_RIF_DECON_B1_R5_Circle_Br8325A/03102022_RIF_DECON_B1_R5_Circle_Br8325A/HALO/HALO_B1_R5_Circle_Br8325A_Fused.tif",
  "raw-data/RNAscope_images/03102022_RIF_DECON_B1_R5_Star/03102022_RIF_DECON_B1_R5_Star/HALO_section3_Br8667A/HALO_B1_R5_Star_Br8667A_Fused.tif",
  "raw-data/RNAscope_images/03102022_RIF_DECON_B1_R5_Star/03102022_RIF_DECON_B1_R5_Star/HALO_section1_Br8325A/HALO_B1_R5_Star_Br8325A_Fused.tif",
  "raw-data/RNAscope_images/03102022_RIF_DECON_B1_R5_Star/03102022_RIF_DECON_B1_R5_Star/HALO_section2_Br8325M/HALO_B1_R5_Star_Br8325M_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B1_R2_Rerun_Circle/031822_RIF_DECON_B1_R2_Rerun_Circle/HALO/HALO_section1_Br6522M_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B1_R2_Rerun_Circle/031822_RIF_DECON_B1_R2_Rerun_Circle/HALO/HALO_section2_Br6522P_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B1_R2_Rerun_Circle/031822_RIF_DECON_B1_R2_Rerun_Circle/HALO/HALO_section3_Br6471P_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B1_R2_Rerun_Circle/031822_RIF_DECON_B1_R2_Rerun_Circle/HALO/HALO_section4_Br8492P_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/03102022_RIF_DECON_B1_R5_Circle_BR8325M/03102022_RIF_DECON_B1_R5_Circle_BR8325M/HALO/HALO_B1_R5_Circle_Br8325M_Fused.tif",
  "raw-data/RNAscope_images/02242022_RIF_DECON_B1_R1_RERUN_STAR/HALO/Fused_HALO_section3_Br6471_Ant_Star.tif",
  "raw-data/RNAscope_images/02242022_RIF_DECON_B1_R1_RERUN_STAR/HALO/Fused_HALO_section2_Br6432_Post_Star.tif",
  "raw-data/RNAscope_images/02242022_RIF_DECON_B1_R1_RERUN_STAR/HALO/Fused_HALO_section1_Br2720_Mid_Star.tif",
  "raw-data/RNAscope_images/02242022_RIF_DECON_B1_R1_RERUN_STAR/HALO/Fused_HALO_section4_Br6432_Ant_Star.tif",
  "raw-data/RNAscope_images/02242022_RIF_DECON_B1_R1_RERUN_CIRCLE/HALO/Fused_HALO_section1_Br6432_Post_Circle.tif",
  "raw-data/RNAscope_images/02242022_RIF_DECON_B1_R1_RERUN_CIRCLE/HALO/Fused_HALO_section1_Br6432_Ant_Circle.tif",
  "raw-data/RNAscope_images/02242022_RIF_DECON_B1_R1_RERUN_CIRCLE/HALO/Fused_HALO_section1_Br2720_Mid_Circle.tif",
  "raw-data/RNAscope_images/02242022_RIF_DECON_B1_R1_RERUN_CIRCLE/HALO/Fused_HALO_section3_Br6471_Ant_Circle.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_Circle/031822_RIF_DECON_B2_R3_Circle/HALO/HALO_section3_Br8492M_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_Circle/031822_RIF_DECON_B2_R3_Circle/HALO/HALO_section2_Br2720P_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_Circle/031822_RIF_DECON_B2_R3_Circle/HALO/HALO_section4_Br3942P_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_Circle/031822_RIF_DECON_B2_R3_Circle/HALO/HALO_section1_Br2723A_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_BO_Circle/031822_RIF_DECON_B2_R3_BO_Circle/HALO/HALO_section1_Br3942A_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_BO_Circle/031822_RIF_DECON_B2_R3_BO_Circle/REDO_HALO/HALO_Round3_BO_AKT3_Circle_Fused.tif",
  "raw-data/RNAscope_images/32422_RIF_DECON_B2_R4_Star/32422_RIF_DECON_B2_R4_Star/HALO/HALO_section2_Br6423P_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/32422_RIF_DECON_B2_R4_Star/32422_RIF_DECON_B2_R4_Star/HALO/HALO_section4_Br8667M_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/32422_RIF_DECON_B2_R4_Star/32422_RIF_DECON_B2_R4_Star/HALO/HALO_section3_Br3942M_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/32422_RIF_DECON_B2_R4_Star/32422_RIF_DECON_B2_R4_Star/HALO/HALO_section1_Br6423A_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_Star/031822_RIF_DECON_B2_R3_Star/HALO/HALO_section4_Br3942P_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_Star/031822_RIF_DECON_B2_R3_Star/HALO/HALO_section3_Br8492M_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_Star/031822_RIF_DECON_B2_R3_Star/HALO/HALO_section2_Br2720P_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_Star/031822_RIF_DECON_B2_R3_Star/HALO/HALO_section1_Br2723A_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B1_R2_Rerun_Star/031822_RIF_DECON_B1_R2_Rerun_Star/HALO/HALO_section4_Br8492P_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B1_R2_Rerun_Star/031822_RIF_DECON_B1_R2_Rerun_Star/HALO/HALO_section1_Br6522M_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B1_R2_Rerun_Star/031822_RIF_DECON_B1_R2_Rerun_Star/HALO/HALO_section3_Br6471M_AKT3_Star_Fused.tif",
  "raw-data/RNAscope_images/031822_RIF_DECON_B1_R2_Rerun_Star/031822_RIF_DECON_B1_R2_Rerun_Star/HALO/HALO_section2_Br6522P_AKT3_Star_Fused.tif"
)

potential_matches <- apply(
  info,
  1,
  function(row) {
    section <- grepl(row[["Section"]], jhpce, ignore.case = TRUE)
    if (!any(section)) {
      if (grepl("A", row[["Section"]])) {
        section <- grepl(gsub("A", "", row[["Section"]]), jhpce) &
          grepl("Ant", jhpce, ignore.case = TRUE)
      } else if (grepl("M", row[["Section"]])) {
        section <- grepl(gsub("M", "", row[["Section"]]), jhpce) &
          grepl("Mid", jhpce, ignore.case = TRUE)
      } else if (grepl("P", row[["Section"]])) {
        section <- grepl(gsub("P", "", row[["Section"]]), jhpce) &
          grepl("Post", jhpce, ignore.case = TRUE)
      }
    }
    section & grepl(row[["Combination"]], jhpce, ignore.case = TRUE)
  },
  simplify = FALSE
)

## Check that there is at most 1 potential match. That is, that there is
## no issue due to "REDO" or things like that.
stopifnot(max(sapply(potential_matches, sum)) == 1)

## Find the actual match
info$jhpce_path <- sapply(potential_matches, function(x) {
  if (length(which(x)) > 0) {
    jhpce[which(x)]
  } else {
    NA
  }
})

## Check the entries where this code couldn't find a JHPCE match
as.data.frame(info[is.na(info$jhpce_path), ])
#   Section Slide Round Combination Average Copy Intensity for AKT3 from initial analysis
# 1   6432M  2, B     1        Star                                                  1.82
# 2   6432M  2, A     1      Circle                                                 10.97
# 3   6471M  4, C     2      Circle                                                 23.90
#                                            Comments/Issues Excluded Re HALO completed?
# 1 High Oligo/low SLC; EXTREME GLARE= calls everything OLIG      Yes           Excluded
# 2                                                     <NA>      Yes           Excluded
# 3                           Pretty section, Figure worthy?       No           Complete
#   Confidence in overall quality Average AKT3 Copy Intensity from initial analysis
# 1                      Excluded                                                NA
# 2                      Excluded                                                NA
# 3                          High                                                NA
#                                                         HALO Folder Path (on Neurolucida)
# 1 KDM\\DECON_RIF_Batch1_Round1_RERUN_AKT3\\STAR\\RESCAN_Decon_RIF_6432M_B1R1_BO_AKT3_Star
# 2  KDM\\DECON_RIF_Batch1_Round1_RERUN_AKT3\\CIRCLE\\Decon_RIF_6432_Mid_B1_BO_RERUN_Circle
# 3         KDM\\DECON_RIF_Batch1_round2_AKT3_RERUN\\Reunmix_RIF_Batch1_Round2_AKT3\\Circle
#                      HALO Algorithm used for "Final" analysis
# 1 Decon_AKT3_Star_Cyto_Nuc_Copyintensity_SEVEREGLARE_FINAL_R1
# 2      Decon_AKT3_Circle_Cyto_Nuc_Copyintensity_FINAL_R1_SMC2
# 3      Decon_AKT3_Circle_Cyto_Nuc_Copyintensity_FINAL_R2_SMC3
#                                                                                     Path to "Final" csv File
# 1            Z:\\Kelsey\\Deconvolution HALO Analysis\\STAR\\R1_HA_Final_Analysis\\HA_R1_6432M_Star_Final.csv
# 2 Z:\\Kelsey\\Deconvolution HALO Analysis\\CIRCLE\\R1_HA_Circle_Final_Analysis\\HA_R1_6432M_Circle_Final.csv
# 3        Z:\\Kelsey\\Deconvolution HALO Analysis\\CIRCLE\\R2_HA_Final_Analysis\\HA_R2_6471M_Circle_Final.csv
#                                             File Path to Raw Image    ...15
# 1           Z:\\Kelsey\\Polaris\\3422_DECON_RIF_STAR_BOR1S2\\Scan1     <NA>
# 2 Z:\\Kelsey\\Polaris\\02242022_RIF_DECON_B1_R1_BO_Circle_L\\Scan1     <NA>
# 3         Z:\\Kelsey\\Polaris\\031822_RIF_DECON_B1_R2_Rerun_Circle Use redo
#                                                                           Path to Fused/Stitched Images
# 1                  \\\\10.17.9.46\\Neural_Plasticity\\Kelsey\\Polaris\\3422_DECON_RIF_STAR_BOR1S2\\HALO
# 2 \\\\10.17.9.46\\Neural_Plasticity\\Kelsey\\Polaris\\02242022_RIF_DECON_B1_R1_BO_Circle_L\\Scan1\\HALO
# 3    \\\\10.17.9.46\\Neural_Plasticity\\Kelsey\\Polaris\\031822_RIF_DECON_B1_R2_Rerun_Circle\\REDO_HALO
#   jhpce_path
# 1       <NA>
# 2       <NA>
# 3       <NA>

## JHPCE paths that were not matched to
jhpce_found <- unlist(sapply(potential_matches, which))
## Check that none were duplicated
stopifnot(!any(duplicated(jhpce_found)))
jhpce[!seq_len(length(jhpce)) %in% jhpce_found]
# [1] "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_BO_Star/031822_RIF_DECON_B2_R3_BO_Star/REDO_HALO/HALO_round3_BO_AKT3_Star_Fused.tif"
# [2] "raw-data/RNAscope_images/031822_RIF_DECON_B1_R2_Rerun_Circle/031822_RIF_DECON_B1_R2_Rerun_Circle/HALO/HALO_section3_Br6471P_AKT3_Circle_Fused.tif"
# [3] "raw-data/RNAscope_images/031822_RIF_DECON_B2_R3_BO_Circle/031822_RIF_DECON_B2_R3_BO_Circle/REDO_HALO/HALO_Round3_BO_AKT3_Circle_Fused.tif"

## To print the info so I can add it to Google Sheets
cat(paste0(info$jhpce_path, "\n"))
