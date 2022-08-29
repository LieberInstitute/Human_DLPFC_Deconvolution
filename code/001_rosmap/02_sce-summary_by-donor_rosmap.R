#!/usr/bin/env R

# Author: Sean Maden
#
# Summarize sce data by donor for ROSMAP sce object.
#
#

#----------
# load data
#----------
sce.fpath <- file.path("rosmap_snrnaseq", "sce_all_rosmap-original.rda")
sce <- get(load(sce.fpath))
# get metadata
cd <- colData(sce)

#------------------------------------------------------
# number of cells by type and braak score, within donor
#------------------------------------------------------
dfp <- do.call(rbind, lapply(unique(cd$donor), function(di){
  cdf <- cd[cd$donor==di,]; dfp <- as.data.frame(table(cdf$broad_class))
  dfp$donor <- di; dfp$braak_score <- unique(cdf$braaksc); dfp
}))
head(dfp)
#   Var1 Freq               donor braak_score
# 1 Astr 1574 MFC-B1-S1-Cdx1-pAD0           3
# 2 Endo  100 MFC-B1-S1-Cdx1-pAD0           3
# 3  Exc 4252 MFC-B1-S1-Cdx1-pAD0           3
# 4  Inh 1147 MFC-B1-S1-Cdx1-pAD0           3
# 5 Micr  417 MFC-B1-S1-Cdx1-pAD0           3
# 6 None    6 MFC-B1-S1-Cdx1-pAD0           3

#-------------------------
# summaries by braak score
#-------------------------
# number of donors by braak score
dfp <- do.call(rbind, lapply(unique(cd$braaksc), 
                             function(bi){
  cdf <- cd[cd$braaksc==bi,]; sidv <- length(unique(cdf$specimenID))
  data.frame(braak_score = bi, donors = sidv)
}))
dfp
#   braak_score donors
# 1           3      7
# 2           1      5
# 3           5      5
# 4           4      6
# 5           2      1

# number of cells by braak score
dfp <- as.data.frame(table(cd$braaksc))
dfp
# Var1  Freq
# 1    1 31006
# 2    2  7043
# 3    3 53150
# 4    4 37046
# 5    5 34522

# number of cells by braak score, type
dfp <- as.data.frame(table(cd$braaksc, cd$broad_class))
dfp
#    Var1 Var2  Freq
# 1     1 Astr  5691
# 2     2 Astr  1415
# 3     3 Astr  9562
# 4     4 Astr  7532
# 5     5 Astr  5878
# 6     1 Endo   622
# 7     2 Endo    86
# 8     3 Endo   557
# 9     4 Endo   344
# 10    5 Endo   379
# 11    1  Exc  9579
# 12    2  Exc  2157
# 13    3  Exc 18196
# 14    4  Exc 15895
# 15    5  Exc 12532
# 16    1  Inh  4051
# 17    2  Inh  1225
# 18    3  Inh  8510
# 19    4  Inh  5678
# 20    5  Inh  5341
# 21    1 Micr   696
# 22    2 Micr   247
# 23    3 Micr  1763
# 24    4 Micr   399
# 25    5 Micr   811
# 26    1 None    68
# 27    2 None     1
# 28    3 None    42
# 29    4 None    43
# 30    5 None    51
# 31    1 Olig  8254
# 32    2 Olig  1640
# 33    3 Olig 11736
# 34    4 Olig  5330
# 35    5 Olig  7630
# 36    1  OPC  1716
# 37    2  OPC   243
# 38    3  OPC  2380
# 39    4  OPC  1645
# 40    5  OPC  1679
# 41    1 Peri   329
# 42    2 Peri    29
# 43    3 Peri   404
# 44    4 Peri   180
# 45    5 Peri   221

#-------------------------------------------
# number of cells by braak score, specimenID
#-------------------------------------------
cd$donor <- cd$specimenID
dfp <- as.data.frame(table(cd$donor, cd$braaksc))
dfp[dfp$Freq>0,]
#                     Var1 Var2 Freq
# 2    MFC-B1-S2-Cdx1-pAD0    1 8492
# 6    MFC-B1-S6-Cdx4-pAD0    1 6256
# 9   MFC-B2-10-Cog1-Path0    1 6100
# 13  MFC-B2-14-Cog4-Path0    1 4777
# 16   MFC-B2-9-Cog1-Path0    1 5381
# 46  MFC-B3-22-Cog4-Path0    2 7043
# 49   MFC-B1-S1-Cdx1-pAD0    3 9632
# 51   MFC-B1-S3-Cdx1-pAD1    3 6417
# 53   MFC-B1-S5-Cdx4-pAD0    3 9233
# 56   MFC-B1-S8-Cdx4-pAD1    3 7813 
# 59  MFC-B2-12-Cog1-Path1    3 5785
# 65  MFC-B3-17-Cog1-Path0    3 6314
# 66  MFC-B3-18-Cog1-Path0    3 7956
# 82  MFC-B2-11-Cog1-Path1    4 5166
# 84  MFC-B2-13-Cog4-Path0    4 7532
# 91  MFC-B3-19-Cog1-Path1    4 5868
# 92  MFC-B3-20-Cog1-Path1    4 5884
# 93  MFC-B3-21-Cog4-Path0    4 6500
# 95  MFC-B3-23-Cog4-Path1    4 6096
# 100  MFC-B1-S4-Cdx1-pAD1    5 6541
# 103  MFC-B1-S7-Cdx4-pAD1    5 6049
# 110 MFC-B2-15-Cog4-Path1    5 6888
# 111 MFC-B2-16-Cog4-Path1    5 9254
# 120 MFC-B3-24-Cog4-Path1    5 5790
