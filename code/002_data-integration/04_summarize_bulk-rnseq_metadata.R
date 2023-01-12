#!/usr/bin/env R

# Author: Sean Maden
#
# Description: 
# 
# Notes:

#----------
# load data
#----------

bulk.dpath <- file.path("dcs04/lieber/lcolladotor/deconvolution_LIBD4030",
                        "Human_DLPFC_Deconvolution/processed-data",
                        "01_SPEAQeasy/")
list.files(bulk.dpath)
# [1] "data_info.csv"         "library_type.csv"      "round1_2021-10-19"
# [4] "round2_v25_2022-07-06" "round2_v40_2022-07-06"

list.files(file.path(bulk.dpath, "round2_v40_2022-07-06"))

data.fpath <- file.path(bulk.dpath, "data_info.csv")
lib.fpath <- file.path(bulk.dpath, "library_type.csv")
rse.fpath <- file.path(bulk.dpath, "round2_v40_2022-07-06",
                       "rse", "rse_gene.Rdata")


# read csvs
lib <- read.csv(lib.fpath)
data <- read.csv(data.fpath)

#------------------
# dataset variables
#------------------
colnames(data)
# [1] "SAMPLE_ID"    "Dataset"      "BrNum"        "location"     "library_prep"
# [6] "library_type" "round"        "fastq1"       "fastq2"

lib
#         Dataset library_type round
# 1 2107UNHS-0291        polyA     1
# 2 2107UNHS-0293 RiboZeroGold     1
# 3    AN00000904        polyA     2
# 4    AN00000906 RiboZeroGold     2

table(data$library_type, data$Dataset, data$round)
# , ,  = 1
#              2107UNHS-0291 2107UNHS-0293 AN00000904 AN00000906
# polyA                   12             0          0          0
# RiboZeroGold             0            12          0          0
# , ,  = 2
#              2107UNHS-0291 2107UNHS-0293 AN00000904 AN00000906
# polyA                    0             0         44          0
# RiboZeroGold             0             0          0         45

#------------
# inspect rse
#------------
list.files(file.path(bulk.dpath, "round2_v40_2022-07-06", "rse"))
rse <- get(load(rse.fpath))
dim(rse)
# [1] 61544   113
class(rse)
# [1] "RangedSummarizedExperiment"
# attr(,"package")
# [1] "SummarizedExperiment"

dim(colData(rse))
# [1] 113  70
class(colData(rse))
# [1] "DFrame"
# attr(,"package")
# [1] "S4Vectors"

metadata(rse)
# list()



#-----------------
# dataset mappings
#-----------------


# 

# compare gene counts
# compare gene count levels

# perform ttest?
# perform other de test?
