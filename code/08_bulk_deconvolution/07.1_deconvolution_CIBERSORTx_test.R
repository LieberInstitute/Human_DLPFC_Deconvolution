
library(tidyverse)
library(here)
library(jaffelab)

data_dir <- here("processed-data" , "08_bulk_deconvolution", "07_deconvolution_CIBERSORTx_prep")

list.files(data_dir)

#### Example with sce_test ####
first_line <- readLines(file(here(data_dir, "sce_counts_test.txt"), "r"),n=1)
first_line_top25 <- readLines(file(here(data_dir, "DLPFC_sc_counts_marker25.txt"), "r"),n=1)
first_line_top25
sce_test <- read.delim(here(data_dir, "sce_counts_test.txt"))
dim(sce_test)
corner(sce_test)
table(ss(colnames(sce_test), "\\."))
# Excit Inhib Oligo 
# 1107   443   450

list.files(here(data_dir, "webtool_test"))
# [1] "CIBERSORTx_Job17_cell_type_sourceGEP.txt"                                                                              
# [2] "CIBERSORTx_Job17_sce_counts_test_inferred_phenoclasses.CIBERSORTx_Job17_sce_counts_test_inferred_refsample.bm.K999.txt"
# [3] "CIBERSORTx_Job17_sce_counts_test_inferred_phenoclasses.txt"                                                            
# [4] "CIBERSORTx_Job17_sce_counts_test_inferred_refsample.txt"                                                               
# [5] "CIBERSORTx_Job18_Results.txt" 

signature_matrix <- read.delim(here(data_dir, "webtool_test", "CIBERSORTx_Job17_cell_type_sourceGEP.txt"), row.names = 1)
dim(signature_matrix)
# [1] 1999    9

head(signature_matrix)
#                         X2       X0       X1         X3        X5         X4         X7         X6
# ENSG00000000005     1.0000    1.000    1.000     1.0000     1.000     1.0000     1.0000     1.0000
# ENSG00000001617     1.0000    1.000    1.000     1.0000     1.000     1.0000     1.0000     1.0000
# ENSG00000002746 11870.0623 5884.588 9254.981 10356.6844 11651.760 12658.1506 16374.6176 10759.8148
# ENSG00000004864   858.7974 2307.973 1234.036   636.9737   314.056   533.1481   294.2172   755.7802
# ENSG00000004939     1.0000    1.000    1.000     1.0000     1.000     1.0000     1.0000     1.0000
# ENSG00000004948     1.0000    1.000    1.000     1.0000     1.000     1.0000     1.0000     1.0000

refsample <- read.delim(here(data_dir, "webtool_test","CIBERSORTx_Job17_sce_counts_test_inferred_refsample.txt"), row.names = 1)
dim(refsample)
# [1] 1999   40
corner(refsample)
# X2     X2.1     X2.2     X2.3     X2.4       X0
# ENSG00000229649    0.000   0.0000    0.000    0.000    0.000    0.000
# ENSG00000273360    0.000   0.0000    0.000    0.000    0.000    0.000
# ENSG00000286885    0.000   0.0000    0.000    0.000    0.000    0.000
# ENSG00000255250    0.000   0.0000    0.000    0.000    0.000    0.000
# ENSG00000263923    0.000   0.0000    0.000    0.000    0.000    0.000
# ENSG00000085788 1055.991 863.3275 1240.982 1197.509 1093.333 1031.526

pheno_classes <- scan(here(data_dir, "webtool_test","CIBERSORTx_Job17_sce_counts_test_inferred_phenoclasses.txt"), what="", sep="\t")
table(pheno_classes)
pheno_classes
# 0   1   2   3   4   5   6   7 
# 1  41 281   1   1   1   1   1


#### check example data ####
scRNA_tutorial <- read.delim(here(data_dir, "tutorial_data", "scRNA-Seq_reference_melanoma_Tirosh_Fig2b-d.txt"))
first_line_tutorial <- readLines(file(here(data_dir, "tutorial_data", "scRNA-Seq_reference_melanoma_Tirosh_Fig2b-d.txt"), "r"),n=1)
first_line_tutorial
dim(scRNA_tutorial)
scRNA_tutorial[1:5,1:5]

## 
list.files((here(data_dir, "tutorial_data", "Fig2ab-NSCLC_PBMCs")))
cs_sn_test <- read.delim(here(data_dir, "tutorial_data", "Fig2ab-NSCLC_PBMCs", "Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt"))
corner(cs_sn_test)
# gene T.cells.CD8 T.cells.CD8.1 T.cells.CD8.2 Monocytes Monocytes.1
# 1  RP11.34P13.7           0             0             0         0           0
# 2    AL627309.1           0             0             0         0           0
# 3    AP006222.2           0             0             0         0           0
# 4 RP4.669L17.10           0             0             0         0           0
# 5  RP5.857K21.3           0             0             0         0           0
# 6 RP11.206L10.3           0             0             0         0           0

cs_bulk_test <- read.delim(here(data_dir, "tutorial_data", "Fig2ab-NSCLC_PBMCs", "Fig2b-WholeBlood_RNAseq.txt"))
corner(cs_bulk_test)
# GeneSym W070517001156 W070517001157 W070517001159 W070517001160 W070517001161
# 1 5_8S_rRNA    0.00000000     0.0000000    0.00000000     0.0000000    0.00000000
# 2   5S_rRNA    0.00000000     0.0000000    0.00000000     0.0000000    0.00000000
# 3       7SK    0.00000000     0.0000000    0.00000000     0.0000000    0.00000000
# 4      A1BG    1.52458880     1.1982094    2.28110130     2.5109631    1.75268600
# 5  A1BG-AS1    0.21001956     0.2630730    0.41086500     0.5714840    0.13972470
# 6      A1CF    0.00707868     0.1093625    0.02942712     0.0129127    0.02266773


DLPFC_bulk <- read.delim(here(data_dir, "DLPFC_bulk_counts.txt"))
corner(DLPFC_bulk)
# GeneSymbol X2107UNHS.0291_Br2720_Mid_Bulk X2107UNHS.0291_Br2720_Mid_Cyto X2107UNHS.0291_Br2720_Mid_Nuc
# 1 ENSG00000227232                             42                             72                           128
# 2 ENSG00000278267                              6                              3                            23
# 3 ENSG00000268903                              2                              7                             2
# 4 ENSG00000269981                              1                              1                             3
# 5 ENSG00000279457                            131                            210                           365
# 6 ENSG00000228463                            101                             46                            56
# X2107UNHS.0291_Br6432_Ant_Bulk X2107UNHS.0291_Br6432_Ant_Cyto
# 1                            146                             70
# 2                             24                              5
# 3                             50                              7
# 4                             15                              1
# 5                            354                            193
# 6                            277                             96

DLPFC_sc <- read.delim(here(data_dir, "sce_counts_test.txt"))
corner(DLPFC_sc)
# Excit Oligo Oligo.1 Inhib Excit.1 Inhib.1
# ENSG00000237765     2     0       0     1       1       1
# ENSG00000229649     0     0       0     0       0       0
# ENSG00000273360     0     0       0     0       0       0
# ENSG00000286885     0     0       0     0       1       0
# ENSG00000255250     0     0       0     0       0       0
# ENSG00000263923     0     0       0     0       0       0
