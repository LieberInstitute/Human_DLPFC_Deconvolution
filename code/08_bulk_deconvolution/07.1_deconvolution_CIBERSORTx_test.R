
library(tidyverse)
library(here)
library(jaffelab)

data_dir <- here("processed-data" , "08_bulk_deconvolution", "07_deconvolution_CIBERSORTx_prep")

#### Example with sce_test ####
first_line <- readLines(file(here(data_dir, "sce_counts_test.txt"), "r"),n=1)
first_line
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
first_line
dim(scRNA_tutorial)
scRNA_tutorial[1:5,1:5]



