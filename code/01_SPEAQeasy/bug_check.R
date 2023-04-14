
library(here)
library(SummarizedExperiment)
library(tidyverse)

#### confirm issue in SPEAQeasy output ####
load(here("processed-data","01_SPEAQeasy","round2_v40_2023-04-05","count_objects","rse_gene_Human_DLPFC_Deconvolution_n113.Rdata" ), verbose = TRUE)
# load(here("processed-data","01_SPEAQeasy","round2_v40_2022-07-06","count_objects","rse_gene_Human_DLPFC_Deconvolution_n113.Rdata" ), verbose = TRUE)
# load(here("processed-data","01_SPEAQeasy","round2_v25_2022-07-06","count_objects","rse_gene_Human_DLPFC_Deconvolution_n113.Rdata" ), verbose = TRUE)

## issue exists in v40 + v25 data

rse_cyto <- rse_gene[,grepl("Cyto",rse_gene$SAMPLE_ID)]
colnames(rse_cyto) <- gsub("_Cyto","",colnames(rse_cyto))

assays(rse_cyto)$counts[1:5,1:5]

rse_bulk <- rse_gene[,grepl("Bulk",rse_gene$SAMPLE_ID)]
# rse_bulk <- rse_gene[,!grepl("Cyto|Nuc",rse_gene$SAMPLE_ID)]
colnames(rse_bulk) <- gsub("_Bulk","",colnames(rse_bulk))
identical(colnames(rse_bulk), colnames(rse_cyto))

## counts are identical !!!
identical(assays(rse_bulk)$counts,assays(rse_cyto)$counts)
# [1] TRUE

## Fixed in 2023-04-05 run :)
# [1] FALSE

## colData is not
identical(colData(rse_bulk), colData(rse_cyto))
# [1] FALSE
all.equal(colData(rse_bulk), colData(rse_cyto))
# [1] "Attributes: < Component “listData”: Component “SAMPLE_ID”: 38 string mismatches >"                       
# [2] "Attributes: < Component “listData”: Component “per_base_sequence_content”: 6 string mismatches >"        
# [3] "Attributes: < Component “listData”: Component “per_sequence_gc_content”: 16 string mismatches >"         
# [4] "Attributes: < Component “listData”: Component “sequence_duplication_levels”: 16 string mismatches >"     
# [5] "Attributes: < Component “listData”: Component “overrepresented_sequences”: 19 string mismatches >"       
# [6] "Attributes: < Component “listData”: Component “percentGC_R1”: 36 string mismatches >"                    
# [7] "Attributes: < Component “listData”: Component “Adapter65-69_R1”: Mean relative difference: 0.4898548 >"  
# [8] "Attributes: < Component “listData”: Component “Adapter100-104_R1”: Mean relative difference: 0.5981347 >"
# [9] "Attributes: < Component “listData”: Component “Adapter140_R1”: Mean relative difference: 0.4194887 >"    
# [10] "Attributes: < Component “listData”: Component “percentGC_R2”: 37 string mismatches >"                    
# [11] "Attributes: < Component “listData”: Component “Adapter65-69_R2”: Mean relative difference: 0.441441 >"   
# [12] "Attributes: < Component “listData”: Component “Adapter100-104_R2”: Mean relative difference: 0.6823818 >"
# [13] "Attributes: < Component “listData”: Component “Adapter140_R2”: Mean relative difference: 0.5136472 >"    
# [14] "Attributes: < Component “listData”: Component “ERCCsumLogErr”: Mean relative difference: 0.1919732 >"    
# [15] "Attributes: < Component “listData”: Component “numReads”: Mean relative difference: 0.1190214 >"         
# [16] "Attributes: < Component “listData”: Component “numMapped”: Mean relative difference: 0.1399605 >"        
# [17] "Attributes: < Component “listData”: Component “numUnmapped”: Mean relative difference: 1.955913 >"       
# [18] "Attributes: < Component “listData”: Component “overallMapRate”: Mean relative difference: 0.04354653 >"  
# [19] "Attributes: < Component “listData”: Component “concordMapRate”: Mean relative difference: 0.04213651 >"  
# [20] "Attributes: < Component “listData”: Component “bamFile”: 38 string mismatches >"                         
# [21] "Attributes: < Component “listData”: Component “totalMapped”: Mean relative difference: 0.2318836 >"      
# [22] "Attributes: < Component “listData”: Component “mitoMapped”: Mean relative difference: 0.3180629 >"       
# [23] "Attributes: < Component “listData”: Component “mitoRate”: Mean relative difference: 0.2811975 >"




#### Check manifest ####
manifest_check <- read.table(here("raw-data", "bulkRNA", "samples.manifest"))
head(manifest_check)
any(duplicated(manifest_check$V1))
any(duplicated(manifest_check$V3))
any(duplicated(manifest_check$V5))

manifest_check |>
  filter(grepl("Br6471_Ant", V5))

# "raw-data/bulkRNA/2107UNHS-0291/Br6471_Ant_Cyto/2107UNHS-0291_Br6471_Ant_Cyto_1.fastq.gz"
# "raw-data/bulkRNA/2107UNHS-0291/Br6471_Ant/2107UNHS-0291_Br6471_Ant_1.fastq.gz"

#### Are fastq's identical ? ####
# diff "raw-data/bulkRNA/2107UNHS-0291/Br6471_Ant_Cyto/2107UNHS-0291_Br6471_Ant_Cyto_1.fastq.gz" "raw-data/bulkRNA/2107UNHS-0291/Br6471_Ant/2107UNHS-0291_Br6471_Ant_1.fastq.gz"
# Binary files raw-data/bulkRNA/2107UNHS-0291/Br6471_Ant_Cyto/2107UNHS-0291_Br6471_Ant_Cyto_1.fastq.gz and raw-data/bulkRNA/2107UNHS-0291/Br6471_Ant/2107UNHS-0291_Br6471_Ant_1.fastq.gz differ


## md5sum
# md5sum "raw-data/bulkRNA/2107UNHS-0291/Br6471_Ant/2107UNHS-0291_Br6471_Ant_1.fastq.gz"
# af05a12e5dd96823faf85bc57bfde77a  raw-data/bulkRNA/2107UNHS-0291/Br6471_Ant/2107UNHS-0291_Br6471_Ant_1.fastq.gz
# md5sum "raw-data/bulkRNA/2107UNHS-0291/Br6471_Ant_Cyto/2107UNHS-0291_Br6471_Ant_Cyto_1.fastq.gz"
# 620291afd86798f6b523650c97a131bb  raw-data/bulkRNA/2107UNHS-0291/Br6471_Ant_Cyto/2107UNHS-0291_Br6471_Ant_Cyto_1.fastq.gz

# ls -l raw-data/bulkRNA/2107UNHS-0291/Br6471_Ant*
#   raw-data/bulkRNA/2107UNHS-0291/Br6471_Ant:
#   total 9621454
# -rw-r----- 1 lhuuki lieber_lcolladotor 4830269591 Sep 22  2021 2107UNHS-0291_Br6471_Ant_1.fastq.gz
# -rw-r----- 1 lhuuki lieber_lcolladotor 5028162704 Sep 22  2021 2107UNHS-0291_Br6471_Ant_2.fastq.gz
# -rw-r----- 1 lhuuki lieber_lcolladotor        112 Sep 22  2021 2107UNHS-0291_Br6471_Ant.md5
# 
# raw-data/bulkRNA/2107UNHS-0291/Br6471_Ant_Cyto:
#   total 12484984
# -rw-r----- 1 lhuuki lieber_lcolladotor 6296921725 Sep 22  2021 2107UNHS-0291_Br6471_Ant_Cyto_1.fastq.gz
# -rw-r----- 1 lhuuki lieber_lcolladotor 6495572608 Sep 22  2021 2107UNHS-0291_Br6471_Ant_Cyto_2.fastq.gz
# -rw-r----- 1 lhuuki lieber_lcolladotor        122 Sep 22  2021 2107UNHS-0291_Br6471_Ant_Cyto.md5




