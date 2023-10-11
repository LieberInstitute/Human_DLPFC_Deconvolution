
library("SummarizedExperiment")
library("tidyverse")
library("here")
library("sessioninfo")

#### Set up ####

## dirs
plot_dir <- here("plots", "09_bulk_DE", "01_bulk_data_exploration")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data", "00_data_prep", "bulk_colors.Rdata"), verbose = TRUE)
# library_combo_colors
# library_prep_colors
# library_type_colors

#### load data ####
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

addmargins(table(rse_gene$library_prep, rse_gene$library_type))
#       polyA RiboZeroGold Sum
# Bulk    19           19  38
# Cyto    18           19  37
# Nuc     18           17  35
# Sum     55           55 110

colnames(colData(rse_gene))
# [1] "SAMPLE_ID"                          "Sample"                             "BrNum"                             
# [4] "pos"                                "Position"                           "library_prep"                      
# [7] "library_type"                       "library_combo"                      "batch"                             
# [10] "round"                              "fastq1"                             "fastq2"                            
# [13] "strandness"                         "basic_statistics"                   "per_base_sequence_quality"         
# [16] "per_tile_sequence_quality"          "per_sequence_quality_scores"        "per_base_sequence_content"         
# [19] "per_sequence_gc_content"            "per_base_n_content"                 "sequence_length_distribution"      
# [22] "sequence_duplication_levels"        "overrepresented_sequences"          "adapter_content"                   
# [25] "kmer_content"                       "SeqLength_R1"                       "percentGC_R1"                      
# [28] "phred15.19_R1"                      "phred65.69_R1"                      "phred115.119_R1"                   
# [31] "phred150.151_R1"                    "phredGT30_R1"                       "phredGT35_R1"                      
# [34] "Adapter65.69_R1"                    "Adapter100.104_R1"                  "Adapter140_R1"                     
# [37] "SeqLength_R2"                       "percentGC_R2"                       "phred15.19_R2"                     
# [40] "phred65.69_R2"                      "phred115.119_R2"                    "phred150.151_R2"                   
# [43] "phredGT30_R2"                       "phredGT35_R2"                       "Adapter65.69_R2"                   
# [46] "Adapter100.104_R2"                  "Adapter140_R2"                      "trimmed"                           
# [49] "ERCCsumLogErr"                      "numReads"                           "numMapped"                         
# [52] "numUnmapped"                        "overallMapRate"                     "concordMapRate"                    
# [55] "bamFile"                            "totalMapped"                        "mitoMapped"                        
# [58] "mitoRate"                           "totalAssignedGene"                  "gene_Assigned"                     
# [61] "gene_Unassigned_Unmapped"           "gene_Unassigned_Read_Type"          "gene_Unassigned_Singleton"         
# [64] "gene_Unassigned_MappingQuality"     "gene_Unassigned_Chimera"            "gene_Unassigned_FragmentLength"    
# [67] "gene_Unassigned_Duplicate"          "gene_Unassigned_MultiMapping"       "gene_Unassigned_Secondary"         
# [70] "gene_Unassigned_NonSplit"           "gene_Unassigned_NoFeatures"         "gene_Unassigned_Overlapping_Length"
# [73] "gene_Unassigned_Ambiguity"          "rRNA_rate"                          "sex"                               
# [76] "age"                                "diagnosis"                          "qc_class"

summary(rse_gene$ERCCsumLogErr)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -43.794 -11.156  -5.260  -7.898  -1.577   7.141      23
