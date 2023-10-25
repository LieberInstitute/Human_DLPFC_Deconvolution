
library("SummarizedExperiment")
library("edgeR")
library("variancePartition")
library("purrr")
library("here")
library("jaffelab")
library("sessioninfo")

#### Set up ####
## dirs
# plot_dir <- here("plots", "09_bulk_DE", "08_DREAM_library-type")
# if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## dirs
data_dir <- here("processed-data", "09_bulk_DE", "03_DREAM_sn_v_bulk")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load Data ####
## bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
# rse_gene

## just bulk
rse_gene <- rse_gene[,rse_gene$library_prep == "Bulk"]


## pseudobulk sce
load(here("processed-data", "10_bulk_vs_sn_DE", "sce_pb_sample.Rdata"), verbose = TRUE)
# sce_pb_sample

#### create combined colData ####
## rm all NA cols
colData(sce_pb_sample) <- colData(sce_pb_sample)[,!apply(apply(colData(sce_pb_sample), 2, is.na), 2, all)]
colData(sce_pb_sample)$data_type <- "snRNA-seq"
colData(sce_pb_sample)$library_combo <- NA
colData(sce_pb_sample)$library_type <- NA
colData(sce_pb_sample)$library_prep <- NA

colnames(colData(sce_pb_sample))
# [1] "Sample"                 "SAMPLE_ID"              "pos"                    "BrNum"                 
# [5] "round"                  "Position"               "age"                    "sex"                   
# [9] "diagnosis"              "high_mito"              "low_sum"                "low_detected"          
# [13] "discard_auto"           "registration_variable"  "registration_sample_id" "ncells"                
# [17] "data_type" 

colData(rse_gene)$data_type <- "Bulk"

(common_colnames <- intersect(colnames(colData(sce_pb_sample)), colnames(colData(rse_gene))))
# [1] "Sample"        "SAMPLE_ID"     "pos"           "BrNum"         "round"         "Position"      "age"          
# [8] "sex"           "diagnosis"     "library_combo" "library_type"  "library_prep" 

colData(sce_pb_sample) <- colData(sce_pb_sample)[,common_colnames]
colData(rse_gene) <- colData(rse_gene)[,common_colnames]


#### Match rowData ####
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

common_genes <- intersect(rownames(sce_pb_sample), rownames(rse_gene))
length(common_genes)
# [1] 17660

## subset to common genes
rse_gene <- rse_gene[common_genes,]
sce_pb_sample <- sce_pb_sample[common_genes,]

##create combined rse objects
rse_list <- map(splitit(rse_gene$library_type), 
                ~SummarizedExperiment(colData = rbind(colData(rse_gene[,.x]), 
                                                      colData(sce_pb_sample)),
                                      rowData=rowData(rse_gene[,.x]),
                                      assays = list(counts = cbind(assays(rse_gene[,.x])$counts, 
                                                                   assays(sce_pb_sample)$counts),
                                                    logcounts = cbind(assays(rse_gene[,.x])$logcounts, 
                                                                      assays(sce_pb_sample)$logcounts)
                                      )
                ))

map(rse_list, ~table(.x$data_type))
# $polyA
# 
# Bulk snRNA-seq 
# 19        19 
# 
# $RiboZeroGold
# 
# Bulk snRNA-seq 
# 19        19 


#### RUN DE ####
# The variable to be tested must be a fixed effect
form <- ~data_type + (1|BrNum)

DREAM_library_type <- map2(rse_list, names(rse_list), function(rse, prep_name){
  message(Sys.time(), ' - Running DE ', prep_name)
  
  dge = DGEList(assays(rse)$counts)
  dge = calcNormFactors(dge)
  metadata = as.data.frame(colData(rse))
  
  param = SnowParam(4, "SOCK", progressbar=TRUE)
  # estimate weights using linear mixed model of dream
  message(Sys.time(), ' - voomWithDreamWeights')
  vobjDream = voomWithDreamWeights(dge, form, metadata, BPPARAM=param)
  
  # Fit the dream model on each gene
  # By default, uses the Satterthwaite approximation for the hypothesis test
  message(Sys.time(), ' - Dream')
  fitmm = dream(vobjDream, form, metadata)
  message(Sys.time(), ' - eBayes')
  fitmm = eBayes(fitmm)
  
  # Examine design matrix
  head(fitmm$design, 3)
  
  outDE <- topTable( fitmm, coef=2, number=Inf , sort.by = "none")
  
  message(Sys.time(), " - Saving")
  try(write.csv(outDE, file = here(data_dir, paste0("DREAM_data-type_gene_",prep_name,".csv")), row.names = TRUE))
  
  return(outDE)
})

save(DREAM_library_type, file = here(data_dir, "DREAM_data-type.Rdata"))

# slurmjobs::job_single(name = "03_DREAM_sn_v_bulk", memory = "10G", cores = 1, create_shell = TRUE, command = "Rscript 03_DREAM_sn_v_bulk.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
