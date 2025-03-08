library("SummarizedExperiment")
library("purrr")
library("dplyr")
library("recount")
library("sessioninfo")
library("here")
library("jaffelab")

## prep dirs ##
plot_dir <- here("plots", "02_quality_control", "04_get_expression_cutoff")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load Data ####
message(Sys.time(), " - Load data")
features <- c("gene", "exon", "jx", "tx")

preQC_paths <- here("processed-data","rse","preQC", paste0("rse_", features, "_preQC.Rdata"))
names(preQC_paths) <- features

rse_list <- lapply(preQC_paths, function(x) get(load(x)))

message("Original Dimesions")
map(rse_list, dim)

#### Add missing colData ####
sample_info <- read.csv(here("processed-data", "00_data_prep", "sample_info.csv"))
sample_info$sample_id <- gsub("_2","",sample_info$sample_id)
rownames(sample_info) <- sample_info$sample_id

rse_samples <- unique(rse_list$gene$Sample)
all(rse_samples %in% sample_info$sample_id)

sample_info_missing <- sample_info[rse_list$gene$Sample, c("sex","age","diagnosis")]

rse_list <- map(rse_list, function(rse){
  colData(rse) <- cbind(colData(rse), sample_info_missing)
  return(rse)
  })

#### drop QC samples ####
message(Sys.time(), " - Drop Samples identified in QC")
qc_tb <- read.csv(file = here("processed-data","02_quality_control","QC_record_DLPFC_bulk.csv")) |>
  dplyr::mutate(drop_pca = SAMPLE_ID == "AN00000906_Br8492_Mid_Nuc")

## three samples to drop
qc_tb |>
  dplyr::filter(qc_class == "drop" | drop_pca)
#                      SAMPLE_ID auto_drop auto_warn qc_class drop_pca
# 1 2107UNHS-0293_Br2720_Mid_Nuc      TRUE     FALSE     drop    FALSE
# 2   AN00000904_Br2743_Ant_Cyto      TRUE      TRUE     drop    FALSE
# 3    AN00000906_Br8492_Mid_Nuc     FALSE      TRUE     warn     TRUE

## save preQC col data
bulk_qc <- as.data.frame(colData(rse_list$gene)) |>
  dplyr::left_join(qc_tb)
write.csv(bulk_qc, here("processed-data", "02_quality_control", "preQC_colData.csv"), row.names = FALSE)


rse_list <- map(rse_list, function(rse){
  stopifnot(identical(qc_tb$SAMPLE_ID, colnames(rse)))
  rse$qc_class <- qc_tb$qc_class 
  
  ## drop 2 poor QC samples
  rse <- rse[,rse$qc_class != "drop"]
  
  ## drop PCA outlier AN00000906_Br8492_Mid_Nuc
  rse <- rse[,rse$SAMPLE_ID != "AN00000906_Br8492_Mid_Nuc"]
  return(rse)
})

message("Number of Samples post-filter")
map_int(rse_list, ncol)

#### check meanExprs ####
#### get expression cutoff ####

## set jxn to length to 100
rowData(rse_list$jx)$Length <- 100 

message(Sys.time(), " - Get RPKM")
exprs <- map(rse_list[c('gene',"exon","jx")], ~recount::getRPKM(.x, "Length"))
exprs[["tx"]] <- assays(rse_list$tx)$tpm

seed <- 20230417
seeds <- seed + 0:3
names(seeds) <- features

message(Sys.time(), " - Get Expression Cutoff")
cutoffs <- map2(exprs, names(exprs), function(expr, type){
  message(toupper(type))
  
  pdf(here(plot_dir, paste0('suggested_expr_cutoffs_', tolower(type), '.pdf')), width = 12)
  
  cuts <- jaffelab::expression_cutoff(expr, seed = seeds[type])
  message(paste(cuts, collapse = ' '))
  cut <- max(cuts)
  
  dev.off()
  
  return(cut)
})

### Filter features in RSEs ####
message(Sys.time(), " - Drop features below cutoff")

# means <- lapply(exprs, rowMeans)
rse_list <- pmap(list(rse = rse_list, cutoff = cutoffs, expr = exprs, n = names(rse_list)), function(rse, cutoff, expr, n){
  means <- rowMeans(expr)
  rowRanges(rse)$meanExprs <- means
  rowRanges(rse)$passExprsCut <- means > cutoff
  message(n, " - Filter out ", sum(!rowRanges(rse)$passExprsCut), " features (", round(sum(!rowRanges(rse)$passExprsCut)/nrow(rse), 2),"%)")
  rse <- rse[rowRanges(rse)$passExprsCut,]
  return(rse)
})

message("Dimensions post-filter")
map(rse_list, dim)

#### Add normalized logcounts & Save ####
message("Object Sizes:")
walk(rse_list, ~print(object.size(.x),units = "auto"))

# walk2(rse_list, names(rse_list), ~{
#   message(Sys.time(), " - Save ", .y)
#   save(rse=.x, file = here("processed-data","rse", paste0("rse_", .y, ".Rdata")))
# })

rse_gene <- rse_list$gene
rse_exon <- rse_list$exon
rse_jx <- rse_list$jx
rse_tx <- rse_list$tx

## add logcounts
message(Sys.time(), " - Calc Logcounts")
assays(rse_gene)$logcounts <- log2(getRPKM(rse_gene, "Length")+1)
assays(rse_exon)$logcounts <- log2(getRPKM(rse_exon, "Length")+1)
assays(rse_jx)$logcounts <- log2(getRPKM(rse_jx, "Length")+1)
assays(rse_tx)$logcounts <- log2(assays(rse_tx)$tpm+1)

message(Sys.time(), " - Save Data")
save(rse_gene, file = here("processed-data","rse", "rse_gene.Rdata"))
save(rse_exon, file = here("processed-data","rse", "rse_exon.Rdata"))
save(rse_jx, file = here("processed-data","rse", "rse_jx.Rdata"))
save(rse_tx, file = here("processed-data","rse", "rse_tx.Rdata"))

# slurmjobs::job_single('04_get_expression_cutoff', create_shell = TRUE, memory = '50G', command = "Rscript 04_get_expression_cutoff.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
