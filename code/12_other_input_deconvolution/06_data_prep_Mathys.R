
library("SingleCellExperiment")
library("Matrix")
library("dplyr")
library("here")
library("sessioninfo")
library("scater")
library("org.Hs.eg.db")
#### load Mathy's data #####
## read in data

mathys_dir <- "/dcs05/lieber/marmaypag/legacySingleCell_Tran_Maynard_Neuron_2021_LIBD001/Mathys/"
list.files(mathys_dir)

pd = read.csv(paste0(mathys_dir, "snRNAseqPFC_BA10_biospecimen_metadata.csv"), as.is=TRUE)
pd

dx_info <- read.csv("/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/mathys/ROSMAP_Clinical_2019-05_v3.csv") |>
  mutate(Dx = factor(ifelse(age_first_ad_dx == "", "Control", "AD"), levels = c("Control", "AD")))

pd2 <- pd |> left_join(dx_info)
table(pd2$Dx)
# Control      AD 
# 31      17 

pheno = read.delim(paste0(mathys_dir, "filtered_column_metadata.txt"), row.names = 1) |>
  mutate(cellType_broad = factor(case_when(grepl("Ast", broad.cell.type) ~ "Astro",
                                    grepl("Ex", broad.cell.type) ~ "Excit",
                                    grepl("In", broad.cell.type) ~ "Inhib",
                                    grepl("Mic", broad.cell.type) ~ "Micro",
                                    grepl("Oli", broad.cell.type) ~ "Oligo",
                                    grepl("Opc", broad.cell.type) ~ "OPC",
                                    grepl("Per|End", broad.cell.type) ~ "EndoMural",
                                    TRUE ~ "Other"), 
                                 levels = c('Astro',"EndoMural", "Micro", "Oligo", "OPC", "Excit", "Inhib"))) |>
  left_join(pd2 |> dplyr::select(projid, specimenID, tissue, BrodmannArea, individualID, Dx, cogdx)) 

dat = Matrix::readMM(paste0(mathys_dir, "filtered_count_matrix.mtx"))
genes = read.delim(paste0(mathys_dir, "filtered_gene_row_names.txt"),header=FALSE,as.is=TRUE)

colnames(genes) <- "Symbol"


ensembl <- mapIds(org.Hs.eg.db, keys = genes$Symbol, keytype = "SYMBOL", column = "ENSEMBL")
length(ensembl) == nrow(genes)
table(is.na(ensembl))
# FALSE  TRUE 
# 16549  1377 

genes$gene_id <- ensembl

## add names
rownames(dat) = genes$V1
colnames(dat) = rownames(pheno)

## bulk sce
sce <- SingleCellExperiment(
  assays = list(counts = dat),
  colData = pheno,
  rowData = genes
)

table(sce$cellType_broad)
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib 
# 3392       288      1920     18235      2627     34976      9196

## add logcounts
sce <- logNormCounts(sce)

## filter to genes with ensemblIDs
sce <- sce[!is.na(rowData(sce)$gene_id), ]
rownames(sce) <- rowData(sce)$gene_id

## save data as sce for downstream compatiabilty 
save(sce, file = here("processed-data", "12_other_input_deconvolution", "sce_Mathys.Rdata"))

## cell type prop
ct_prop <- colData(sce) |>
  as.data.frame() |>
  dplyr::group_by(individualID, cellType_broad, Dx) |>
  dplyr::count() |>
  dplyr::group_by(individualID, Dx) |>
  dplyr::mutate(prop = n/sum(n))

write.csv(ct_prop, file = here("processed-data", "12_other_input_deconvolution", "Mathys_ct_prop.csv"))

# slurmjobs::job_single('06_data_prep_Mathys', create_shell = TRUE, memory = '10G', command = "Rscript 06_data_prep_Mathys.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
