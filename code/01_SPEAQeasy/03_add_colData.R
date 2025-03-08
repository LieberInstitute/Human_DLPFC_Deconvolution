library("SummarizedExperiment")
library("here")
library("tidyverse")
library("jaffelab")
library("sessioninfo")

#### Round 2 version 40 ####
output_dir <- here("processed-data", "01_SPEAQeasy", "round2_v40_2023-04-05")

rse_fn_v40 <- list.files(here(output_dir, "count_objects"),
    pattern = "rse*", full.names = TRUE
)

names(rse_fn_v40) <- ss(basename(rse_fn_v40), "_", 2)
basename(rse_fn_v40)

# [1] "rse_exon_Human_DLPFC_Deconvolution_n113.Rdata" "rse_gene_Human_DLPFC_Deconvolution_n113.Rdata"
# [3] "rse_jx_Human_DLPFC_Deconvolution_n113.Rdata"   "rse_tx_Human_DLPFC_Deconvolution_n113.Rdata"

load(rse_fn_v40[["gene"]], verbose = TRUE)

#### Load sample data ####
pos_df <- data.frame(
    Position = c("Anterior", "Middle", "Posterior"),
    pos = c("ant", "mid", "post")
)

data_info <- read.csv(here("processed-data", "01_SPEAQeasy", "data_info.csv")) |>
    mutate(
        pos = tolower(location),
        Sample = paste0(BrNum, "_", pos),
        library_combo = paste0(library_type, "_", library_prep),
        library_combo2 = paste0(library_type, "_", rna_extract) ##library_combo2 matches language in paper
    ) |>
    left_join(pos_df) |>
    column_to_rownames("SAMPLE_ID") |>
    select(Sample, BrNum, pos, Position, library_prep, library_type, rna_extract, library_combo,library_combo2, batch = Dataset, round, fastq1, fastq2)

#### add to rse objects ####
dim(data_info)
data_info <- data_info[colnames(rse_gene), ]
identical(rownames(data_info), colnames(rse_gene))

colData(rse_gene) <- cbind(data_info, colData(rse_gene)) %>%
    select(SAMPLE_ID, everything()) %>%
    DataFrame()

table(rse_gene$batch)
# # 2107UNHS-0291 2107UNHS-0293    AN00000904    AN00000906
# #            12            12            44            45
#
table(rse_gene$library_type)
# # polyA RiboZeroGold
# # 56           57


#### Add colData ####
add_colData <- function(rse_fn) {
    rse <- get(load(rse_fn))

    data_info <- data_info[colnames(rse), ]
    stopifnot(identical(rownames(data_info), colnames(rse)))

    colData(rse) <- cbind(data_info, colData(rse)) %>%
        select(SAMPLE_ID, everything()) %>%
        DataFrame()

    ## there was no ERCC spike in for round1 so those values are not "real" replace with NA
    rse$ERCCsumLogErr[rse$round == 1] <- NA

    return(rse)
}

#### add colData to all 4 rse objects & save ####

(new_rse_dir <- here("processed-data", "rse", "preQC"))
if (!dir.exists(new_rse_dir)) dir.create(new_rse_dir)

message(Sys.time(), "- Gene")
rse_gene <- add_colData(rse_fn_v40[["gene"]])
dim(rse_gene)
# [1] 61544   113
save(rse_gene, file = here(new_rse_dir, "rse_gene_preQC.Rdata"))

message(Sys.time(), "- Exon")
rse_exon <- add_colData(rse_fn_v40[["exon"]])
dim(rse_exon)
# [1] 657702    113
save(rse_exon, file = here(new_rse_dir, "rse_exon_preQC.Rdata"))

message(Sys.time(), "- Jx")
rse_jx <- add_colData(rse_fn_v40[["jx"]])
dim(rse_jx)
# [1] 5275646     113
save(rse_jx, file = here(new_rse_dir, "rse_jx_preQC.Rdata"))

message(Sys.time(), "- Tx")
rse_tx <- add_colData(rse_fn_v40[["tx"]])
dim(rse_tx)
# [1] 246624    113
save(rse_tx, file = here(new_rse_dir, "rse_tx_preQC.Rdata"))

# slurmjobs::job_single('03_add_colData', create_shell = TRUE, memory = '5G', command = "Rscript 03_add_colData.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
