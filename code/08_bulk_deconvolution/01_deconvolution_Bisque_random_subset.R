library("SingleCellExperiment")
library("BisqueRNA")
library("here")
library("sessioninfo")

marker_label <- 'MeanRatio_top25'
marker_file <- here(
    'processed-data', '08_bulk_deconvolution', 'markers_MeanRatio_top25.txt'
)
sce_path = here("processed-data", "sce", "sce_DLPFC.Rdata")
bulk_path = here("processed-data","rse", "rse_gene.Rdata")
out_dir <- here(
    "processed-data", "08_bulk_deconvolution", "01_deconvolution_Bisque"
)

message("Using ", marker_label," marker genes from:", marker_file)
if(!dir.exists(data_dir)) dir.create(data_dir)

#### load data ####
## load bulk data
load(bulk_path, verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce data
load(sce_path, verbose = TRUE)
rownames(sce) <- rowData(sce)$gene_id

#   Drop ambiguous cells
sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

## find common genes
common_genes <- intersect(rowData(sce)$gene_id, rowData(rse_gene)$ensemblID)
length(common_genes)
# [1] 17804

message("Input Markers:")
markers <- scan(marker_file, what="", sep="\n")
if(!all(markers %in% common_genes)) {
    warning(
        "Markers missing from common genes: ",
        paste(setdiff(markers, common_genes), collapse = ", ")
    )
}
markers <- intersect(markers, common_genes)

message(Sys.time(), " - Prep data with ", length(markers), " genes")

#   Subset to markers and cells with at least some gene expression
nonempty_cells = colSums(assays(sce)$counts[markers,]) > 0
sce = sce[markers, nonempty_cells]
message("Excluding ", sum(!nonempty_cells), " zero-expression cells")

#### Build Expression sets ####
message(Sys.time(), " - Prep Bisque Data (bulk data)")
exp_set_bulk <- ExpressionSet(
    assayData = assays(rse_gene)$counts[markers,],
    phenoData = AnnotatedDataFrame(
        as.data.frame(colData(rse_gene))[c("SAMPLE_ID")]
    )
)

exp_set_sce <- ExpressionSet(
    assayData = as.matrix(assays(sce)$counts),
    phenoData = AnnotatedDataFrame(
    as.data.frame(colData(sce)[,c("key","Sample","BrNum", "cellType_broad_hc", "cellType_hc")]))
)

a = Sys.time()
est_prop_bisque <- ReferenceBasedDecomposition(
    bulk.eset = exp_set_bulk,
    sc.eset = exp_set_sce,
    cell.types = "cellType_broad_hc",
    subject.names = "Sample",
    use.overlap = FALSE
)

Sys.time() - a

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
