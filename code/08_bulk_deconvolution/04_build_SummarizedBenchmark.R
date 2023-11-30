
library("SummarizedBenchmark")
library("SingleCellExperiemnt")
library("here")

library("MuSiC")
library("BisqueRNA")

#### load data ####
## load bulk data
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

rownames(sce) <- rowData(sce)$gene_id

## find common genes
common_genes <- intersect(rowData(sce)$gene_id, rowData(rse_gene)$ensemblID)
length(common_genes)
# [1] 17804

## load marker gene data
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)
# marker_stats

marker_stats |>
  dplyr::filter(gene %in% common_genes,
                rank_ratio <= 25) |>
  dplyr::count(cellType.target)

# cellType.target     n
# <fct>           <int>
# 1 Astro              24
# 2 EndoMural          24
# 3 Micro              17
# 4 Oligo              25
# 5 OPC                18
# 6 Excit              23
# 7 Inhib              20

markers <- marker_stats |> 
  dplyr::filter(gene %in% common_genes, rank_ratio <= 25) |>
  dplyr::pull(gene)

## expression set 
exp_set_bulk <- ExpressionSet(assayData = assays(rse_gene)$counts[markers,],
                              phenoData=AnnotatedDataFrame(
                                as.data.frame(colData(rse_gene))[c("SAMPLE_ID")]))

exp_set_sce <- ExpressionSet(assayData = as.matrix(assays(sce)$counts[markers,]),
                             phenoData=AnnotatedDataFrame(
                               as.data.frame(colData(sce)[,c("key","Sample","BrNum", "cellType_broad_hc", "cellType_hc")])))

exprs(exp_set_sce)

sage(Sys.time(), " - Bisque Prep")
exp_set_sce <- exp_set_sce[markers,]
zero_cell_filter <- colSums(Biobase::exprs(exp_set_sce)) != 0
message("Exclude ",sum(!zero_cell_filter), " cells")
# Exclude 33 cells
exp_set_sce <- exp_set_sce[,zero_cell_filter]

#### Build Summarized Benchmark #####
deconvo_bench <- BenchDesign(data = list(bulk = rse_gene, 
                                         bulk.eset = exp_set_bulk, 
                                         sc.eset = exp_set_sce))


#### Bisque
deconvo_bench <- addMethod(deconvo_bench,
                           label = "Bisque",
                           func = BisqueRNA::ReferenceBasedDecomposition,
                           params = rlang::quos(bulk.eset = bulk.eset,
                                                sc.eset = sc.eset,
                                                cell.types = "cellType_broad_hc",
                                                subject.names = "Sample",
                                                use.overlap = FALSE),
                           post = function(x){x$bulk.props})

printMethods(deconvo_bench)


deconvo_bench2 <- buildBench(deconvo_bench)

head(assay(deconvo_bench2))
dim(assay(deconvo_bench2))
# Bisque
# [1,] 0.03431589
# [2,] 0.02021306
# [3,] 0.01102649
# [4,] 0.19167913
# [5,] 0.05537558
# [6,] 0.38401795

colData(deconvo_bench2)

