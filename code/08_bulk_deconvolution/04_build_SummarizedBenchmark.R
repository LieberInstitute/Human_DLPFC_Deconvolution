
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
                                                use.overlap = FALSE))

printMethods(deconvo_bench)


deconvo_bench2 <- buildBench(deconvo_bench)
