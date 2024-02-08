
# devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")

library("BayesPrism")
library("SingleCellExperiment")
library("here")
library("sessioninfo")

plot_dir <- here("plots" , "08_bulk_deconvolution", "06_deconvolution_BayesPrism")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load DLPFC data ####
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]

# ## to test
# sce <- sce[,sce$cellType_broad_hc %in% c("Oligo", "Excit", "Inhib")]
# sce <- sce[sample(rownames(sce), 2000), sample(colnames(sce), 2000)]

sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

rownames(sce) <- rowData(sce)$gene_id

#### convert to bayPrism matricies ####
# sample-by-gene raw count matrix of bulk RNA-seq expression
bk.dat <- t(assays(rse_gene)$counts)

dim(bk.dat)
# [1]   110 21745
bk.dat[1:5,1:5]
class(bk.dat)
# [1] "matrix" "array" 

#cell-by-gene raw count matrix of single cell RNA-seq expression
sc.dat <- t(as.matrix(counts(sce)))
dim(sc.dat)
# [1] 23793 60294
class(sc.dat)
# [1] "matrix" "array"

# cell.type.labels is a character vector of the same length as nrow(sc.dat)
cell.type.labels <- as.character(sce$cellType_broad_hc)
length(cell.type.labels)
# [1] 23793

table(cell.type.labels)
rm(sce)

message(Sys.time(), "- Filter outlier genes")

sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE, #return the data used for plotting. 
  pdf.prefix = here(plot_dir, "BayesPrism_DLPFC") #specify pdf.prefix
)

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE,
  pdf.prefix = here(plot_dir, "BayesPrism_DLPFC_bulk.outlier") #specify pdf.prefix
)

sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)


plot.bulk.vs.sc(sc.input = sc.dat.filtered,
                 bulk.input = bk.dat,
                 pdf.prefix= here(plot_dir, "BayesPrism_DLPFC")
)

## marker finding 
sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered,
                                       gene.type = "protein_coding")

## slow - needs lots of memory
message(Sys.time(), "- get.exp.stat")
diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],# filter genes to reduce memory use
                              cell.type.labels=cell.type.labels,
                              cell.state.labels=cell.type.labels, ## test if this works for just cell type
                              pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                              cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                              n.cores=1 #number of threads
)

## debug get.exp.stat

# cell.type.labels <- as.character(cell.type.labels)
# cell.state.labels <- cell.type.labels
# 
# ct.to.cst <- unique(cbind(cell.type=cell.type.labels, cell.state=cell.state.labels))
# cst.count.table <- table(cell.state.labels)
# cell.count.cutoff=50
# low.count.cst <- names(cst.count.table)[cst.count.table < cell.count.cutoff]
# 
# #normalize ref.dat to prepare input for findMarker	
# lib.size <- rowSums(sc.dat)
# lib.size <- lib.size / median(lib.size)
# dat.tmp <- sc.dat/lib.size
# pseudo.count=0.1
# dat.tmp <- log2(dat.tmp + pseudo.count) - log2(pseudo.count)
# 
# fit.up <- scran::pairwiseTTests(x= t(dat.tmp), 
#                          groups= cell.state.labels, 
#                          direction="up")
# 
# pairs.celltype.first <- ct.to.cst[match(fit.up$pairs$first, ct.to.cst[,"cell.state"]),"cell.type"]
# pairs.celltype.second <- ct.to.cst[match(fit.up$pairs$second, ct.to.cst[,"cell.state"]),"cell.type"]
# 
# filter.idx <- pairs.celltype.first != pairs.celltype.second & ! fit.up$pairs$second %in% low.count.cst
# 
# fit.up[[1]] <- fit.up[[1]][filter.idx]
# fit.up[[2]] <- fit.up[[2]][filter.idx,]
# 
# #get the maxmimum pvalue
# output.up <- scran::combineMarkers(fit.up$statistics, fit.up$pairs, pval.type="all", min.prop=NULL, 
#                             log.p.in=F, log.p.out=F, full.stats=F, pval.field="p.value", 
#                             effect.field="logFC", sorted=F)
# 
# all.ct <- unique(ct.to.cst[,"cell.type"])
# 
# ct.stat.list <- lapply(all.ct,function(ct.i){
#   
#   #subset on the subtypes associated with celltype i (ct.i)
#   cst.i <- ct.to.cst[ct.to.cst[,"cell.type"]==ct.i,"cell.state"]
#   output.up.i <- output.up[cst.i]
#   
#   #take the minimum pvalue over all cst.i
#   pval.up.i <- do.call(cbind,lapply(output.up.i, '[', "p.value"))
#   pval.up.min.i <- apply(pval.up.i,1,min)
#   
#   #take the max lfc over the min lfc of cst.i over cst.j in other ct.j (same as the pvalue=min over "all" type)
#   lfc.i <- apply(do.call(cbind,lapply(output.up.i, function(output.up.i.j) {
#     apply(output.up.i.j[,grepl("logFC",colnames(output.up.i.j)),drop=F],1,min)
#   } )),1,max)
#   
#   data.frame(pval.up.min = pval.up.min.i, 
#              min.lfc = lfc.i)	
# })		
# 
# names(ct.stat.list) <- all.ct

# #get the maxmimum pvalue
# output.up <- combineMarkers(fit.up$statistics, fit.up$pairs, pval.type="all", min.prop=NULL, 
#                             log.p.in=F, log.p.out=F, full.stats=F, pval.field="p.value", 
#                             effect.field="logFC", sorted=F)
# 
# all.ct <- unique(ct.to.cst[,"cell.type"])


## fit markers

sc.dat.filtered.pc.sig <- select.marker(sc.dat=sc.dat.filtered.pc,
                                         stat=diff.exp.stat,
                                         pval.max=0.01,
                                         lfc.min=0.1)


## Construct a prism object
message(Sys.time(), "- Build prism object")
myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.type.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

## run prism 
message(Sys.time(), "- Run Prism")
est_prop_BayesPrisim <- run.prism(prism = myPrism, n.cores=50)

message(Sys.time(), "- Saving")
save(est_prop_BayesPrisim, diff.exp.stat, file = here("processed-data","08_bulk_deconvolution","est_prop_BayesPrisim.Rdata"))


# slurmjobs::job_single(name = "06_deconvolution_BayesPrism", memory = "100G", cores = 1, create_shell = TRUE, command = "Rscript 06_deconvolution_BayesPrism.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_in