
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
cell.type.labes <- sce$cellType_broad_hc
length(cell.type.labels)
# [1] 23793

table(cell.type.labels)

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

# number of genes filtered in each category: 
#   Rb      Mrp other_Rb     chrM   MALAT1     chrX     chrY 
# 89       78     1011       37        1     2464      594 
# A total of  4214  genes from Rb Mrp other_Rb chrM MALAT1 chrX chrY  have been excluded 
# A total of  24343  gene expressed in fewer than  5  cells have been excluded 

plot.bulk.vs.sc(sc.input = sc.dat.filtered,
                 bulk.input = bk.dat,
                 pdf.prefix= here(plot_dir, "BayesPrism_example")
)

## marker finding 
sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered,
                                       gene.type = "protein_coding")

## slow - needs lots of memory
message(Sys.time(), "- get.exp.stat")
diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],# filter genes to reduce memory use
                              cell.type.labels=cell.type.labels,
                              cell.state.labels=cell.state.labels,
                              pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                              cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                              n.cores=1 #number of threads
)

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
  cell.state.labels = cell.state.labels,
  key="tumor",
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
session_info()

