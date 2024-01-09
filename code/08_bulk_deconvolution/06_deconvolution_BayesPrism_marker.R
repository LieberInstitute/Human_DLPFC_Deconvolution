
# devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")

library("BayesPrism")
library("SingleCellExperiment")
library("here")
library("sessioninfo")

plot_dir <- here("plots" , "08_bulk_deconvolution", "06_deconvolution_BayesPrism_marker")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load DLPFC data ####
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)
dim(rse_gene)
# [1] 21745   110

rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)
load(here("processed-data","08_bulk_deconvolution", "markers_top25.Rdata"), verbose = TRUE)

rownames(sce) <- rowData(sce)$gene_id
sce <- sce[markers_top25, sce$cellType_broad_hc != "Ambiguous"]

# ## to test
# sce <- sce[,sce$cellType_broad_hc %in% c("Oligo", "Excit", "Inhib")]
# sce <- sce[sample(rownames(sce), 2000), sample(colnames(sce), 2000)]

sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)


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

sc.dat.filtered <- cleanup.genes(input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)


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
est_prop_BayesPrisim_marker <- run.prism(prism = myPrism, n.cores=50)

message(Sys.time(), "- Saving")
save(est_prop_BayesPrisim_marker, diff.exp.stat, file = here("processed-data","08_bulk_deconvolution","est_prop_BayesPrisim_marker.Rdata"))


# slurmjobs::job_single(name = "06_deconvolution_BayesPrism_marker", memory = "50G", cores = 1, create_shell = TRUE, command = "Rscript 06_deconvolution_BayesPrism_marker.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_in
