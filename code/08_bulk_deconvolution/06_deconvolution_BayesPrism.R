
# devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")

library("BayesPrism")
library("here")

#### BayesPrism Example ####

plot_dir <- here("plots" , "08_bulk_deconvolution", "06_deconvolution_BayesPrism")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
load(here("processed-data", "08_bulk_deconvolution","example_data","tutorial.gbm.rdata"), verbose = TRUE)
# [1] "bk.dat"            "cell.state.labels" "cell.type.labels"  "sc.dat"

dim(bk.dat)
# [1]   169 60483
bk.dat[1:5,1:5]
class(bk.dat)
# [1] "matrix" "array" 

#cell-by-gene raw count matrix of single cell RNA-seq expression
dim(sc.dat)
# [1] 23793 60294
class(sc.dat)
# [1] "matrix" "array"

# cell.type.labels is a character vector of the same length as nrow(sc.dat)
length(cell.type.labels)
# [1] 23793

# plot.cor.phi (input=sc.dat, 
#               input.labels=cell.type.labels, 
#               title="cell type correlation",
#               #specify pdf.prefix if need to output to pdf
#               #pdf.prefix="gbm.cor.ct",
#               cexRow=0.5, cexCol=0.5,
# )

sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE, #return the data used for plotting. 
  pdf.prefix = here(plot_dir, "BayesPrism_example") #specify pdf.prefix
)

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE,
  pdf.prefix = here(plot_dir, "BayesPrism_example_bulk.outlier") #specify pdf.prefix
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
diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],# filter genes to reduce memory use
                              cell.type.labels=cell.type.labels,
                              cell.state.labels=cell.state.labels,
                              pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                              cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                              n.cores=1 #number of threads
)


## Construct a prism object
