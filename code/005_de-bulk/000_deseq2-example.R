#BiocManager::install("DESeq2")
BiocManager::install("airway")
library(DESeq2)
library(airway)
library(GEOquery)
library(tximport)
library("readr")
library("tximportData")
library(tximeta)


####################
# deseq example data
####################
dds <- makeExampleDESeqDataSet(m=4)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dispersionFunction(dds)

# de
dds <- DESeq(dds)
res <- results(dds)
summary(res)
# out of 985 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results



#############
# airway data
#############
#------------------
# load example data
#------------------
data("airway")
se <- airway

# convert se -> deseqdataset
ddsSE <- DESeqDataSet(se, design = ~ cell + dex)
# check results
dds <- DESeq(ddsSE)
ddsSE <- DESeq(ddsSE)
identical(dds, ddsSE) # TRUE
# inspect results
res1 <- results(dds)
res2 <- results(ddsSE)
identical(res1, res2) # TRUE


dir <- system.file("extdata",package="airway")
geofile <- file.path(dir, "GSE52778_series_matrix.txt")
gse <- getGEO(filename=geofile)
class(gse)
# [1] "ExpressionSet"
# attr(,"package")
# [1] "Biobase"
dim(gse)
# Features  Samples 
#       0       16
#------
# do de
#------
cts <- counts(gse)
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata,
                              design= ~ batch + condition)

###############
# tximport data
###############
#-------
# params
#-------

#------
# paths
#------
dir <- system.file("extdata", package="tximportData")

#-----
# load
#-----
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
# inspect
samples
# pop center                assay    sample experiment       run
# 1 TSI  UNIGE NA20503.1.M_111124_5 ERS185497  ERX163094 ERR188297
# 2 TSI  UNIGE NA20504.1.M_111124_7 ERS185242  ERX162972 ERR188088
# 3 TSI  UNIGE NA20505.1.M_111124_6 ERS185048  ERX163009 ERR188329
# 4 TSI  UNIGE NA20507.1.M_111124_7 ERS185412  ERX163158 ERR188288
# 5 TSI  UNIGE NA20508.1.M_111124_2 ERS185362  ERX163159 ERR188021
# 6 TSI  UNIGE NA20514.1.M_111124_4 ERS185217  ERX163062 ERR188356

#------------
# format data
#------------
samples$condition <- factor(rep(c("A","B"),each=3))
rownames(samples) <- samples$run
samples[,c("pop","center","run","condition")]

# get files
files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
names(files) <- samples$run
samples$files <- files
samples$name <- samples$run
se <- tximeta(samples)

ddsTxi <- DESeqDataSet(se, design = ~ condition)

#------
# do de
#------