#!/usr/bin/env R

#
# Learning DESeq2 analysis and formatting scripts to do evaluation of
# DE results, figure generation, etc.
#

# BiocManager::install("DESeq2")
BiocManager::install("airway")
library(DESeq2)
# indep hyp weighting
library(IHW) 
# example datasets
library(airway)
library(GEOquery)
library(tximport)
library("readr")
library("tximportData")
library(tximeta)
library("ggplot2")

#----------------------
# make some random data
#----------------------
# note: this is meant to simulate incoming bulk rnaseq data
num.cond <- 6; num.rep <- 5
lab <- 
dds <- makeExampleDESeqDataSet(n = 1000, m = num.cond*num.rep)
dds$







####################
# deseq example data
####################
dds <- makeExampleDESeqDataSet(m=4)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dispersionFunction(dds)

# de
dds <- DESeq(dds)

#----------------
# using results()
#----------------
# notes:
# this function does limited behind scenes filtering (mean of norm ct.)
# alpha default = 0.1, can set as arg (e.g. results(dds, alpha = num))
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

#---------------------
# hypothesis weighting
#---------------------
# notes:
# objective to weight independent hypotheses to maximize overall 
# expt power
# independent hypo weighting (IHW) can be used to weight pvalues similar
# to multiple test adj, only weighting by the hypothesis tested
res.ihw <- results(dds, filterFun = ihw)
# Only 1 bin; IHW reduces to Benjamini Hochberg (uniform weights)
metadata(res.ihw)$ihwResult
# ihwResult object with 1000 hypothesis tests 
# Nominal FDR control level: 0.1 
# Split into 1 bins, based on an ordinal covariate

#--------
# ma plot
#--------
# axes: mean norm. counts (all samples) vs. lfc
png("ma_test.png", width = 5, height = 5, units = "in", res = 400)
plotMA(res, ylim = c(-2,2))
dev.off()

#----------------------
# scaled counts by gene
#----------------------
# specify gene
genei <- which.min(res$padj)
# varname of groups to compare
varname <- "condition"
# with deseq2 
plotCounts(dds, gene = genei, intgroup = varname)

# or just get the dfp to plot with ggplot
dfpi <- plotCounts(dds, gene = genei, 
                intgroup="condition", 
                returnData=TRUE)

# make new ggplot
ggplot(dfpi, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#------------------------------
# scaled counts for gene vector
#------------------------------

# get gene vector
resf <- res[!is.na(res$padj),]
resf <- resf[order(resf$padj),]
genev <- rownames(resf)[seq(10)]
# get dfp
dfpv <- do.call(rbind, lapply(genev, function(genei){
  dfpi <- plotCounts(dds, gene = genei, 
                     intgroup="condition", returnData=TRUE)
  dfpi$sample <- rownames(dfpi)
  return(dfpi)
}))

# use pipes to aggregate across genes
dfp %>% group_by(sample) %>% 
  summarize_all(groups = sample, funs(median))


#--------------------------------
# multi-factor experiment designs
#--------------------------------

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