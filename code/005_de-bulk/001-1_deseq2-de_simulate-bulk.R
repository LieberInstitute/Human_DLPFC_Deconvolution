#!/usr/bin/env R

#
# Simulate the bulk RNA-seq dataset and write functions for analysis 
# of differential expression results.
#

library(DESeq2)
library("ggplot2")

#-----------------
# helper functions
#-----------------

random_bulkdata <- function(design.str, num.lib = 2, num.prep = 3, 
                            num.rep = 19, num.genes = 1000, 
                            seed.val = 0, verbose = TRUE){
  # note: this is meant to simulate incoming bulk rnaseq data
  # 
  # design.str: this is the string of the variable for the design of 
  #   the de experiment (e.g."prep + lib + rep" for: `~ prep + lib + rep`).
  # num.lib: this is either polya or ribozero
  # num.prep: either cyt, nuc, or bulk
  # num.rep: number of distinct samples (i.e. samples!=donor)
  # num.genes: num genes to simulate
  #
  #
  set.seed(seed.val); total.sample <- num.rep*num.prep*num.lib
  # manage coldata labels
  vprep <- rep(c("cyt", "nuc", "bulk"), each = num.rep*num.lib)
  vlib <- rep(c("polya", "rrna"), each = num.rep*num.prep)
  vrep <- rep(seq(num.rep), each = num.prep*num.lib)
  vfinal <- paste0(label.prep, "_", label.lib, "_", label.rep)
  # get random data
  num.dat <- length(label.final)
  # make dds with design default
  dds <- makeExampleDESeqDataSet(n = num.genes, m = num.dat)
  # append coldata
  dds$prep <- factor(vprep);dds$lib <- factor(vlib)
  dds$rep <- factor(vrep);dds$label <- factor(vfinal)
  # define design
  design(dds) <- eval(parse(text=paste0("~",design.str)))
  eval.str <- paste0("makeExampleDESeqDataSet(n = num.genes, ",
                     "m = num.dat, design = ",design,")")
  dds <- eval(parse(text = eval.str))
  
  if(verbose){message("evaluating the following for dds: ")}
  
  
  return(dds)
}

get_de_expt <- function(dds){
  # get a list of de experiment objects from deseq2 outputs
  # dds: a DESeqDataSet, inc. design attribute for de experiment
  dei <- DESeq(dds)
  resi <- results(dei)
  # handle se object (versus dds)
  ## convert se -> deseqdataset
  #ddsSE <- DESeqDataSet(se, design = ~ cell + dex)
  ## check results
  #dds <- DESeq(ddsSE)
}

compare_geneset <- function(){
  # compare normalized counts across groups for a set of genes
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
}

#----------------
# do simulated de
#----------------
dds.rand <- random_bulkdata()
# get de results outputs
lde <- list()
