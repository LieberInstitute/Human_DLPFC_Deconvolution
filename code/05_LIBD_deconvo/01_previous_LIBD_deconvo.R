library("tidyverse")
library("sessioninfo")
library("here")

## Set up plotting

#### Compile data ####

deconvo_data_paths <- list(
    BipSeq = "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/deconvolution/est_prop_Bisque.Rdata",
    BSP3 = "/dcl01/lieber/RNAseq/Datasets/BrainSeq_hg38/Phase3/caudate_BS3P/est_prop.Rdata",
    MDDseq = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/deconvolution/data/est_prop_Bisque.Rdata",
    # AANRI = "/dcl01/beer/lieber/kynon_temp/control_ancestry/analysis/deconvolution",
    GTEx = "/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/GTEx_est_prop.Rdata"
    # Suicide = #local
    #  BSP1 =
    #  BSP2 =
    #  Degradation =
)


est_prop <- lapply(deconvo_data_paths, function(x) get(load(x)))
