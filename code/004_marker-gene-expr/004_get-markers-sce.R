# get the sce comprising just identified marker genes

# load sce
base.path <- paste0("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/",
                    "DLPFC_snRNAseq/processed-data/sce")
sce.fname <- "sce_DLPFC.Rdata"
sce <- get(load(file.path(base.path, sce.fname)))

# load marker genes
# load ratio table
rt.fpath <- "/users/smaden/ratio-table_dlpfc-ro1.rda"
# get marker genes vector
mgv <- rt.fpath


# subset sce

# recast sce

# save recast sce