#   This script was run interactively to investigate whether MuSiC involved
#   randomness (and thus must be run with a seed for reproducibility). It does
#   not involve randomness.

library(MuSiC)
library(SingleCellExperiment)

GSE50244.bulk.eset = readRDS(url('https://xuranw.github.io/MuSiC/data/GSE50244bulkeset.rds'))
EMTAB.sce = readRDS(url('https://xuranw.github.io/MuSiC/data/EMTABsce_healthy.rds'))
XinT2D.sce = readRDS(url('https://xuranw.github.io/MuSiC/data/XinT2Dsce.rds'))

bulk.mtx = exprs(GSE50244.bulk.eset)
bulk.meta = exprs(GSE50244.bulk.eset)

set.seed(0)
est_prop1 = music_prop(
    bulk.mtx = bulk.mtx, sc.sce = EMTAB.sce, clusters = 'cellType',
    samples = 'sampleID', verbose = FALSE,
    select.ct = c('alpha', 'beta', 'delta', 'gamma', 'acinar', 'ductal')
)
set.seed(50)
est_prop2 = music_prop(
    bulk.mtx = bulk.mtx, sc.sce = EMTAB.sce, clusters = 'cellType',
    samples = 'sampleID', verbose = FALSE,
    select.ct = c('alpha', 'beta', 'delta', 'gamma', 'acinar', 'ductal')
)

#   TRUE!
identical(est_prop1, est_prop2)
