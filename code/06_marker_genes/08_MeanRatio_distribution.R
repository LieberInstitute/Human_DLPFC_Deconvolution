
# library("SummarizedExperiment")
library("tidyverse")
library("here")
library("sessioninfo")
# library("jaffelab")

#### Load data ####
## marker_stats
load(here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"), verbose = TRUE)

marker_stats |>
  filter(ratio > 1) |>
  group_by(cellType.target) |>
  summarize(n1 = n(), 
            n2 = sum(ratio > 2),
            median = median(ratio),
            mad = mad(ratio),
            mad3 = median + 3*mad,
            n_mad = sum(ratio > mad3),
            sd = sd(ratio),
            max = max(ratio))

# cellType.target     n1      n2 median   mad  mad3 n_mad    sd   max
# 1 Astro             206    61   1.41 0.470  2.82    37  2.04   13.9
# 2 EndoMural         262    76   1.35 0.455  2.71    51  5.97   62.4
# 3 Micro             260   100   1.69 0.831  4.18    48 12.1    87.2
# 4 Oligo             287    81   1.42 0.497  2.91    50  3.35   21.6
# 5 OPC               392    76   1.23 0.281  2.08    68  3.64   47.1
# 6 Excit            1664   145   1.22 0.237  1.94   160  0.944  16.1
# 7 Inhib            1460    59   1.15 0.137  1.56   126  0.979  18.1

marker_stats |>
  filter(ratio > 1) |>
  ggplot(aes(x = ratio)) +
  geom_density() +
  facet_wrap(~cellType.target)

marker_stats |>
  filter(ratio > 2) |>
  ggplot(aes(x = ratio, y = cellType.target)) +
  geom_boxplot()
