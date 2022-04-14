library("tidyverse") 
library("here")
library("sessioninfo")

halo_files <- list.files(here("raw-data","HALO"), full.names = TRUE)
names(halo_files) <- gsub(".csv", "", basename(halo_files))


halo_anno <- map2(halo_files, names(halo_files), function(fn, samp){
  output <- read_csv(fn) %>%
    mutate(cell_type = case_when(OLIG2 == 1 ~ "Oligo",
                                 SLC17A7 == 1 ~ "Excit",
                                 TMEM119 == 1 ~ "Micro",
                                 TRUE~"Other"),
           Sample_name = samp) %>%
    separate(Sample_name , into = c(NA, "Round", "Sample", "Combo", NA), sep = "_")
    ## maybe add combo to "Other" label?
  return(output)
})

halo_all <- do.call("rbind", halo_anno)

(sample_n <- halo_all %>% count(Sample, Combo) %>% rename(sample_n = n))
# Sample Combo sample_n
# <chr>  <chr>    <int>
# 1 2720M  Star     22677
# 2 6432A  Star     54921
# 3 6432M  Star     21017
# 4 6432P  Star     31633
# 5 6471A  Star     46238

sample_prop <- halo_all %>% 
  group_by(Sample, cell_type) %>%
  count() %>%
  left_join(sample_n) %>%
  mutate(prop = n/sample_n)

(sample_prop_wide <- sample_prop %>%
  select(Sample, cell_type, prop) %>%
  pivot_wider(names_from = cell_type, values_from = prop))
# Sample    Excit    Micro Oligo Other
# <chr>     <dbl>    <dbl> <dbl> <dbl>
# 1 2720M   0.280   0.0386   0.125 0.557
# 2 6432A   0.00275 0.00120  0.674 0.322
# 3 6432M  NA       0.000666 0.868 0.131
# 4 6432P   0.116   0.0602   0.269 0.555
# 5 6471A   0.156   0.0319   0.257 0.555

