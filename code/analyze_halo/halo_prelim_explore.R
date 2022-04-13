library("tidyverse") 
library("here")
library("sessioninfo")

halo_files <- list.files(here("raw-data","HALO"), full.names = TRUE)
names(halo_files) <- gsub(".csv", "", basename(halo_files))

halo_test <- read_csv(halo_files[[1]])
colnames(halo_test)

halo_test %>% count(`Analysis Region`)

halo_test %>% filter(`OLIG2 (Alexa 647) Positive` == 1) %>%
  count(`OLIG2 (Alexa 647) Positive Cytoplasm`, 
        `OLIG2 (Alexa 647) Positive Nucleus`, 
        OLIG2)

halo_test %>% count(`TMEM119 (Alexa 555) Copies` > 0, TMEM119)

halo_test %>% count(OLIG2, SLC17A7, TMEM119)
# OLIG2 SLC17A7 TMEM119     n
# <dbl>   <dbl>   <dbl> <int>
# 1     0       0       0 12631 Other/multi
# 2     0       0       1   876  Micorglia
# 3     0       1       0  6341 Excit
# 4     1       0       0  2829 Oligo

halo_test %>% count(`OLIG2 (Alexa 647) Positive`,
                    `SLC17A7 (Opal 520) Positive`,
                    `TMEM119 (Alexa 555) Copies` > 0)
# `OLIG2 (Alexa 647) Positive` `SLC17A7 (Opal 520) Positive` `\`TMEM119 (Alexa 555) Copies\` > 0`     n
# <dbl>                         <dbl> <lgl>                                <int>
#   1                            0                             0 FALSE                                10016 
# 2                            0                             0 TRUE                                   876 * 
# 3                            0                             1 FALSE                                 6341 *
# 4                            0                             1 TRUE                                   300
# 5                            1                             0 FALSE                                 2829 * 
# 6                            1                             0 TRUE                                  1430
# 7                            1                             1 FALSE                                  530
# 8                            1                             1 TRUE                                   355


halo_test %>% select(starts_with("TMEM119")) %>% summary()

halo_test %>%
  ggplot(aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax)) +
  geom_rect() +
  theme_void() +
  coord_equal()

halo_test2 <- halo_test %>%
  mutate(cell_type = case_when(OLIG2 == 1 ~ "Oligo",
                              SLC17A7 == 1 ~ "Excit",
                              TMEM119 == 1 ~ "Micro",
                              TRUE~"Other"))

halo_test2 %>% count(cell_type) %>%
  mutate(prop = n/nrow(halo_test2))

# cell_type     n   prop
# <chr>     <int>  <dbl>
# 1 Excit      6341 0.280 
# 2 Micro       876 0.0386
# 3 Oligo      2829 0.125 
# 4 Other     12631 0.557

dlpfc_cell_types <- read_csv("code/analyze_halo/DLPFC_cell_types.csv")
total_dlpfc <- sum(dlpfc_cell_types$DLPFC)
dlpfc_cell_types %>% mutate(prop = DLPFC/total_dlpfc)
