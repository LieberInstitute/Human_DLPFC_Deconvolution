library("tidyverse") 
library("here")
library("sessioninfo")

halo_files_prelim <- list.files(path = here("raw-data","HALO","prelim"), full.names = TRUE, recursive = TRUE)
halo_files <- list.files(path = here("raw-data","HALO","Deconvolution_HALO_analysis"), full.names = TRUE, recursive = TRUE)
names(halo_files) <- gsub("HA_|_Final\\.csv","",basename(halo_files))

## All prelim data in final analysis csv
all(basename(halo_files_prelim) %in% basename(halo_files))
# [1] TRUE

kelsey_notes <- read.csv(here("processed-data", "03_HALO","Deconvolution_HALO_Analysis_kelsey_googlesheet.csv"))

slide_tab <- kelsey_notes %>% 
  select(Section, Slide, Combo = Combination) %>%
  separate(Slide, into = c("Slide","subslide"), sep = ', ') %>%
  group_by(Slide, subslide, Combo) %>%
  mutate(i = row_number(),
         subslide2 = factor(paste0(subslide, i)),
         Slide = factor(Slide)) %>%
  ungroup()
  
slide_tab %>% count(Slide, Combo)
slide_tab %>% count(Combo, Slide, subslide2)

meta_data <- tibble(file = names(halo_files)) %>%
  separate(file, into = c("Round","Section", "Combo"), sep = "_", remove = FALSE)%>%
  separate(Section, into = c("BrNum", "Position"), sep = "(?<=[0-9])(?=[A-Z])", remove = FALSE) %>%
  mutate(BrNum = paste0("Br", BrNum)) %>%
  left_join(slide_tab)


meta_data %>% count(Round)
meta_data %>% count(Combo)
# Combo     n
# <chr>      <int>
# 1 Circle        21
# 2 Star          21
meta_data %>% count(BrNum)
# BrNum      n
# <chr>  <int>
# 1 Br2720     4
# 2 Br2723     2
# 3 Br3942     6
# 4 Br6423     4
# 5 Br6432     6
# 6 Br6471     4
# 7 Br6522     4
# 8 Br8325     4
# 9 Br8492     4
# 10 Br8667    4

meta_data %>% count(Position)
# Position     n
# <chr>    <int>
# 1 A           14
# 2 M           16
# 3 P           12

meta_data %>% count(Combo, Slide, subslide2) %>% count(n)

#### Read csv files ####

halo_tables <- map(halo_files, read.csv)
map_int(halo_tables, nrow)

colnames(halo_tables$R1_2720M_Circle)
# [1] "Image.Location"                      "Analysis.Region"                     "Algorithm.Name"                     
# [4] "Object.Id"                           "XMin"                                "XMax"                               
# [7] "YMin"                                "YMax"                                "DAPI.AKT3"                          
# [10] "GAD1.AKT3"                           "GFAP.AKT3"                           "CLDN5.AKT3"                         
# [13] "GAD1"                                "GFAP"                                "CLDN5" 
colnames(halo_tables$R1_2720M_Star)
# [1] "Image.Location"                         "Analysis.Region"                       
# [3] "Algorithm.Name"                         "Object.Id"                             
# [5] "XMin"                                   "XMax"                                  
# [7] "YMin"                                   "YMax"                                  
# [9] "DAPI.AKT3"                              "OLIG2.AKT3"                            
# [11] "SLC17A7.AKT3"                           "TMEM119.AKT3"                          
# [13] "OLIG2"                                  "SLC17A7"                               
# [15] "TMEM119"                                "SLC17A7..Opal.520..Positive"     

tibble(cellType = c('Endo',"Astro", "Inhib", "Excit", "Micro", "Oligo"),
       marker = c('CLDN5','GFAP','GAD1','SLC17A7','TMEM119','OLIG2'),
       Combo = rep(c("Circle", "Star"), each = 3))

# cellType marker  Combo 
# <chr>    <chr>   <chr> 
# 1 Endo     CLDN5   Circle
# 2 Astro    GFAP    Circle
# 3 Inhib    GAD1    Circle
# 4 Excit    SLC17A7 Star  
# 5 Micro    TMEM119 Star  
# 6 Oligo    OLIG2   Star

meta_data <- meta_data %>% add_column(n_nuc = map_int(halo_tables, nrow))
sum(meta_data$n_nuc)
# [1] 1936064write_csv(meta_data, file = here("processed-data","HALO","HALO_meta_data.csv"))

## Exploratory n nuc plots
summary(meta_data$n_nuc)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 22677   40586   46286   46097   53010   63293 

n_nuc_col <- ggplot(meta_data, aes(x = Section, y = n_nuc, fill = Combo)) +
  geom_col(position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(n_nuc_col, filename = here("plots", "03_HALO","n_nuc_col.png"), width = 10)


n_nuc_box <- ggplot(meta_data, aes(x = Slide, y = n_nuc, fill = Combo)) +
  geom_boxplot(position = "dodge") 

ggsave(n_nuc_box, filename = here("plots", "03_HALO","n_nuc_box.png"))


n_nuc_box_round <- ggplot(meta_data, aes(x = Combo, y = n_nuc, fill = Combo)) +
  geom_boxplot(position = "dodge") +
  facet_wrap(~round)

ggsave(n_nuc_box_round, filename = here("plots", "HALO","n_nuc_box_round.png"))


slide_tile_nuc <-  meta_data %>%
  ggplot(aes(x = Slide, y = subslide2, fill = n_nuc)) +
  geom_tile() +
  geom_text(aes(label = Section))+
  facet_wrap(~Combo) +
  theme_bw()

ggsave(slide_tile_nuc, filename = here("plots", "03_HALO","n_nuc_slide_tile.png"), width = 10)


slide_tile_round <-  meta_data %>%
  ggplot(aes(x = Slide, y = subslide2, fill = Round)) +
  geom_tile() +
  geom_text(aes(label = Section))+
  facet_wrap(~Combo) +
  theme_bw()

ggsave(slide_tile_round, filename = here("plots", "03_HALO","round_slide_tile.png"), width = 10)


