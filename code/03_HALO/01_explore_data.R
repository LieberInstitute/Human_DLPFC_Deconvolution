library("tidyverse") 
library("scales")
library("here")
library("sessioninfo")
library("DeconvoBuddies")

#### Plot Set-up ####
plot_dir <- here("plots", "03_HALO",  "01_explore_data")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

load(here("processed-data","00_data_prep","cell_colors.Rdata"), verbose = TRUE)
# cell_type_colors_halo

## Halo files 
halo_files_prelim <- list.files(path = here("raw-data","HALO","prelim"), full.names = TRUE, recursive = TRUE)
halo_files <- list.files(path = here("raw-data","HALO","Deconvolution_HALO_analysis"), full.names = TRUE, recursive = TRUE)
names(halo_files) <- gsub("HA_|_Final\\.csv","",basename(halo_files))

## All prelim data in final analysis csv
all(basename(halo_files_prelim) %in% basename(halo_files))
# [1] TRUE

kelsey_notes <- read.csv(here("processed-data", "03_HALO","Deconvolution_HALO_Analysis_kelsey_googlesheet.csv")) %>%
  mutate(Glare = grepl("GLARE", Comments.Issues))

slide_tab <- kelsey_notes %>% 
  select(Section, Slide, Combo = Combination, Glare) %>%
  separate(Slide, into = c("Slide","subslide"), sep = ', ') %>%
  group_by(Slide, subslide, Combo) %>%
  mutate(i = row_number(),
         subslide2 = factor(paste0(subslide, i)),
         Slide = factor(Slide)) %>%
  ungroup()
  
slide_tab %>% count(Slide, Combo)
slide_tab %>% count(Combo, Slide, subslide2)

## Position Names
position_codes <- data.frame(Pos = c("A","M","P"),
                             Position = factor(c("Anterior","Middle","Posterior")))

## Sample_Combo - maybe need a better name?
metadata <- tibble(Sample_Combo = names(halo_files)) %>%
  separate(Sample_Combo, into = c("Round","Section", "Combo"), sep = "_", remove = FALSE)%>%
  separate(Section, into = c("BrNum", "Pos"), sep = "(?<=[0-9])(?=[A-Z])", remove = FALSE) %>%
  left_join(position_codes) %>%
  mutate(BrNum = paste0("Br", BrNum),
         Sample = paste0(BrNum, "_", tolower(substr(Position, 1, 3)))) %>%
  left_join(slide_tab)

head(metadata)
# Sample_Combo    Round Section BrNum  Pos   Combo  Position  Sample     Slide subslide     i subslide2
# <chr>           <chr> <chr>   <chr>  <chr> <chr>  <fct>     <chr>      <fct> <chr>    <int> <fct>    
#   1 R1_2720M_Circle R1    2720M   Br2720 M     Circle Middle    Br2720_mid 5     A            1 A1       
# 2 R1_6432A_Circle R1    6432A   Br6432 A     Circle Anterior  Br6432_ant 5     D            1 D1       
# 3 R1_6432M_Circle R1    6432M   Br6432 M     Circle Middle    Br6432_mid 2     A            1 A1       
# 4 R1_6432P_Circle R1    6432P   Br6432 P     Circle Posterior Br6432_pos 5     B            1 B1       
# 5 R1_6471A_Circle R1    6471A   Br6471 A     Circle Anterior  Br6471_ant 5     C            1 C1       
# 6 R2_6471M_Circle R2    6471M   Br6471 M     Circle Middle    Br6471_mid 4     C            1 C1  


metadata %>% count(Round)
metadata %>% count(Combo)
# Combo     n
# <chr>      <int>
# 1 Circle        21
# 2 Star          21
metadata %>% count(BrNum)
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

metadata %>% count(Position)
# Position      n
# <fct>     <int>
# 1 Anterior     14
# 2 Middle       16
# 3 Posterior    12

metadata %>% count(Combo, Slide, subslide2)

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


metadata <- metadata %>% add_column(n_nuc = map_int(halo_tables, nrow))
sum(metadata$n_nuc)
# [1] 1936064

write_csv(metadata, file = here("processed-data","03_HALO","HALO_metadata.csv"))

#### Exploratory n nuc plots ####
summary(metadata$n_nuc)
sd(metadata$n_nuc)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 22677   40586   46286   46097   53010   63293 

n_nuc_col <- ggplot(metadata, aes(x = Sample, y = n_nuc, fill = Combo)) +
  geom_col(position = "dodge") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(n_nuc_col, filename = here(plot_dir, "n_nuc_col.png"), width = 10)


n_nuc_box <- ggplot(metadata, aes(x = Slide, y = n_nuc, fill = Combo)) +
  geom_boxplot(position = "dodge") + 
  theme_bw()

ggsave(n_nuc_box, filename = here(plot_dir, "n_nuc_box.png"))

## checkout distribution by combo
n_nuc_box_combo <- ggplot(metadata, aes(x = Combo, y = n_nuc, fill = Combo)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2)+
  theme_bw() +
  theme(legend.position = "None")

ggsave(n_nuc_box_combo, filename = here(plot_dir, "n_nuc_box_combo.png"))

n_nuc_box_round <- ggplot(metadata, aes(x = Round, y = n_nuc, fill = Combo)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(color = Glare), position=position_jitterdodge(jitter.width = 0.2))+
  theme_bw() 

ggsave(n_nuc_box_round, filename = here(plot_dir, "n_nuc_box_round.png"))

n_nuc_box_pos <- ggplot(metadata, aes(x = Position, y = n_nuc, fill = Position)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2)+
  facet_wrap(~Combo) +
  theme_bw() +
  theme(legend.position = "None")

ggsave(n_nuc_box_pos, filename = here(plot_dir, "n_nuc_box_position.png"))

## Explore Slide + Splide postion of Samples
slide_tile_nuc <-  metadata %>%
  ggplot(aes(x = Slide, y = subslide2, fill = n_nuc)) +
  geom_tile() +
  geom_text(aes(label = Sample), size = 2.5)+
  facet_wrap(~Combo) +
  scale_fill_viridis()+
  theme_bw()

ggsave(slide_tile_nuc, filename = here(plot_dir, "n_nuc_slide_tile.png"), width = 10)


slide_tile_round <-  metadata %>%
  ggplot(aes(x = Slide, y = subslide2, fill = Round)) +
  geom_tile() +
  geom_text(aes(label = Sample), size = 2.5)+
  facet_wrap(~Combo) +
  theme_bw()

ggsave(slide_tile_round, filename = here(plot_dir, "round_slide_tile.png"), width = 10)

#### Cell Type Annotation ####

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

## CIRCLE COMBO
halo_circle <- do.call("rbind", halo_tables[grep("Circle", names(halo_tables))])
dim(halo_circle)
# [1] 1011916      45

halo_circle %>% count(GAD1,GFAP,CLDN5)
# GAD1 GFAP CLDN5      n
# 1    0    0     0 683079
# 2    0    0     1  38385
# 3    0    1     0 188200
# 4    1    0     0 102252

## cell cell types
halo_circle <- halo_circle %>% 
  rownames_to_column("Sample") %>%
  separate(Sample, into = c("Sample_Combo",NA), sep = "\\.") %>%
  mutate(n_marker = GAD1 + GFAP + CLDN5,
         cell_type = case_when(
           n_marker > 1 ~ "Multi",
           GAD1 == 1 ~ "Inhib",
           GFAP == 1 ~ "Astro",
           CLDN5 == 1 ~ "Endo",
           TRUE ~ "Other"
         ))

head(halo_circle)

## ct marker copies
table(halo_circle$CLDN5..Alexa.488..Copies == 0)
summary(halo_circle$CLDN5..Alexa.488..Copies)
summary(halo_circle$GFAP..Alexa.594..Copies)

## TREG copies
summary(halo_circle$AKT3..Opal.570..Copies)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   2.000   4.000   6.017   8.000  64.000 


## no GAD1 copies?
head(halo_circle[,grep("Copies", colnames(halo_circle))])
# GFAP..Alexa.594..Copies AKT3..Opal.570..Copies CLDN5..Alexa.488..Copies


halo_circle %>% count(cell_type)
# cell_type      n
# 1     Astro 188200
# 2      Endo  38385
# 3     Inhib 102252
# 4     Other 683079

circle_n_nuc <- halo_circle %>%
  group_by(Sample_Combo) %>%
  summarize(n_nuc = n())

circle_prop <- halo_circle %>% 
  count(Sample_Combo, cell_type) %>%
  right_join(metadata, .) %>% 
  mutate(prop = n/n_nuc)

circle_prop_wide <- circle_prop %>%
  mutate(prop = round(prop, 2)) %>%
  select(Sample_Combo, cell_type, prop) %>%
  pivot_wider(names_from = "cell_type", values_from = "prop") 

#### Star combo ####

## need to fix R5 colnames SLC17A -> SLC17A7
cn <- colnames(halo_tables$R1_2720M_Star)

## fix typos
colnames(halo_tables$R5_8325A_Star) <- gsub("SLC17A\\.","SLC17A7",colnames(halo_tables$R5_8325A_Star))
colnames(halo_tables$R5_8325M_Star) <- gsub("SLC17A\\.","SLC17A7",colnames(halo_tables$R5_8325M_Star))
colnames(halo_tables$R5_8667A_Star) <- gsub("SLC17A\\.","SLC17A7",colnames(halo_tables$R5_8667A_Star))

map_lgl(halo_tables[grep("Star", names(halo_tables))], ~all(colnames(.x) == cn))

## rbind all star tables
halo_star <- do.call("rbind", halo_tables[grep("Star", names(halo_tables))])
dim(halo_star)
halo_star %>% count(SLC17A7,TMEM119,OLIG2)
# SLC17A7 TMEM119 OLIG2      n
# 1       0       0     0 623519
# 2       0       0     1  66579
# 3       0       1     0  34504
# 4       1       0     0 199546

halo_star <- halo_star %>% 
  rownames_to_column("Sample") %>%
  separate(Sample, into = c("Sample_Combo",NA), sep = "\\.") %>%
  mutate(n_marker = SLC17A7 + TMEM119 + OLIG2,
         cell_type = case_when(
           n_marker > 1 ~ "Multi",
           SLC17A7 == 1 ~ "Excit",
           TMEM119 == 1 ~ "Micro",
           OLIG2 == 1 ~ "Oligo",
           TRUE ~ "Other"
         ))

halo_star %>% count(cell_type)
halo_star %>% count(Sample_Combo, cell_type)
# cell_type      n
# 1     Excit 199546
# 2     Micro  34504
# 3     Oligo  66579
# 4     Other 623519

colnames(halo_tables$R5_8325A_Star)
halo_tables$R5_8325A_Star %>% count(SLC17A7,TMEM119,OLIG2)
halo_tables$R5_8325A_Star %>% count(SLC17A7..Opal.520..Positive.Cytoplasm)

halo_star %>% count(Sample_Combo)

star_prop <- halo_star %>% 
  count(Sample_Combo, cell_type) %>%
  right_join(metadata, .) %>% 
  mutate(prop = n/n_nuc)

star_prop_wide <- star_prop %>%
  mutate(prop = round(prop, 2)) %>%
  select(Sample_Combo, cell_type, prop) %>%
  pivot_wider(names_from = "cell_type", values_from = "prop") 

star_prop %>% filter(cell_type == 'Other', prop == 1) %>% select(Sample, Sample_Combo, Round)


save(halo_star, halo_circle, file = here("processed-data","03_HALO","HALO_Data.Rdata"))

#### Composition Plots ####
circle_colors <- cell_type_colors_halo[c("Astro", "Endo", "Inhib", "Other")]
star_colors <- cell_type_colors_halo[c("Excit", "Micro", "Oligo", "Other")]

circle_prop_bar <- plot_composition_bar(circle_prop, sample_col = "Sample",x_col = "Sample") +
  scale_fill_manual(values = circle_colors) +
  labs(title = "Circle Combo") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(circle_prop_bar, filename = here(plot_dir, "prop_bar_circle.png"), width = 12)

star_prop_bar <- plot_composition_bar(star_prop, sample_col = "Sample", x_col = "Sample") +
  scale_fill_manual(values = star_colors) +
  labs(title = "Star Combo") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(star_prop_bar, filename = here(plot_dir, "prop_bar_star.png"), width = 12)

## Combine ##
prop_all <- rbind(circle_prop, star_prop)

prop_boxplots <- prop_all %>%
  ggplot(aes(x = Round, y  = prop, color = cell_type)) +
  geom_boxplot() +
  scale_color_manual(values = cell_type_colors_halo) +
  facet_wrap(~cell_type, scales = "free_y")

ggsave(prop_boxplots, filename = here(plot_dir, "prop_boxplots.png"))

prop_boxplots <- prop_all %>%
  filter(cell_type != "Other") %>%
  ggplot(aes(x = cell_type, y  = prop, fill = cell_type)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) + 
  geom_jitter(width = 0.2, colour="black",pch=21) +
  scale_fill_manual(values = cell_type_colors_halo) +
  # scale_color_manual(values = cell_type_colors_broad) +
  theme_bw() +
  theme(legend.position = "None")+
  labs(title = "RNAscope Cell Type Proportions")

ggsave(prop_boxplots, filename = here(plot_dir, "prop_boxplots.png"))


prop_boxplot_position <- prop_all %>%
  ggplot(aes(x = Position, y = prop, fill = cell_type)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) + 
  geom_jitter(width = 0.2, colour="black",pch=21) +
  scale_fill_manual(values = cell_type_colors_halo) +
  facet_wrap(~cell_type, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(title = "RNAscope Cell Type Proportions")

ggsave(prop_boxplot_position, filename = here(plot_dir, "prop_boxplot_position.png"))


prop_other_adj <- prop_all %>%
  filter(cell_type != "Other") %>%
  group_by(Sample) %>%
  summarize(cell_type = "Other_est", 
            prop = 1 - sum(prop)) 

prop_all_adj <- prop_all %>%
  filter(cell_type != "Other") %>%
  select(Sample, cell_type, prop) %>%
  rbind(prop_other_adj)

prop_all_adj %>% group_by(Sample) %>% summarize(sum(prop))

adj_prop_bar <- plot_composition_bar(prop_all_adj, sample_col = "Sample",x_col = "Sample",min_prop_text = .01) +
  scale_fill_manual(values = cell_type_colors_halo) +
  labs(title = "Estimated Sample Compositions") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(adj_prop_bar, filename = here(plot_dir, "prop_bar_adj.png"), width = 12)


prop_all %>% 
  filter(prop != 1) %>%
  group_by(Combo, cell_type) %>% 
  summarize(min = min(prop),
            mean = mean(prop),
            median = median(prop),
            max = max(prop))

write_csv(prop_all, file = here("processed-data","03_HALO","HALO_cell_type_proportions.csv"))
