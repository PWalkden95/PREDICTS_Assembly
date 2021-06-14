rm(list = ls())

require(tidyverse)
require(magrittr)

####### Load in PREDICTS, traits & Taxonomy data

Forage <- readRDS("../../Datasets/GBD/Trophic_Foraging_Niche.rds")
PREDICTS <- readRDS("../../Datasets/PREDICTS/diversity-2021-02-24-03-32-59.rds")
PREDICTS_Taxonomy <- readRDS("../../Datasets/PREDICTS/PREDICTS_AVES_Updated_Taxonomy.rds")
Jetz_Traits <- read.csv("../../Datasets/GBD/GBD_Jetz_averages_11_Nov_2020.csv")

#### Change colname for simplicity

colnames(Jetz_Traits)[1] <- "Jetz_Name" 

#### Filter PREDICTS to birds that are resolved to the species level, join updated taxonomy to be compatable to Jetz phylogeny, join traits

PREDICTS_Aves <- PREDICTS %>% dplyr::filter(Rank == "Species", Class == "Aves") %>% 
  dplyr::left_join(PREDICTS_Taxonomy[,c("Taxon", "Jetz_Name")], by = "Taxon") %>%
  dplyr::left_join(Forage[,c(34:46)], by = c("Jetz_Name" = "Taxon")) %>%
  dplyr::left_join(Jetz_Traits[,c(1,8:18)], by = "Jetz_Name")


### For the study we would  like to have abundance recorded sites with greater than 1 species & is focussed on birds or surveys the
### entire community to ensure that survey is capturing the bird community to the fullest extent - those focussed on specific species
### will not be suitable for analyses. 

PRED_Assem <- PREDICTS_Aves %>% dplyr:: filter(Diversity_metric_type == "Abundance", Measurement > 0) %>% dplyr::group_by(SS) %>%
  dplyr::mutate(spp_num = n_distinct(Jetz_Name)) %>% dplyr::filter(spp_num > 1, grepl(Study_name, pattern = "bird", ignore.case = TRUE)|Sampling_target == "Entire community") %>%
  droplevels()


drop <- data.frame(study = c("VK1_2011__Zimmerman", "DI1_2004__Naidoo", "JB3_2016__Gardner", "JB3_2017__Kajtoch",
                             "MJ1_2009__Lehouck", "HP1_2010__Bicknell" , "JB3_2019__Leso", "SE2_2010__McCarthy", "SE1_2011__Arbelaez",
                             "HB1_2009__Parry", "DL1_2010__Sosa", "VK1_2007__StLaurent"),
                   reason= c("songbirds", "songbirds", "group_abundance", "orchard birds",
                             "frugivores","targeted sampling", "targeted sampling", "camera_traps", "Mixed Bird Flocks", "not entire community",
                             "couldnt find paper", "passerines"))
