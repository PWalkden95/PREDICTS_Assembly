rm(list = ls())

require(tidyverse)
require(magrittr)
require(gtools)

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


## after manually checking that the study surveyed the complete community of birds (Some studies marked as surveying Entire communities)
## had occassionally targeted specfic species, unsuitable for commuity assembly analyses. Espc. when working with observational data

drop <- data.frame(study = c("VK1_2011__Zimmerman", "DI1_2004__Naidoo", "JB3_2016__Gardner", "JB3_2017__Kajtoch",
                             "MJ1_2009__Lehouck", "HP1_2010__Bicknell" , "JB3_2019__Leso", "SE2_2010__McCarthy", "SE1_2011__Arbelaez",
                             "HB1_2009__Parry", "DL1_2010__Sosa", "VK1_2007__StLaurent", "DI1_2011__Dawson", "DL1_2013__Bartolommei"),
                   reason= c("songbirds", "songbirds", "group_abundance", "orchard birds",
                             "frugivores","targeted sampling", "targeted sampling", "camera_traps", "Mixed Bird Flocks", "not entire community",
                             "couldnt find paper", "passerines", "proportion of plots", "Owls and nocturnal birds"))

## filter out these studies 

refine_data <- PRED_Assem %>% filter(!(Source_ID %in% drop$study)) %>% filter(SS != "GN1_2010__Hvenegaard 2")

#### Correct for sampling effort 

refine_data <- refine_data %>% 
  
  
  dplyr::group_by(SS) %>% ## group by study 
  
  dplyr::mutate(Max_Sampling_Effort = max(Sampling_effort)) %>%  # max samplin effort in each site
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(Rescaled_Sampling_Effort = Sampling_effort/Max_Sampling_Effort) %>%  ## rescale sampling effort so that the maximum effort with a study is 1
  
  dplyr::mutate(Effort_Corrected_Measurement = ifelse(Diversity_metric_is_effort_sensitive == TRUE,
                                                      Measurement/Rescaled_Sampling_Effort,
                                                      Measurement)) %>%
  
  droplevels()



## Lasky study is split into surveys conducted in the wet and dry season and they should be combined 


#### enter PREDICTS data and Source_ID of studies to be combined 
source("Functions/Combine_studies.R")




studies <- as.character(unique(refine_data$SS))
reduce <- substr(studies,1,(nchar(studies)-2))
duplicate <- reduce[which(duplicated(reduce))]
duplicate <- duplicate[-8]


for(dup in duplicate){
refine_data <- combine_study(refine_data, dup) %>% droplevels()
}



#### Visulise the sites geographically

wm<-map_data("world") %>% filter(region != "Antartica" ) %>% fortify()

## site coords

site_points <- refine_data %>% distinct(SSBS,Longitude,Latitude)

# generate and plot map

site_plot<-ggplot()+ coord_fixed()+
  geom_map(data =wm, map = wm,
           aes(group = group, map_id= region),
           fill = "darkgrey")+
  geom_point(data = fortify(site_points), aes(Longitude, Latitude),
             colour = "blue", size = 1)+
  theme_classic()

plot(site_plot)


## Save 

write_rds(refine_data, "Outputs/refined_predicts.rds")
  

#### also calculate trait values - mean for morpho traits and the first four axes of the dietary strategy principal coordinate analysis

source("Functions/trait_scores.R")

traits <- trait_scores(refine_data)


write_rds(traits, "Outputs/assembly_trait_scores.rds")
