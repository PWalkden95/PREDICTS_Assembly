rm(list = ls())

require(tidyverse)

### load in refined dataset with combined studies 

PREDICTS <- readRDS("Outputs/refined_predicts.rds")
traits <- readRDS("Outputs/assembly_trait_scores.rds")




rao_data <- PREDICTS %>%  
  
    #### group by Site and Species get abundance if there are some sites that the same species is recorded multiple times 
  dplyr::group_by(SSBS,Jetz_Name) %>% dplyr::mutate(SpeciesSiteAbundance = sum(Effort_Corrected_Measurement), n_spp = n()) %>%
  
  filter(!duplicated(n_spp) | n_spp == 1) %>%
  
  ungroup() %>%   droplevels() %>%  
  
  group_by(SSBS) %>% dplyr::mutate(TotalSiteAbundance = sum(SpeciesSiteAbundance), site_spp = n_distinct(Jetz_Name)) %>% ungroup () %>%
  
  dplyr::filter(site_spp > 1) %>%
  
    ### relative abundance of each species at each site SpeciesSitelevel abundance/TotalSite abundance
  dplyr::mutate(RelativeAbundance = SpeciesSiteAbundance/TotalSiteAbundance) %>%
  
  droplevels()



lala <- test %>% group_by(Block) %>% mutate(species_num = n_distinct(Jetz_Name)) 

test <- rao_data %>% filter(grepl(Source_ID, pattern = "Lasky"))
as.character(unique(rao_data$Source_ID))



wm<-map_data("world") %>% filter(region != "Antartica" ) %>% fortify()

## site coords

site_points <- test %>% distinct(SSBS,Longitude,Latitude)

# generate and plot map

site_plot<-ggplot()+ coord_fixed()+
  geom_map(data =wm, map = wm, 
           aes(group = group, map_id= region),
           fill = "darkgrey")+
  geom_point(data = fortify(site_points), aes(Longitude, Latitude),
             colour = "blue", size = 1)+
  theme_classic()

plot(site_plot)



## 85 studies 

studies <- distinct(PREDICTS,SS, .keep_all = FALSE) %>% pull %>% as.character

study <- studies[1]
### pull species 

source("../Functional_Intactness_Index/Functions/Site_Rao_Q.R")

Rao_test <- Rao_Q_Func_bias(rao_data, traits = traits[["both"]]) 


Rao <- data.frame(Site = Rao_test$call , Bias = as.numeric(Rao_test$FunRao))


hist(Rao)
