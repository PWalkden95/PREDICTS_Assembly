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
  
    ### relative abundance of each species at each site SpeciesSitelevel abundance/TotalSite abundance
  dplyr::mutate(RelativeAbundance = SpeciesSiteAbundance/TotalSiteAbundance) %>%
  
  droplevels()


test <- rao_data %>% filter(grepl(Source_ID, pattern = "Ndanganga"))
as.character(unique(rao_data$Source_ID))

## 85 studies 

studies <- distinct(PREDICTS,SS, .keep_all = FALSE) %>% pull %>% as.character

study <- studies[1]
### pull species 

for(study in studies){
  

  
}