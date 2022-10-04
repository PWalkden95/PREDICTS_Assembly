rm(list = ls())

require(tidyverse)
require(gtools)
require(ggrepel)
require(gawdis)

source("Functions/assembly_rao_func.R")



### load in refined dataset with combined studies 

PREDICTS <- readRDS("Outputs/refined_predicts.rds")
traits <- readRDS("Outputs/assembly_trait_scores.rds")
drop_spp <- readRDS("Outputs/assembly_drop_spp.rds")
dist_mat <- readRDS("Outputs/study_species_pool_dist_ninety.rds")


rao_data <- PREDICTS %>% dplyr::filter(!(Birdlife_Name %in% drop_spp)) %>%  
  
    #### group by Site and Species get abundance if there are some sites that the same species is recorded multiple times 
  dplyr::group_by(SSBS,Birdlife_Name) %>% dplyr::mutate(SpeciesSiteAbundance = sum(Effort_Corrected_Measurement), n_spp = n()) %>%
  
  filter(!duplicated(n_spp) | n_spp == 1) %>%
  
  ungroup() %>%   droplevels() %>%  
  
  group_by(SSBS) %>% dplyr::mutate(TotalSiteAbundance = sum(SpeciesSiteAbundance), site_spp = n_distinct(Birdlife_Name)) %>% ungroup () %>%
  
  dplyr::filter(site_spp > 1) %>%
  
    ### relative abundance of each species at each site SpeciesSitelevel abundance/TotalSite abundance
  dplyr::mutate(RelativeAbundance = SpeciesSiteAbundance/TotalSiteAbundance) %>%
  
  droplevels()


site_data <- rao_data %>% dplyr::distinct(SSBS)
site_rao <- Rao_Q_Func_bias(rao_data, traits = traits[["both"]])

site_data$Rao <- site_rao$FunRao


