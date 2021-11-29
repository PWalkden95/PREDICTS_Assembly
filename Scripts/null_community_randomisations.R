rm(list = ls())

require(tidyverse)


PREDICTS <- readRDS("Outputs/refined_predicts.rds") 
species_pool <- readRDS("Outputs/predicts_sites_species_pools.rds")
dispersal_prob <- readRDS("Outputs/sp_dispersal_probabilities.rds")
drop_spp <- readRDS("Outputs/assembly_drop_spp.rds")


#############################################
#############################################

prep_data <- PREDICTS %>% dplyr::filter(!(Birdlife_Name %in% drop_spp)) %>%  
  
  #### group by Site and Species get abundance if there are some sites that the same species is recorded multiple times 
  dplyr::group_by(SSBS,Birdlife_Name) %>% dplyr::mutate(SpeciesSiteAbundance = sum(Effort_Corrected_Measurement), n_spp = n()) %>%
  
  filter(!duplicated(n_spp) | n_spp == 1) %>%
  
  ungroup() %>%   droplevels() %>%  
  
  group_by(SSBS) %>% dplyr::mutate(TotalSiteAbundance = sum(SpeciesSiteAbundance), site_spp = n_distinct(Birdlife_Name)) %>% ungroup () %>%
  
  dplyr::filter(site_spp > 1) %>%
  
  ### relative abundance of each species at each site SpeciesSitelevel abundance/TotalSite abundance
  dplyr::mutate(RelativeAbundance = SpeciesSiteAbundance/TotalSiteAbundance) %>%
  
  droplevels()

################################################################################################################################################



randomisation_func <- function(site,pool){
  
  
  data <- prep_data %>% dplyr::filter(SSBS == site) %>% dplyr::select(RelativeAbundance, Birdlife_Name)
  
  
  sp_pool <- species_pool[[site]][[pool]]
  
  dis_p <- dispersal_prob[site,sp_pool]
  
  random_probs <- data.frame(Birdlife_Name = names(dis_p), prob = dis_p)
  
  
  for(i in 1:1000){
    randomisation <- data.frame(randomised_spp = sample(random_probs$Birdlife_Name, replace = FALSE, prob = random_probs$prob, size = nrow(data)))
    colnames(randomisation) <- paste(site,"random", paste(i), sep = "_")
    
    data <- cbind(data,randomisation)
  }
  
  return(data)
}




require(future)
require(future.apply)

plan(multicore(workers = 8))

predicts_sites <- as.character(unique(prep_data$SSBS))

randomisations <- future_lapply(predicts_sites,randomisation_func, pool = "0.9", future.seed = TRUE)


write_rds(x = randomisations, file = "Outputs/site_randomisations.rds")
