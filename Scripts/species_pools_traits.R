#########################
######################## This script is getting all the species that are present within the entire spectrum of the species pools and join their traits
########################  then calculating teh distance matrix using gawdis to equally weight all the traits 




rm(list = ls())

require(tidyverse)
require(gawdis)

PREDICTS <- readRDS("Outputs/refined_predicts.rds")

traits_df <- readRDS("../PREDICTS_Taxonomy/PREDICTS_imputed_BL_traits.rds")
Forage <- readRDS("../PREDICTS_Taxonomy/PREDICTS_imputed_BL_traits_forage.rds")
AVONET <- read.csv("../../Datasets/GBD/AVONET1_BirdLife.csv")
drop_spp <- readRDS("Outputs/assembly_drop_spp.rds")

species_pool <- readRDS("Outputs/predicts_sites_species_pools.rds")


assembly_species <- unique(c(unlist(unlist(species_pool))))

#### need to add the pseudo spps 

pseudospp <- PREDICTS %>% filter(pseudospp) %>% distinct(Birdlife_Name) %>% pull() %>% as.character()

colnames(traits_df)[1] <- "Birdlife_Name"
#### also calculate trait values - mean for morpho traits and the first four axes of the dietary strategy principal coordinate analysis

assembly_traits <- data.frame(Birdlife_Name = c(assembly_species,pseudospp))



assembly_traits <- assembly_traits %>%
  dplyr::left_join(Forage) %>% dplyr::left_join(traits_df) %>% filter(!duplicated(Birdlife_Name))

source("Functions/trait_scores.R")

colnames(assembly_traits)[44:53] <- c("Bill.TotalCulmen", "Bill.Nares", "Bill.Width", "Bill.Depth", "Tarsus.Length", 
                                   "Wing.Chord","Kipps.Distance", "Secondary1","Hand_Wing_Index","Tail.Length")


site_studies <- data.frame(SSBS = names(species_pool)) %>% dplyr::left_join(PREDICTS[,c("SSBS","SS")]) %>% distinct(SSBS,SS)
studies <- site_studies %>% distinct(SS) %>% pull()


####################################################################


assembly_dist_func <- function(study,traits, pool){
  
  list <- species_pool[which(site_studies$SS == study)]
  PRED_sp <- PREDICTS %>% dplyr::filter(SS == study, !(Birdlife_Name %in% drop_spp)) %>% dplyr::distinct(Birdlife_Name) %>% pull()
  
  
  sp <- c()
  for(i in 1:length(list)){
    
  sp <- unique(c(sp,unlist(list[[i]][[pool]]),PRED_sp))
  }
  
  tr_val <- traits %>% dplyr::filter(Birdlife_Name %in% sp)
  
  tr <- trait_scores(tr_val, traits = "morpho")
  
  diet_traits <- c("Birdlife_Name","Invertivore","Aquatic.predator","Vertivore","Scavenger","Nectarivore","Frugivore",       
                   "Granivore","Herbivore_A","Herbivore_T")
  
#  for_traits <- c("Birdlife_Name","IASC","ISS","ISG", "IGA","IGB","IGG","AQGR","AQPE", "AQAI","AQPL","AQSU","AQDI",
 #                 "FAE","FGL","FGR","NAE","NGL","GRA","GRG","HBTA","HBTG","HBAG","HBAS","HBAD","VASC","VAS","VSS","VGA","VGG","SCA","SCG")
    
  
  
  
#  mf_tr <- tr[["morpho_traits"]][["PC_Scores"]] %>% dplyr::left_join(tr_val[,for_traits]) %>%
 #   set_rownames(tr_val[,"Birdlife_Name"]) %>% dplyr::select(-Birdlife_Name)
  
  md_tr <- tr[["morpho_traits"]][["PC_Scores"]] %>% dplyr::left_join(tr_val[,diet_traits]) %>%
    set_rownames(tr_val[,"Birdlife_Name"]) %>% dplyr::select(-Birdlife_Name)

  
  
  
  
  
  
 # mf_dist <- gawdis(mf_tr, groups = c(1,2,3,rep(4,ncol(mf_tr)-3)),fuzzy = TRUE)
  
  md_dist <- gawdis(md_tr, groups = c(1,2,3,rep(4,ncol(md_tr)-3)), fuzzy = c(4), w.type = "optimized")
  
return(md_dist)
  
  }




study_distance_mats_surround <- lapply(studies,assembly_dist_func,traits = assembly_traits,pool = "surround")
names(study_distance_mats_surround) <- studies
write_rds(study_distance_mats_surround, file = "Outputs/study_species_pool_dist_surround.rds")


study_distance_mats_ninety <- lapply(studies,assembly_dist_func, traits = assembly_traits, pool = "0.9")
names(study_distance_mats_ninety) <- studies 
write_rds(study_distance_mats_ninety, file = "Outputs/study_species_pool_dist_ninety.rds")






