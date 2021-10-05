rm(list = ls())

require(tidyverse)
require(TPD)
require(magrittr)

PREDICTS <- readRDS("Outputs/refined_predicts.rds")
TPD_traits <- readRDS("Outputs/full_morpho_traits_list.rds")
drop_spp <- readRDS("Outputs/assembly_drop_spp.rds")
for_traits <- readRDS("Outputs/predicts_foraging_pcoa.rds")

################### format PREDICTS data for TPD calculations

TPD_data <- data.frame(PREDICTS) %>% filter(!(Birdlife_Name %in% drop_spp)) %>% 
  
  dplyr::group_by(SSBS,Birdlife_Name) %>% dplyr::mutate(SpeciesSiteAbundance = sum(Effort_Corrected_Measurement), n_spp = n()) %>%
  
  filter(!duplicated(n_spp) ) %>%
  
  ### calculating a hypervolume in 3 dimensions with fewer than 21 species on result in inaccurcies 
  group_by(SSBS) %>% dplyr::mutate(Site_spp = n_distinct(Birdlife_Name),TotalSiteAbundance = sum(SpeciesSiteAbundance)) %>%
  
  ungroup() %>%
  
  dplyr::mutate(RelativeAbundance = SpeciesSiteAbundance/TotalSiteAbundance) %>%
  
  ## how many studies have at least one site of primary minimal for comparisons to be made and make sure the sites have more than a single site. 
  
  droplevels() %>%
  
  data.frame()



trait_range_calc <- function(range, traits){
  
  n_range_2 <- dist(c(min(traits[,2]),max(traits[,2])))[1]
  n_range_3 <- dist(c(min(traits[,3]),max(traits[,3])))[1]
  n_range_4 <- dist(c(min(traits[,4]),max(traits[,4])))[1]
  
  trait_ranges <- list(c(min(traits[,2]) -(range * n_range_2),max(traits[,2]) + (range * n_range_2)),
                       c(min(traits[,3]) -(range * n_range_3),max(traits[,3]) + (range * n_range_3)),
                       c(min(traits[,4]) -(range * n_range_4),max(traits[,4]) + (range * n_range_4)))
  return(trait_ranges)
}


trait_ranges <- trait_range_calc(range = 0.15,traits = TPD_traits$complete_traits)


#### create list of community data



PRED_sites <- rep(list(NA),length(unique(TPD_data$SSBS)))
i <- 1
for(sit in as.character(unique(TPD_data$SSBS))){
  comm_dat <- TPD_data %>% dplyr::filter(SSBS == sit) %>% dplyr::select(Birdlife_Name, RelativeAbundance)
  PRED_sites[[i]] <- comm_dat
  names(PRED_sites)[i] <- sit
  i <- i +1
}



PREDICTS_TPD <- function(site){



comm_sp <- site$Birdlife_Name

mtpd <- FALSE
if(any(comm_sp %in% c(TPD_traits$partial_traits$Birdlife_Name,TPD_traits$single_traits$Birdlife_Name))){
  mtpd <- TRUE
  Mean_sp <- comm_sp[which(comm_sp %in% c(TPD_traits$partial_traits$Birdlife_Name,TPD_traits$single_traits$Birdlife_Name))]

  all_partial <- rbind(TPD_traits$partial_traits,TPD_traits$single_traits)
  
  mean_TPD_dat <- all_partial %>% dplyr::filter(Birdlife_Name %in% Mean_sp) %>% data.frame()
  
  
  mean_TPD <- TPDsMean(species = mean_TPD_dat[,1], means = mean_TPD_dat[,c(2,4,6)], sds = mean_TPD_dat[,c(3,5,7)], trait_ranges = trait_ranges)
  
  }

comm_traits <- TPD_traits[["complete_traits"]] %>% dplyr::filter(Birdlife_Name %in% comm_sp)

trait_density <- TPDs(species = comm_traits[,1], traits = comm_traits[,c(2:4)], trait_ranges = trait_ranges)


if(mtpd){
  trait_density$data$species <- c(trait_density$data$species,mean_TPD$data$species)
  trait_density$TPDs <- c(trait_density$TPDs,mean_TPD$TPDs)
  trait_density$data$traits <- rbind(trait_density$data$traits, mean_TPD$data$means)
  }


comm <- site %>% set_rownames(site$Birdlife_Name) %>% dplyr::select(RelativeAbundance)


Comm_tpd <- TPDc(TPDs = trait_density, sampUnit = t(comm))

Comm_tpd$TPDc <- Comm_tpd$TPDc$TPDc

return(Comm_tpd)

}

require(future)
require(future.apply)

plan(multicore(workers = 8))


PREDICTS_tpds <- future_lapply(PRED_sites,PREDICTS_TPD)



write_rds(file = "Outputs/PREDICTS_sites_tpds.rds", x = PREDICTS_tpds)

######################################################
#####################################################

## Foraging TPDs 

####### get all the standard deviations for the TPDsMean


sds <- sqrt(diag(Hpi.diag(for_traits[["foraging_traits"]][["PCoA_Scores"]][,c(2:4)])))

trait_ranges <- trait_range_calc(range = 0.15, traits = for_traits[["foraging_traits"]][["PCoA_Scores"]])


PREDICTS_TPD_forage <- function(site){
  
  
  
  
  comm_sp <- site$Birdlife_Name
  
  for_TPD_dat <- for_traits[["foraging_traits"]][["PCoA_Scores"]] %>% dplyr::filter(Birdlife_Name %in% comm_sp)
    
  mean_TPD <- TPDsMean(species = for_TPD_dat[,1], means = for_TPD_dat[,c(2:4)], sds = matrix(rep(sds,nrow(site)), ncol = 3, byrow = TRUE),
                       trait_ranges = trait_ranges)
  
  
  comm <- site %>% set_rownames(site$Birdlife_Name) %>% dplyr::select(RelativeAbundance)
  
  
  Comm_tpd <- TPDc(TPDs = mean_TPD, sampUnit = t(comm))
  
  Comm_tpd$TPDc <- Comm_tpd$TPDc$TPDc
  
  return(Comm_tpd)
  
}


For_PREDICTS_tpds <- future_lapply(PRED_sites[1:8],PREDICTS_TPD_forage)



write_rds(For_PREDICTS_tpds, file = "Outputs/PREDICTS_sites_for_tpds.rds")



closeAllConnections()

