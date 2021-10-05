rm(list = ls())

require(tidyverse)
require(TPD)
require(doParallel)
require(foreach)



randomisations <- readRDS("Outputs/site_randomisations.rds")
PREDICTS <- readRDS("Outputs/refined_predicts.rds")
TPD_traits <- readRDS("Outputs/full_morpho_traits_list.rds")
for_traits <- readRDS("Outputs/predicts_foraging_pcoa.rds")
PREDICTS_TPD <- readRDS("Outputs/PREDICTS_sites_tpds.rds")


randomisations <- randomisations[which(names(randomisations) %in% names(PREDICTS_TPD))]

trait_range_calc <- function(range, traits){
  
  n_range_2 <- dist(c(min(traits[,2]),max(traits[,2])))[1]
  n_range_3 <- dist(c(min(traits[,3]),max(traits[,3])))[1]
  n_range_4 <- dist(c(min(traits[,4]),max(traits[,4])))[1]
  
  #######################################
  #####################################
  trait_ranges <- list(c(min(traits[,2]) -(range * n_range_2),max(traits[,2]) + (range * n_range_2)),
                       c(min(traits[,3]) -(range * n_range_3),max(traits[,3]) + (range * n_range_3)),
                       c(min(traits[,4]) -(range * n_range_4),max(traits[,4]) + (range * n_range_4)))
  return(trait_ranges)
}


#################################



# registerDoParallel(cores = 2)
# 
# TPD_randomisations_list <- foreach(list = randomisations[1:2],
#                                                .combine = "c",
#                                                .packages = c("tidyverse", "TPD"),
#                                                .inorder = FALSE) %dopar% {} 


TPD_randomisations_func <- function(list, traits){
  
  site_name <- substr(colnames(list)[3],1,nchar(colnames(list)[3])-9)
  
  comm_sp <- unique(c(as.matrix(list[,-1])))
  
  if(traits == "morpho"){
  
    trait_ranges <- trait_range_calc(range = 0.15, traits = TPD_traits$complete_traits)
    
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
  
  } 
  
  if(traits == "foraging"){
    trait_ranges <- trait_range_calc(for_traits[["foraging_traits"]][["PCoA_Scores"]], range = 0.15)
    
    
    
    comm_traits <- for_traits[["foraging_traits"]][["PCoA_Scores"]] %>% dplyr::filter(Birdlife_Name %in% comm_sp)
  
    trait_density <- TPDsMean(species = comm_traits[,1], means = comm_traits[,c(2:4)], sds = matrix(rep(sds,nrow(comm_traits)), ncol = 3, byrow = TRUE),
                         trait_ranges = trait_ranges)
    
    }
  
  
  
  comm_mat <-  matrix(rep(0,length(comm_sp)*100), nrow = length(comm_sp),ncol = 1000)
  
  rownames(comm_mat) <- comm_sp
  colnames(comm_mat) <- colnames(list)[3:1002]
  
  for(i in 3:ncol(list)){
    species <- list[,i]
    names(species) <- list$RelativeAbundance
    for(j in 1:length(species)){
             comm_mat[species[j],colnames(list)[i]] <- as.numeric(names(species)[j])
      }

  }
  
  
  comm_mat <- as.matrix(rowSums(comm_mat)/1000)
  colnames(comm_mat) <- "Randomisation_Comm"
  
  ran_com_TPD <- TPDc(trait_density,sampUnit = t(comm_mat))
  
  #ran_com_TPD$TPDc$TPDc <- c(ran_com_TPD$TPDc$TPDc,PREDICTS_TPD[[which(names(PREDICTS_TPD) == site_name)]][["TPDc"]])

  #test <- dissim(ran_com_TPD)
  
    
mean_TPDc_mat <- as.numeric(ran_com_TPD$TPDc$TPDc$Randomisation_Comm)

return(mean_TPDc_mat)
}


require(future)
require(future.apply)

plan(multicore(workers = 8))

####morpho list

mean_TPD_randomisations_morpho <- list()

trait_ranges <- trait_range_calc(range = 0.15, traits = TPD_traits$complete_traits)

eval_grid <- TPDs(TPD_traits$complete_traits[c(1:14),1], TPD_traits$complete_traits[c(1:14),c(2:4)], trait_ranges = trait_ranges)

mean_TPD_randomisations_morpho$evaluation_grid <- eval_grid$data$evaluation_grid


TPD_randomisations_list <- future_lapply(randomisations[1:8],TPD_randomisations_func, traits = "morpho")

TPD_randomisations_morpho <- c(mean_TPD_randomisations_morpho,TPD_randomisations_list)

##### foraging list

rm(TPD_randomisations_list)

mean_TPD_randomisations_foraging <- list()

trait_ranges <- trait_range_calc(range = 0.15, traits = for_traits[["foraging_traits"]][["PCoA_Scores"]])
sds <- sqrt(diag(Hpi.diag(for_traits[["foraging_traits"]][["PCoA_Scores"]][,c(2:4)])))


eval_grid <- TPDsMean(species = for_traits[["foraging_traits"]][["PCoA_Scores"]][c(1:14),1], 
                      means = for_traits[["foraging_traits"]][["PCoA_Scores"]][c(1:14),c(2:4)], 
                      sds = matrix(rep(sds,14), ncol = 3, byrow = TRUE),
                      trait_ranges = trait_ranges)

mean_TPD_randomisations_foraging$evaluation_grid <- eval_grid$data$evaluation_grid

TPD_randomisations_list <- future_lapply(randomisations[1:8],TPD_randomisations_func, traits = "foraging")

TPD_randomisations_foraging <- c(mean_TPD_randomisations_foraging,TPD_randomisations_list)

write_rds(file = "Outputs/randomisations_TPD_for.rds", TPD_randomisations_foraging)
write_rds(file = "Outputs/randomisations_TPD_morpho.rds", TPD_randomisations_morpho)
  
closeAllConnections()

