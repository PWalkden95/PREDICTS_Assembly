require(tidyverse)
require(magrittr)
require(SYNCSA)



Rao_Q_Func_bias <- function(data, traits){
  
  
  ### get the list of uncique species across teh whole dataset
  if(class(traits) == "dist"){
    Species_abundance <- data.frame(Birdlife_Name = attr(traits,which = "Labels"))
  } else {
    Species_abundance <- data.frame(Birdlife_Name = traits[,"Birdlife_Name"])
  }
  ### loop over every site in the dataset to collate the relative abundance of each species
  
  sites <- as.character(unique(data$SSBS))
  
  for(site in sites){
    
    Spp_abd <- data %>% filter(SSBS == site) %>% droplevels() %>% data.frame()
    
    
    ### Join species withining site to dataframe and rename column as site name
    
    Species_abundance <- Species_abundance %>% left_join(Spp_abd[,c("Birdlife_Name", "RelativeAbundance")], by = "Birdlife_Name")
    colnames(Species_abundance)[which(colnames(Species_abundance) == "RelativeAbundance")] <- paste(site)
  }
  
  ## rename rows as species and drop column from dataset 
  
  rownames(Species_abundance) <- Species_abundance$Birdlife_Name
  Species_abundance <- as.matrix(Species_abundance[,-1])
  
  ## Nas to zeros
  
  Species_abundance[is.na(Species_abundance)] <- 0  
  
  ### Join all species in datasets traits scores
  
  if(class(traits) == "dist"){
    spp_traits <- traits  
  } else {
    
    spp_traits <- traits %>% filter(Birdlife_Name %in% rownames(Species_abundance))
    rownames(spp_traits) <- spp_traits$Birdlife_Name
    spp_traits <- spp_traits[,-1] %>% as.matrix()
  }
  
  source("C:/Users/patri/OneDrive - Imperial College London/Work/Biology/(2020-_Natural_History_Museum/PhD_Code/PREDICTS_Assembly/Functions/rao_diversity_gaw.R")
  
  Rao_Bias <- rao_diversity_gaw(comm = t(Species_abundance),traits =  spp_traits)#THIS is using the package SYNCSA that calcuates Rao's using gowdis 
  
  return(Rao_Bias)
}
######################################################################
########## Here we are also going to calculate an "unbiased" Raos Q###
######################################################################

Rao_Q_Func_unbias <- function(data,traits){  
  
  
  if(class(dist) == "dist"){
    Species_abundance_2 <- data.frame(Birdlife_Name = attr(traits, which = "Labels"))
  } else{
    Species_abundance_2 <- data.frame(Birdlife_Name = traits[,"Birdlife_Name"])
  }
  ### loop over every site in the dataset to collate the relative abundance of each species
  
  for(site in levels(data$SSBS)){
    
    Spp_abd <- data %>% filter(SSBS == site) %>% droplevels() %>% data.frame()
    
    
    ### Join species withining site to dataframe and rename column as site name
    
    Species_abundance_2 <- Species_abundance_2 %>% left_join(Spp_abd[,c("Birdlife_Name", "SpeciesSiteAbundance")], by = "Birdlife_Name")
    colnames(Species_abundance_2)[which(colnames(Species_abundance_2) == "SpeciesSiteAbundance")] <- paste(site)
  }
  
  ## rename rows as species and drop column from dataset 
  
  rownames(Species_abundance_2) <- Species_abundance_2$Birdlife_Name
  Species_abundance_2 <- as.matrix(Species_abundance_2[,-1])
  
  ## Nas to zeros
  
  Species_abundance_2[is.na(Species_abundance_2)] <- 0  
  
  
  comm <- t(Species_abundance_2)
  
  #### species traits 
  if(class(traits == "dist")){
    spp_traits <- traits 
  } else {
    
    spp_traits <- traits %>% filter(Birdlife_Name %in% rownames(Species_abundance_2))
    rownames(spp_traits) <- spp_traits$Birdlife_Name
    spp_traits <- spp_traits[,-1]
  }
  
  ##Load in altered SYNCSA function to calculate Rao's Q unbias 
  source("../Functional_Intactness_Index/Functions/Rao_Diversity_2.R")
  
  
  Rao_Unbias <- rao_diversity_2(comm = comm, traits = spp_traits)
  
  
  return(Rao_Unbias)
}
