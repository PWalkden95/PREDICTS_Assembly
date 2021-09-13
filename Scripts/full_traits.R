####################################################################################################
# This script generates the trait axes used in the trait probability density functions, which utilises AVONET to it's fullest extent as it accounts for
# intraspecific trait variation to calculate a probabilistic hypervolume in trait space. This entails performing the two step PCA again but with all
# traits considered for the generation of the morphological trait space. For the generation of the axes for Dietary trait space, unfortunately there,
# is no data available for the differences between individuals diet at the species level so foraging scores are used and variation within species is 
# estimated using plug-in bandwidth estimators in the ks package. 

# The latter part of the script is running some checks with the TPD package as occasionally the output will generate an out consisting entirely of 
# NaNs and by altering the traits ranges and number of divisions can resolve the problem - However, when I tried to resolve all species it wasn't possible
# so to save computational power and for simplicity those species that are unable to create TPDs at the default (50 divisions) they are passed onto the
# TPDsMean function with the species that have insufficient specimens and the bandwidth of the now multivariate normal distirbution is the determinded in
# a couple of ways i) if there are greater than specimen the bandwidth is the sd ii) otherwise bandwidth is estimated using the plug-in bandwidth estimator
# based on all the specimens

## The outputs of this script are two lists the first that contains three data frames with the necessary scores to calculate the TPDs for the three 
## different cases and the second is a dataframe that contains the PCoA scores for the dietary guild traits.



rm(list = ls())



require(tidyverse) # data manipulation
require(magrittr) # pipiing
require(vegan) ## ordination
require(TPD) ## TPD checking
require(ks) # bandwidth estimation
require(ggbiplot)

## load in databases and species pools

AVONET_full <- read.csv("../../Datasets/GBD/GBD_BiometricsRaw_combined_11_Aug_2021_MASTER.csv")
species_pools <- read_rds("Outputs/predicts_sites_species_pools.rds")

## extract all species from the species pool list 

assembly_species <- unique(c(unlist(unlist(species_pools)))) 


###################################################
#### MORPHOLOGICAL TRAIT AXES #####################
###################################################


### Join all the morpho traits 

assem_traits <- AVONET_full %>% dplyr::filter(Species1_BirdLife %in% assembly_species) %>%
  
  dplyr::select(Species1_BirdLife,Beak.Length_Culmen,Beak.Length_Nares,Beak.Width,Beak.Depth,Tarsus.Length,Wing.Length,
                Kipps.Distance,Secondary1,Hand.wing.Index,Tail.Length)  %>% data.frame()

drop_rows <- c()
for(i in 1:nrow(assem_traits)){
  if(all(is.na(assem_traits[i,c(2:10)]))){
  drop_rows <- c(drop_rows,i)  
  }
}

assem_traits <- assem_traits %>% slice(-drop_rows)

##################################################
# some species traits contain NAs so we don't want to ditch the rest of the data that could be useful so the NAs are going to be imputed as the mean
# of the rest of the specimens


for(trait in colnames(assem_traits)[2:11]){
  
  data <- assem_traits %>% dplyr::select(Species1_BirdLife,all_of(trait))
  sp <- unique(data[which(is.na(data[,trait])),"Species1_BirdLife"])

  
  for(s in sp){
    mean <- mean(assem_traits[which(assem_traits[,"Species1_BirdLife"] == s),trait], na.rm = TRUE)
    assem_traits[which(assem_traits[,"Species1_BirdLife"] == s & is.na(assem_traits[,trait])),trait] <- mean
  }
  
}

assem_traits <- assem_traits %>% na.omit() %>% 
  
  group_by(Species1_BirdLife) %>%
  
  dplyr::mutate(num = n()) %>% ungroup()

##### 7265 species that have enough traits in AVONET (4 or greater specimens)


enough_sp <- assem_traits %>% filter(num > 3)

e_sp <- unique(enough_sp$Species1_BirdLife)

length(e_sp)

######### 199  sp 

insuf_sp <- assem_traits %>% filter(num < 4, num > 1)

in_sp <- unique(insuf_sp$Species1_BirdLife)

length(in_sp)



########################### 88 sp 

sing_sp <- assem_traits %>% filter(!(Species1_BirdLife %in% c(in_sp,e_sp)))

s_sp <- unique(sing_sp$Species1_BirdLife)

length(s_sp)

#################
################ 252 no specimens

no_sp <- assembly_species[!(assembly_species %in% c(e_sp,in_sp,s_sp))]


mean_no_sp <- data.frame(Birdlife_Name = c(s_sp,no_sp)) 


mean_traits <- readRDS("../PREDICTS_Taxonomy/PREDICTS_imputed_BL_traits.rds") 
pseudospp <- mean_traits %>% dplyr::filter(grepl(uniqueIDS, pattern = "_")) %>% dplyr::rename(Birdlife_Name = uniqueIDS)


mean_no_sp <- mean_no_sp %>% dplyr::left_join(mean_traits, by = c("Birdlife_Name" = "uniqueIDS")) %>% rbind(pseudospp) %>% dplyr::mutate(num = 1) 

############################################################################################################################################
#############################################################################################################################################



colnames(enough_sp)[1] <- "Birdlife_Name"
colnames(insuf_sp)[1] <- "Birdlife_Name"
colnames(mean_no_sp)[10] <- colnames(enough_sp)[10]

PCA_traits <- rbind(enough_sp[,c(1:11)],insuf_sp[,c(1:11)],mean_no_sp[,c(1:11)]) %>% na.omit()



PCA_traits[,c(2:11)] <- scale(log(PCA_traits[,c(2:11)]))

For_Data <- data.frame(PCA_traits[,c("Birdlife_Name", "Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Width", "Beak.Depth")])

for_pca <- prcomp(For_Data[,-1],center = TRUE,scale. = TRUE)

for_pca
summary(for_pca)

ggbiplot(for_pca)
#######################################

Loco_Data <- data.frame(PCA_traits[,c("Birdlife_Name", "Tarsus.Length","Wing.Length","Secondary1","Tail.Length")])

loco_pca <- prcomp(Loco_Data[,-1], center = TRUE,scale. = TRUE)

loco_pca
summary(loco_pca)

ggbiplot(loco_pca)

###############################

Body_Data <- data.frame(forage = for_pca$x[,1],loco = loco_pca$x[,1])

body_pca <- prcomp(Body_Data,center = TRUE,scale. = TRUE)

body_pca
summary(body_pca)

ggbiplot(body_pca)

###################################################

Full_PCA_Score <- data.frame(Birdlife_Name = PCA_traits$Birdlife_Name, locomotion = loco_pca$x[,2],
                             foraging = for_pca$x[,2],
                             body = body_pca$x[,1])


###########################################

trait_range_calc <- function(range = 0.15, traits){
  
  n_range_2 <- dist(c(min(traits[,2]),max(traits[,2])))[1]
  n_range_3 <- dist(c(min(traits[,3]),max(traits[,3])))[1]
  n_range_4 <- dist(c(min(traits[,4]),max(traits[,4])))[1]
  
  trait_ranges <- list(c(min(traits[,2]) -(range * n_range_2),max(traits[,2]) + (range * n_range_2)),
                       c(min(traits[,3]) -(range * n_range_3),max(traits[,3]) + (range * n_range_3)),
                       c(min(traits[,4]) -(range * n_range_4),max(traits[,4]) + (range * n_range_4)))
  return(trait_ranges)
}


###################### Function to find out which species a TPD couldn't be calculated for with the range and division combination given

TPD_fucked <- function(TPD){

  Trait_density <- TPD
  
  fucked_species <- c()
  for(i in 1:length(Trait_density$TPDs)){
    if(any(is.nan(Trait_density$TPDs[[i]]))){
      fucked_species <- c(fucked_species,names(Trait_density$TPDs)[[i]]) 
    }
  }
 return(fucked_species)
}

##########################################
### function that iterates through different range and division combinations to find a solution that satisfies the TPDs for the problematic species


TPD_opt <- function(sp_traits, start_values = c(range = 0.01, division = 50), by = 0.01){
  
  range <- start_values["range"]
  division <- start_values["division"]
  

  
  
  fucked_species_traits <- sp_traits
  
  any_nan <- TRUE

  
  
  
  
  while(any_nan){
    print(paste("trait range equals", as.character(range), "number of divisions are", as.character(division)))
    T_d <- TPDs(fucked_species_traits[,1], fucked_species_traits[,c(2:4)], trait_ranges = trait_range_calc(range), n_divisions = division, alpha = 0.95)
    
    any_nan <- FALSE
    for(k in 1:length(T_d$TPDs)){
      if(any(is.nan(T_d$TPDs[[k]]))){
        any_nan <- TRUE
        print(paste(names(T_d$TPDs)[[k]]))
      }
    }
    
    
    if(any_nan){
      range <- round(c(range  + by),5)
      trait_ranges <- trait_range_calc(range)
    
    }
    
    if(range > 0.50){
      range <- start_values["range"]
      division <- division +1
    }
    
    
  }
  
  
  return(c(range = as.numeric(range), division = as.numeric(division)))
  
}

#################################################################
################################################################
### Loop to find a TPD combination that satisfies the solution for all species 

 
####### This didnt't really work couldn;t whittle it down to sp thay every species was able to have a feasible TPD so going to run it on default
### and those species I will then calculate their TPD using tPDsmean using a multivariate normal distribution based on teh trait SDs and mean values.  

####################### so heres some redundant code....

# 0.06 /52 2 species

# while(!is.null(f_sp)){
#   
#   rm(Trait_density)
#   
#   problem_species <- unique(c(problem_species,f_sp))
#   
#   fucked_traits <- e_pca %>% dplyr::filter(Birdlife_Name %in% problem_species)
#   
#   opt_values <- TPD_opt(fucked_traits, start_values = opt_values, by = 0.01)
#   
#   Trait_density <- TPDs(e_pca[,1],e_pca[,c(2:4)], n_division = opt_values["division"], trait_ranges = trait_range_calc(opt_values["range"]),alpha = 0.95)
#   
#   
#   f_sp <- TPD_fucked(Trait_density)
#   
# }

################################################################
################################################################
################################################################

## get the traits scores for all the species that a TPD for the species can be calculated 

e_pca <- Full_PCA_Score %>% dplyr::filter(Birdlife_Name %in% e_sp) %>% na.omit()

## run an initial TPD to flag some potentially problematic species

trait_ranges <- trait_range_calc(range = 0.15, traits = e_pca)

Trait_density <- TPDs(e_pca[,1], e_pca[,c(2:4)], trait_ranges = trait_ranges, alpha = 0.95)


f_sp <- TPD_fucked(TPD = Trait_density)



e_pca <- e_pca %>% dplyr::filter(!(Birdlife_Name %in% f_sp))





#########################################
###########################################
## collate all the traits together to be in a workable format for the computation of TPDs

## first is the full traits that can be used to calculate TPD

f_traits <- list()
f_traits$complete_traits <- e_pca


### for the other species their occupancy for functional trait space is determined by a multivariate normal distribution centered on the traits mean 
### and the "bandwidth" can be decided in a couple of ways depending on how many specimens that trait values are available for the species

# 1) if there are more than one specimens than the bandwidth is the SD of the available trait values

p_traits <- Full_PCA_Score %>% dplyr::filter(Birdlife_Name %in% c(f_sp,in_sp)) %>% dplyr::group_by(Birdlife_Name) %>% 
  dplyr::summarise(loco_mean = mean(locomotion), loco_sd = sd(locomotion),
                   for_mean = mean(foraging), for_sd = sd(foraging),
                   body_mean = mean(body), body_sd = sd(body)) %>% data.frame()


### however also need to check here whether any of the mean calculations with the informed SDs throw up a problem as well so....

mean_trait_density <- TPDsMean(species = p_traits[,1], 
                               means = p_traits[,c(2,4,6)],
                               sds = p_traits[,c(3,5,7)], trait_ranges = trait_ranges, alpha = 0.95)

m_f_sp <- TPD_fucked(mean_trait_density)

p_traits <- p_traits %>% dplyr::filter(!(Birdlife_Name %in% m_f_sp))

## 2) If there is only a single specimen then the standard deviation used is 0.5 * the global SD.

mean_trait_sds <- sqrt(diag(Hpi.diag(Full_PCA_Score[,c(2:4)])))

s_traits <- Full_PCA_Score %>% dplyr::filter(Birdlife_Name %in% c(mean_no_sp$Birdlife_Name, m_f_sp)) %>% dplyr::group_by(Birdlife_Name) %>% 
  dplyr::summarise(loco_mean = mean(locomotion), loco_sd = mean_trait_sds[1],
                   for_mean = mean(foraging), for_sd =  mean_trait_sds[2],
                   body_mean = mean(body), body_sd = mean_trait_sds[3]) %>% data.frame()



###################################
################# join these traits to the list 

f_traits$partial_traits <- p_traits

f_traits$single_traits <- s_traits

######## write the list for later use 

write_rds(file = "Outputs/full_morpho_traits_list.rds", x = f_traits)

###############################################################################################
##### FORAGING TRAITS FOR FULL MORPHO + TPDs #####################
############################

source("Functions/trait_scores.R")


Forage <- readRDS("../PREDICTS_Taxonomy/PREDICTS_imputed_BL_traits_forage.rds")

For_assembly_traits <- data.frame(Birdlife_Name = c(assembly_species,pseudospp$Birdlife_Name)) %>% dplyr::left_join(Forage)

predicts_for_traits <- trait_scores(data = For_assembly_traits, traits = "foraging")

write_rds(predicts_for_traits, file = "Outputs/predicts_foraging_pcoa.rds")

