rm(list = ls())

#### Load in packages 


library(tidyverse)
library(maptools)
library(rgdal)
library(sp)
library(raster)
library(rgeos)
require(sf)
require(doParallel)
require(foreach)
require(fasterize)



###############################
##### ADF species #############
###############################

adf <- readRDS("Outputs/predicts_site_cell_matrix.rds")
bird_diversity <- readRDS("Outputs/raster_species_matrix.rds")
blank <- readRDS("Outputs/blank_map_1deg.rds")

plot(blank)
for(i in 1:50){
blank@data@values <- adf[i,]

plot(blank, main = paste(i))

}
test <- map_data[,c(adf_cells)]

spp <- data.frame(Birdlife_Name = names(which(rowSums(test, na.rm = TRUE) > 0))) %>% 
  dplyr::left_join(AVONET[,c("Birdlife_Name", "Birdlife_family")], by = "Birdlife_Name") %>%
  dplyr::filter(!(Birdlife_family %in% drop))



###############################################
####### SPECIES TO DROP #######################
###############################################


### Swallows, Martins, Swifts, Accipiter, Vultures, Waterbirds, freshwater birds, nocturnak birds etc should all be dropped

AVONET <- read.csv("../../Datasets/GBD/GBD_2021_BirdLife_Taxo_19 April.csv")


Nocturnal_spp <- AVONET %>% dplyr::filter(Nocturnal == 1) %>% pull(Birdlife_Name)
marine_spp <- AVONET %>% filter(Primary_Habitat_updated == "Marine") %>% pull(Birdlife_Name)
swift_spp <- AVONET %>% filter(Birdlife_family %in%  c("Apodidae","Hemiprocnidae")) %>% pull(Birdlife_Name)
swallow_spp <- AVONET %>% filter(Birdlife_family %in% c("Hirundinidae","Artamidae")) %>% pull(Birdlife_Name)
raptor_spp <- AVONET %>% filter(Birdlife_family %in% c("Accipitridae","Cathartidae")) %>% pull(Birdlife_Name)
coastal_spp <- AVONET %>% filter(Primary_Habitat_updated == "Coastal") %>% pull(Birdlife_Name)

drop_spp <- unique(c(Nocturnal_spp,marine_spp,swift_spp,swallow_spp,raptor_spp,coastal_spp))



