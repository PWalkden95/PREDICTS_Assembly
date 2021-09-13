rm(list = ls())

require(tidyverse)
require(sf)
require(sp)
require(rgeos)


island_polys <- readRDS("Outputs/assembly_islands.rds")

species_maps <- list.files("../../Datasets/Birdlife_Maps/Shapefiles/PREDICTS_BL/", full.names = TRUE)


island_spp <- rep(list(NA),4)
names(island_spp) <- names(island_polys)

j <- 1

for(range in species_maps){
  
  print(j)
  
  load(range)
  
  data <- data %>% dplyr::filter(presence %in% c(1,2,3), origin %in% c(1,2,3), seasonal %in% c(1,2,3))
  
  if(nrow(data) == 0){
    next()
  }
  
  if(any(class(data$Shape)[1] == "sfc_MULTISURFACE", class(data$Shape)[1] == "sfc_GEOMETRY")){
    for(k in 1:NROW(data)){
      data$Shape[[k]] <- st_cast(data$Shape[[k]], "MULTIPOLYGON")
    }
  }
  
  shape <- as_Spatial(st_combine(data$Shape))
  
  for(i in 1:length(island_polys)){
  distance <- suppressWarnings(gDistance(shape,island_polys[[i]]))
  
  if(distance == 0){
    island_spp[[i]] <- na.omit(c(island_spp[[i]],data$SCINAME[1]))
  }
  
  }
  
  j <- j +1
}



write_rds(island_spp, file = "Outputs/assembly_island_spp.rds")
