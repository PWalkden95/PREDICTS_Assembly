rm(list = ls())

require(tidyverse)
require(sf)
require(sp)
require(terra)
require(doParallel)


island_polys <- readRDS("Outputs/assembly_islands.rds")

species_maps <- list.files("../../Datasets/Birdlife_Maps/Shapefiles/PREDICTS_BL/", full.names = TRUE)




island_species <- function(island,sp_maps, n_cores){

  isle_poly <- st_as_sf(island)
  
  
  registerDoParallel(cores = n_cores)
  
  island_sp <- foreach(range = sp_maps, 
                       .combine = "c",
                       .packages = c("tidyverse","sf", "sp", "terra")) %dopar% {
  
  
  data <- readRDS(range)
  
  data <- data %>% dplyr::filter(presence %in% c(1,2,3), origin %in% c(1,2,3), seasonal %in% c(1,2,3))
  
  if(nrow(data) == 0){
    sp <- NA
  } else {
  
  
  if(any(class(data$Shape)[1] == "sfc_MULTISURFACE", class(data$Shape)[1] == "sfc_GEOMETRY")){
    for(k in 1:NROW(data)){
      data$Shape[[k]] <- st_cast(data$Shape[[k]], "MULTIPOLYGON")
    }
  }
  
  for(i in 1:nrow(data)){
    data$Shape[i] <- st_make_valid(data$Shape[i])
  }
  

  sf_use_s2(TRUE)
  if(any(!st_is_valid(data$Shape))){
    sf_use_s2(FALSE)
  }
  
  shape <- st_make_valid(st_union(data$Shape))
  st_crs(shape) <- st_crs(isle_poly)
  
  
  if(!st_is_valid(shape)){
    sf_use_s2(FALSE)
  }
  
  
  
  sp <- NA
  if(st_intersects(shape, y = isle_poly, sparse = FALSE)[1,1]){
    sp <- data$SCINAME[1]
  }
  }
  
  return(sp)
                       }

  closeAllConnections()
  registerDoSEQ()
  
  island_sp <- na.omit(island_sp)
  
  return(island_sp)
  
}

island_spp <- lapply(island_polys, island_species, sp_maps = species_maps, n_cores = 8)


write_rds(island_spp, file = "Outputs/assembly_island_spp.rds")
