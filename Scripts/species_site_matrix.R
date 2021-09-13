rm(list = ls())

library(tidyverse)
library(maptools)
library(rgdal)
library(sp)
library(raster)
library(rgeos)
require(sf)
require(doParallel)
require(foreach)

PREDICTS <- readRDS("Outputs/refined_predicts.rds")

sites <- PREDICTS %>% distinct(SSBS,Latitude,Longitude)


species_files <- paste("../../Datasets/Birdlife_Maps/Shapefiles/PREDICTS_BL/",
                       list.files("../../Datasets/Birdlife_Maps/Shapefiles/PREDICTS_BL/"),
                                  sep = "")




site_coords <- SpatialPoints(sites[,c("Longitude","Latitude")], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))




registerDoParallel(cores = 8)

site_spp_matrix <- foreach(range = species_files,
                           .packages = c("tidyverse","maptools","rgdal","sp","raster","rgeos","sf"),
                           .combine = "cbind",
                           .inorder = FALSE) %dopar% {

  load(range)
  
  
  if(any(class(data$Shape)[1] == "sfc_MULTISURFACE", class(data$Shape)[1] == "sfc_GEOMETRY")){
    for(k in 1:NROW(data)){
      data$Shape[[k]] <- st_cast(data$Shape[[k]], "MULTIPOLYGON")
    }
  }
  #combine all polygons together
  
  shape <- as_Spatial(st_combine(data$Shape))
  
mat <- matrix(rep(NA,nrow(sites)),nrow = nrow(sites), ncol = 1)
  colnames(mat) <- as.character(data$SCINAME)[1]
  rownames(mat) <- as.character(sites$SSBS)

  for(j in 1:nrow(sites)){

    Overlap <- suppressWarnings(gDistance(site_coords[j],shape))
  
    if(Overlap == 0){
      mat[j,1] = 1
    }  
    
  }
  
  return(mat)
  
}


registerDoSEQ()
closeAllConnections()


write_rds("Outputs/site_spp_matrix.rds",x =  site_spp_matrix)


sum(site_spp_matrix[1,], na.rm = TRUE)
rowSums(site_spp_matrix, na.rm = TRUE)
sum_site <- rowSums(site_spp_matrix, na.rm = TRUE)


