rm(list = ls())

library(tidyverse) ## for wrangling and general data handling
library(terra) # mapping
library(stars)
library(sp) ## mapping 
library(raster) ## mapping
require(sf) ## mapping 
require(doParallel) ## parallelisation
require(foreach) ## parallelisation
require(fasterize)


species_files <- list.files("../../Datasets/Birdlife_Maps/Shapefiles/PREDICTS_BL/", full.names = TRUE)


registerDoParallel(cores = 8)

BL_sp_polys <- foreach(range = species_files,
                    .combine = "rbind",
                    .packages = c("tidyverse", "sp", "raster","sf","fasterize"),
                    .inorder = FALSE) %dopar% {
                      
                      data <- readRDS(range)
                      
                      
                      sp_name <- data$binomial[1]
                      
                      #### filter data so that distirbution polygons are those that reprsent extant, probably extnat and possibly extant presence,
                      #### native, reintroduced and introduced origin and those that are resident, breeding or non-breeding rnages
                      
                      
                      data <- data %>% dplyr::filter(presence %in% c(1,2,3), origin %in% c(1,2,3), seasonal %in% c(1,2,3))
                      
                      
                      if(nrow(data) == 0){
                      } else {
                        
                        
                        
                        ##### funny thing with some ranges so need to convert back to type MULTIPOLYGON
                        
                        
                        if(any(class(data$Shape)[1] == "sfc_MULTISURFACE", class(data$Shape)[1] == "sfc_GEOMETRY")){
                          for(k in 1:NROW(data)){
                            data$Shape[[k]] <- st_cast(data$Shape[[k]], "MULTIPOLYGON")
                          }
                        }
                        
                        
                        sf_use_s2(TRUE)
                        if(any(!st_is_valid(data$Shape))){
                          sf_use_s2(FALSE)
                        }
                        
                        t <- try(shape <- st_union(data$Shape))
                        
                        if(inherits(t, "try-error")){
                          shape <- st_combine(data$Shape)
                        }
                        
                        
                        st_crs(shape) <- "+proj=longlat +datum=WGS84 +no_defs"
                        
                        shape <- st_as_sf(shape)
                        
                        shape$Species <- sp_name
                        
                        ### convert polygon into a presence raster on the blank world map
                       return(shape)
                      }
} 



registerDoSEQ()
closeAllConnections()

write_rds(file = "Outputs/BL_species_polygons.rds", BL_sp_polys)