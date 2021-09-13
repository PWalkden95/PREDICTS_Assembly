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




########################################################
###### PREDICTS SITES ASSEMBLAGE SIMILARITY MATRIX #####
########################################################

### Load in refined PREDICTS and raster species presence data 

PREDICTS <- readRDS("Outputs/refined_predicts.rds")
map_data <- readRDS("Outputs/raster_species_matrix.rds")
blank <- readRDS("Outputs/blank_map_1deg.rds")

#### extract site coordinates for which we want to calculate an assemblage dispersion field for

sites <- PREDICTS %>% distinct(SS,SSBS,Latitude,Longitude, .keep_all = FALSE) %>% data.frame()


#### Which cell is each site located in 

sites$cell <- raster::cellFromXY(blank, sites[,c(1:2)])


## This loop is going to create a matrix with rows being the cell in which a site is located and columns being all other raster cells
## the value of each cell will the proportion of species shared between each cell

registerDoParallel(cores = 8)

site_site_matrix <- foreach(cell = unique(sites$cell),
                            .combine = "rbind",
                            .packages = c("tidyverse", "raster", "maptools", "rgeos", "sf", "sp"),
                            .inorder = FALSE) %dopar% {
                              
                              ### This bit is kinda redundant when looking just at sites but if we were to do this for all cells then this would create an row of NAs when 
                              ### the target cell had no species present within - this includes marine and oceanic cells.                               
                              
                              if(!any(!is.na(as.numeric(map_data[,cell])))){
                                mat <- matrix(rep(NA,ncol(map_data)), nrow = 1, ncol = ncol(map_data))
                                rownames(mat) <- cell
                              } else {
                                
                                
                                ## extract species that are present within a specified cell 
                                
                                focal_spp <- names(which(map_data[,cell] == 1))
                                
                                ### get all other cells minus the focal cell
                                
                                other_cells <- c(1:(cell-1),(cell+1):ncol(map_data))
                                
                                ### create blank matrix 
                                
                                mat <- matrix(rep(NA,ncol(map_data)), nrow = 1, ncol = ncol(map_data))
                                
                                ## cell shares all species with itself
                                
                                mat[,cell] <- 1
                                rownames(mat) <- cell
                                
                                ### for all other cells extract species and calculate proportion of pseices shared with focal cell 
                                
                                for(other in other_cells){
                                  
                                  if(!any(!is.na(as.numeric(map_data[,other])))){
                                    prop <- NA
                                  } else {
                                    
                                    other_spp <- names(which(map_data[,other] == 1))
                                    
                                    if(is_empty(other_spp)){
                                      other_spp <- "none"
                                    }
                                    
                                    sim <- length(which(focal_spp %in% other_spp)) 
                                    
                                    if(sim == 0 ){
                                      prop <- 0
                                    } else {
                                      prop <- sim/length(focal_spp)
                                    }
                                    
                                  }
                                  
                                  mat[,other] <- prop
                                }
                                
                              }
                              
                              return(mat)
                            }

registerDoSEQ()
closeAllConnections()


write_rds(file = "Outputs/predicts_site_cell_matrix.rds", x = site_site_matrix)