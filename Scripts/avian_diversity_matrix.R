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

## list all file paths that we have species range data for

species_files <- paste("../../Datasets/Birdlife_Maps/Shapefiles/PREDICTS_BL/",
                       list.files("../../Datasets/Birdlife_Maps/Shapefiles/PREDICTS_BL/"),
                       sep = "")

#### get a blank world map at a one degree resolution -- approx 100km2 at the equator 

blank <- raster("../../Datasets/Environmental_Variables/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2000_1_deg.tif")
blank <- reclassify(blank, c(-Inf,Inf,0))

#### write blank raster map for future use

write_rds(x = blank, file = "Outputs/blank_map_1deg.rds")

### This loop is going to rasterize all species ranges create presence absence maps and generally a global map of avian diversity

registerDoParallel(cores = 8)

map_data <- foreach(range = species_files,
                    .combine = "rbind",
                    .packages = c("tidyverse", "maptools", "rgdal", "sp", "raster","rgeos","sf","fasterize"),
                    .inorder = FALSE) %dopar% {
                      
                      
                      load(range)
                      
                      #### filter data so that distirbution polygons are those that reprsent extant, probably extnat and possibly extant presence,
                      #### native, reintroduced and introduced origin and those that are resident, breeding or non-breeding rnages                        
                      
                      
                      data <- data %>% dplyr::filter(any(PRESENCE == c(1,2,3)), any(ORIGIN == c(1,2,3)), any(SEASONAL == c(1,2,3)))                      
                      
                      
                      ##### funny thing with some ranges so need to convert back to type MULTIPOLYGON
                      
                      if(any(class(data$Shape)[1] == "sfc_MULTISURFACE", class(data$Shape)[1] == "sfc_GEOMETRY")){
                        for(k in 1:NROW(data)){
                          data$Shape[[k]] <- st_cast(data$Shape[[k]], "MULTIPOLYGON")
                        }
                      }
                      
                      #combine all polygons together and convert to a sf class polygon
                      
                      shape <- as_Spatial(st_combine(data$Shape))
                      shape <- st_as_sf(shape)
                      
                      ### convert polygon into a presence raster on the blank world map 
                      
                      ras <- fasterize(sf = shape, raster = blank,fun = "count")
                      
                      #### extract cells which overlap with the polygon and exclude marine and oceanic ranges by adding blank map 
                      #### that has non-terrestrial cells as NA
                      
                      ras_data <- (blank + ras)
                      
                  plot(shape)
                      #### create matrix where the row is the species and each column is a cell in the raster - presence indicated by a 1
                      
                      mat <- matrix(ras_data, nrow = 1, ncol = length(ras_data))
                      colnames(mat) <- 1:length(ras_data)
                      rownames(mat) <- data$SCINAME[1]
                      
                      return(mat)
                      
                    }

registerDoSEQ()
closeAllConnections()

## save output map data that can then be put given to the blank raster

write_rds(x= map_data, file = "Outputs/raster_species_matrix.rds")

