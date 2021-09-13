require(tidyverse)
require(gtools)
require(geosphere)
require(sp)
require(geodist)
require(adehabitatHR)
require(rgeos)

#### This is a function that will combine PREDICTS studies when there are multiple that survey the same site - typically this is when the study has conducted
#### sampling over multiple years or in different seasons, a small amount of prio inspection shouuld be required to detremine whether this is a necessary
### step

coordinate_resolution <- function(data){
  
  co_dat <- data
  
  for(site in as.character(unique(co_dat$SSBS))) {
    
    ### filter just coords
    coords <- co_dat %>% dplyr::filter(SSBS == site) %>% distinct(Longitude,Latitude) %>% data.frame()
    
    ### if just one in site no need to change
    if(nrow(coords) > 1){
      
      ### however if more needed calculate centroid
      cat(sprintf("Site %s has %i different coordinates calculating centroid... \n", site, nrow(coords)))
      
      #### if less than four just calculate mean as you cannot form a polygon with less
      if(nrow(coords < 4)){
        Long <- mean(coords[,"Longitude"])  
        Lat <- mean(coords[,"Latitude"])
        
      } else {
        
        ### if four calculate the centroid of the forming polygon
        if(nrow(coords == 4)){
          Long <- centroid(coords)[1]
          Lat <- centroid(coords)[2]
        } else {
          
          ### if there are greater than four points we can also identify points that are likely to be outliers so to do this we can calculate a 
          ### minimum convex polygon around the points and remove outliers that are present 
          
          spat_coords <- SpatialPoints(SpatialCoord<-SpatialPointsDataFrame(coords = coords[,c("Longitude","Latitude")], data = coords,   ##Convert regular data frame into
                                                                            proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")))
          mcp <- mcp(spat_coords, percent = 95)
          
          
          drop_coords <- c()
          for(k in 1:nrow(coords)){
            if(rgeos::gDistance(mcp,spat_coords[k]) > 0){
              drop_coords <- c(drop_coords,k)
            }          }
          
          if(!is.null(drop_coords)){
            cat(sprintf("%i outlier coordinates have been removed",length(drop_coords)))
            coords <- coords[-drop_coords]
          }
          
          Long <- centroid(coords)[1]
          Lat <- centroid(coords)[2]
          
        }
      }
      coords_res <- co_dat %>% dplyr::filter(SSBS == site)
      coords_res$Longitude <- Long
      coords_res$Latitude <- Lat
      
      co_dat <- co_dat %>% dplyr::filter(SSBS != site) %>% rbind(coords_res)
    }
  }
  
  return(co_dat)
  
}




combine_study <- function(data, stud_combine) {
  
  ### for the SOURCE ID that contains multiple studies 
  
  study <- data %>% filter(grepl(SS, pattern = stud_combine)) %>% droplevels()
  
  cat(sprintf("Number of studies in source: %i\n", length(unique(study$SS))))
  
  ### generate a tbale with the site names whether they were sampled in each of the studies 
  
  cross <- as.data.frame(unclass(table(study$SS, study$Site_name)))
  
  #### identify how many sites have been surveyed in each study 
  for(block in rownames(cross)){
    cat(sprintf("Number of sites in study %s: %i\n", block, length(which(cross[block,] > 0))))
  }
  
  
  ## if one no combining needs to take place
  if(nrow(cross) == 1){
    cat("Source only has a single study")
    return(data)
  }  else {
    
    
    #### calculate the different combinations of studies within the source to see how their site composition differs 
    combinations <- gtools::combinations(n = length(rownames(cross)), r = 2, v = rownames(cross))
    
    ### for each combination have a look how many sites are the same and which are different 
    
    for(i in 1:nrow(combinations)){
      same <- which(c(cross[combinations[i,1],] > 0) & c(cross[combinations[i,2],] > 0))
      diff <- as.numeric(which(c(cross[combinations[i,1],] > 0) != c(cross[combinations[i,2],] > 0)))
      
      cat(sprintf("Block %s surveyed %i of the same sites as block %s\n", combinations[i,1], length(same), combinations[i,2]))
      
      if(is_empty(diff)){
        cat("No difference in sites \n")
      } else {
        cat(sprintf("sites different: %i \n", length(diff)))}
    }
    
    #### for those sites which share the same name combine them
    
    cat("Combining the same sites....\n")
    
    
    ### group by site name and identify how how studies that site is surveyed in - this is important when resolving sampling effort across combined sites
    ### the data has previously been effort corrected so that each site has an relative samplng effort of 1 - to ensure that measurements are effort corrected
    ### and each site has the same relative effort each measurement needs to take this into account, therefore the final measurement is diveded by the 
    ### number of studies it occurs in relative it to the maximum number of combined sites 
    
    combine_stud <- study %>% group_by(Site_name) %>% dplyr::mutate(num_combine = n_distinct(SS)) %>% ungroup() %>% 
      
      ## for each site sum the measurements for each species that has been recorded
      
      group_by(Site_name, Jetz_Name) %>% dplyr::mutate(Effort_Corrected_Measurement = sum(Effort_Corrected_Measurement)) %>% ungroup() %>%
      
      #### remove duplicates so each species has a single recording in the site
      
      dplyr::distinct(Site_name, Jetz_Name, .keep_all = TRUE) %>% 
      
      #### acount for number of sites combined
      
      dplyr::mutate(Effort_Corrected_Measurement = Effort_Corrected_Measurement/(num_combine/max(num_combine))) %>%
      
      ## remove column so the data is compatible with the input
      
      dplyr::select(-num_combine)
    
    ###########################################################################
    ###########################################################################
    
    
    ### Site names need to collapsed back into a single site name this is done by standard PREDICTS naming practice, each site is recoded 
    ### so will not be comparable to any previous or other versions of PREDICTS 
    
    cat("Resolving SSBS names...\n")
    
    site_names <- dplyr::distinct(combine_stud, Site_name) %>% pull() %>% as.character()
    
    j <- 1
    for(site in site_names){
      combine_stud <- combine_stud %>% dplyr::mutate(SSBS = ifelse(Site_name == site, paste(stud_combine," 1 ",j, sep = ""), paste(SSBS)))
      j <- j + 1
    }
    
    combine_stud <- combine_stud %>% dplyr::mutate(SSBS = factor(SSBS),
                                                   SS = paste(stud_combine," 1", sep = ""),
                                                   SS = factor(SS),
                                                   SSB = paste(stud_combine," 1 ", "1",sep = ""),
                                                   SSB = factor(SSB),
                                                   SSS = SSBS)
    
    ##################################################################################################
    ##################################################################################################
    
    ### now sites that have been collapsed may have different coordinates so this will need to resolved so that every record has the same site coords
    
    cat("Resolving site coordinates... \n")
    
    combine_stud <- coordinate_resolution(combine_stud)
    
    ### combine the merged and original data
    
    final_dat <- data %>% filter(!(grepl(SS, pattern = stud_combine))) %>% rbind(combine_stud)
    
    return(final_dat)
  }
}




