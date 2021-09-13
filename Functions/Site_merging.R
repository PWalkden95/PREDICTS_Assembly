require(tidyverse)



merge_by <- function(data, by.var){
  
  by.var <- enquo(by.var)
  
  output <- data %>% group_by(!!by.var, Jetz_Name) %>% dplyr::mutate(Effort_Corrected_Measurement = sum(Effort_Corrected_Measurement)) %>% ungroup() %>%
    
    #### remove duplicates so each species has a single recording in the site
    
    dplyr::distinct(!!by.var, Jetz_Name, .keep_all = TRUE) %>% 
    
    #### account for number of sites combined
    
    dplyr::mutate(Effort_Corrected_Measurement = Effort_Corrected_Measurement/(num_combine/max(num_combine)))
   
  
  return(output)
}


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





site_merge <- function(data,merge_data){
  
  fun_data <- data
  
  for(study in merge_data$Study){
  
    cat(sprintf("Checking study %s \n", study))
    
  study_dat <- fun_data %>% dplyr::filter(SS == study)
  merge_code <- merge_data %>% dplyr::filter(Study == study) %>% dplyr::pull(Merge_code)
  
  if(merge_code == "none"){
    cat(sprintf("No sites in %s need to be merged \n",study))
  } else {
    cat(sprintf("Sites in %s need to be merged by %s \n", study, merge_code))
    
    if(merge_code == "LUI"){
    
      ############# create a column by which to merge
      
      LUI_dat <- study_dat %>% dplyr::mutate(LUI = paste(Predominant_habitat, Use_intensity, sep = "_")) %>%
        group_by(LUI) %>% dplyr::mutate(num_combine = n_distinct(SSBS))
      
      
      for(LU in unique(LUI_dat$LUI)){
        
        num_sites <- LUI_dat %>% dplyr::filter(LUI == LU) %>% dplyr::distinct(num_combine) %>% dplyr::pull() %>% as.numeric()
        
        cat(sprintf("Merging %i %s sites \n", num_sites, LU))
      }
        
        
      LUI_com <- merge_by(data = LUI_dat, by.var = LUI)
  
      
      ### going to want to rename the sites and recode them.
      
      cat("Resolving site names... \n")
      
    j <- 1  
      for(LU in unique(LUI_com$LUI)){
        LUI_com <- LUI_com %>% dplyr::mutate(SSBS = ifelse(LUI == LU, paste(study ,j), paste(SSBS)))
        j <- j + 1
      }
    
    combine_dat <- LUI_com %>% dplyr::mutate(SSBS = factor(SSBS),
                                                   SS = factor(SS),
                                                   SSB = factor(SSB),
                                                   SSS = SSBS) %>% 
      dplyr::select(-LUI) %>% dplyr::select(-num_combine)
  
    
    }
    
    #######################################################################
    #######################################################################
    # IF THE STUDY HAS TO BE MERGED BY SITE NAME
    
    if(merge_code == "site_names"){
      
      match <- merge_data %>% dplyr::filter(Study == study) %>% distinct(by) %>% pull() %>% strsplit(match, split = ";") %>% unlist()
      
      site_dat <- study_dat %>% dplyr::mutate(LUI = paste(Predominant_habitat,Use_intensity, sep = "_"))
      site_dat$Code <- NA
      site_dat$num_combine <- NA
      
      
      for(code in match){
        site_dat$Code <- ifelse(grepl(site_dat$Site_name, pattern = code, fixed = TRUE), code, site_dat$Code)
      }
      
    
      for(code in match){
        land <- site_dat %>% dplyr::filter(Code == code) %>% distinct(LUI) %>% pull()
        
        for(LU in land){
          
          num <- site_dat %>% dplyr::filter(Code == code & LUI == LU) %>% distinct(SSBS) %>% nrow() %>% as.numeric()
          
          cat(sprintf("%s named sites in %s study have %i sites of %s LUI that are to be merged \n", code, study, num, LU))
          
          site_dat$num_combine <- ifelse(grepl(site_dat$Code, pattern = code,fixed = TRUE) & grepl(site_dat$LUI, pattern = LU, fixed = TRUE), num,site_dat$num_combine)
          site_dat$Code <- ifelse(grepl(site_dat$Code, pattern = code,fixed = TRUE) & grepl(site_dat$LUI, pattern = LU, fixed = TRUE), paste(site_dat$Code, LU),site_dat$Code)
          
          
        }
        
      }
      
      site_com <- merge_by(site_dat, Code)
      
      ### going to want to rename the sites and recode them.
      
      cat("Resolving site names... \n")
      
      j <- 1  
      for(code in unique(site_com$Code)){
        site_com <- site_com %>% dplyr::mutate(SSBS = ifelse(Code == code, paste(study ,j), paste(SSBS)))
        j <- j + 1
      }
      
      LUI_test <- site_com %>% dplyr::mutate(LUI = paste(Predominant_habitat,Use_intensity,sep = "_")) %>% dplyr::distinct(Code, LUI)
      
      if(any(duplicated(LUI_test$Code))){
        warning(paste("multiple land use classifications within site please merging, may need to be more stringent", study))
      }
      
      combine_dat <- site_com %>% dplyr::mutate(SSBS = factor(SSBS),
                                               SS = factor(SS),
                                               SSB = factor(SSB),
                                               SSS = SSBS) %>% 
        dplyr::select(-Code) %>% dplyr::select(-num_combine) %>% dplyr::select(-LUI)
      
      
      
    } ############## site name merge ends here 
    
    ###############################################
    ###############################################
    ### NEXT IF YOU WANT TO MERGE BY SITE NUMBER 
    
    if(merge_code == "site_number"){
      
      site_dat <- study_dat %>% dplyr::mutate(LUI = paste(Predominant_habitat,Use_intensity, sep = "_"))
      
      match <- merge_data %>% dplyr::filter(Study == study) %>% distinct(by) %>% pull() %>% strsplit(match, split = ";") %>% unlist()
      
      match_list <- c(rep(list(NA), length(match)))
     
       for(i in 1:length(match)){
     
        numbers <- c()
        
        if(grepl(match[i], pattern = "c")){
          match[i] <-  gsub(match[i],pattern = "c",replacement = "")
          match[i] <-  gsub(match[i],pattern = "\\(",replacement = "")
          match[i] <-  gsub(match[i],pattern = ")",replacement = "")
          
          new <- unlist(strsplit(match[i], split = ","))
          
          
          for(num in new){
            if(grepl(num, pattern = ":")){
              chr_num <- unlist(strsplit(num, split = ":"))
              numbers <- c(numbers,c(as.numeric(chr_num[1]):as.numeric(chr_num[2])))
            } else {
              numbers <- c(numbers, as.numeric(num))
            }
          }
        } else {
          if(grepl(match[i], pattern = ":")){
            chr_num <- unlist(strsplit(match[i], split = ":"))
            numbers <- c(numbers,c(as.numeric(chr_num[1]):as.numeric(chr_num[2])))
          } else {
            numbers <- as.numeric(match[i])
          }
        }
        match_list[[i]] <- numbers
      }
  
      site_dat$Code <- NA
      site_dat$num_combine <- NA
      
      for(i in 1:length(match_list)){
        
        site_dat$Code <- ifelse(site_dat$Site_number %in% match_list[[i]], i, site_dat$Code)
        site_dat$num_combine <- ifelse(site_dat$Site_number %in% match_list[[i]], length(match_list[[i]]), site_dat$num_combine)
        
       
      }
      
      for(code in unique(site_dat$Code)){
        land <- site_dat %>% dplyr::filter(Code == code) %>% distinct(LUI) %>% pull()
        
        for(LU in land){
          
          num <- site_dat %>% dplyr::filter(Code == code & LUI == LU) %>% distinct(SSBS) %>% nrow()
          
          cat(sprintf("%s named sites in %s study have %i sites of %s LUI that are to be merged \n", code, study, num, LU))
          
          site_dat$Code <- ifelse(grepl(site_dat$Code, pattern = code) & grepl(site_dat$LUI, pattern = LU, fixed = TRUE), paste(site_dat$Code, LU),site_dat$Code)
          site_dat$num_combine <- ifelse(grepl(site_dat$Code, pattern = code) & grepl(site_dat$LUI, pattern = LU, fixed = TRUE), num,site_dat$num_combine)
          
        }
        
      }
      
      
      site_com <- merge_by(site_dat, by.var = Code)
      
      
      ### going to want to rename the sites and recode them.
      
      cat("Resolving site names... \n")
      
      j <- 1  
      for(code in unique(site_com$Code)){
        site_com <- site_com %>% dplyr::mutate(SSBS = ifelse(Code == code, paste(study ,j), paste(SSBS)))
        j <- j + 1
      }
      
      LUI_test <- site_com %>% dplyr::mutate(LUI = paste(Predominant_habitat,Use_intensity,sep = "_")) %>% dplyr::distinct(Code, LUI)
      
      if(any(duplicated(LUI_test$Code))){
        warning(paste("multiple land use classifications within site please merging, may need to be more stringent", study))
      }
      
      combine_dat <- site_com %>% dplyr::mutate(SSBS = factor(SSBS),
                                                SS = factor(SS),
                                                SSB = factor(SSB),
                                                SSS = SSBS) %>% 
        dplyr::select(-Code) %>% dplyr::select(-num_combine) %>% dplyr::select(-LUI)
                 
      
    } #################### site number merge ends here 
    

    if(merge_code == "blocks"){
      
      block_dat <- study_dat %>% dplyr::mutate(LUI = paste(Predominant_habitat,Use_intensity, sep = "_"))
      block_dat$Code <- NA
      block_dat$num_combine <- NA
      
      blocks <- as.character(unique(block_dat$Block))
      
      for(blo in blocks){
        
      land <- block_dat %>% dplyr::filter(Block == blo) %>% distinct(LUI) %>% pull()
      
      for(LU in land){
        
        num <- block_dat %>% dplyr::filter(Block == blo & LUI == LU) %>% distinct(SSBS) %>% nrow()
        
      cat(sprintf("%s block in %s study has %i sites of %s LUI that are to be merged \n", blo, study, num, LU))
          
      block_dat$Code <- ifelse(grepl(block_dat$Block, pattern = blo) & grepl(block_dat$LUI, pattern = LU, fixed = TRUE), paste(blo,LU,sep = "_"),block_dat$Code)
      block_dat$num_combine <- ifelse(grepl(block_dat$Block, pattern = blo) & grepl(block_dat$LUI, pattern = LU,fixed = TRUE), num,block_dat$num_combine)
      }
      
      }
      
      block_com <- merge_by(block_dat, Code)
      
      
      ### going to want to rename the sites and recode them.
      
      cat("Resolving site names... \n")
      
      
      for(blo in blocks){
      
      codes <- unique(block_com$Code)[grep(unique(block_com$Code), pattern = blo)]
        
      j <- 1  
      for(code in codes){
        block_com <- block_com %>% dplyr::mutate(SSBS = ifelse(Code == code, paste(study ,blo, j), paste(SSBS)))
        j <- j + 1
      }
      
      }
      
      
      
      combine_dat <- block_com %>% dplyr::mutate(SSBS = factor(SSBS),
                                                SS = factor(SS),
                                                SSB = factor(SSB),
                                                SSS = SSBS) %>% 
        dplyr::select(-Code) %>% dplyr::select(-num_combine) %>% dplyr::select(-LUI)
      
      
      
      
    }############################### blocks merge ends here 
    
    if(merge_code == "block_names"){
      
      block_dat <- study_dat  %>% dplyr::mutate(LUI = paste(Predominant_habitat,Use_intensity, sep = "_"))
      
      match <- merge_data %>% dplyr::filter(Study == study) %>% distinct(by) %>% pull() %>% strsplit(match, split = ";") %>% unlist()
    
      
      match_list <- c(rep(list(NA), length(match)))
      
      for(i in 1:length(match)){
        

        if(grepl(match[i], pattern = "c")){
          match[i] <-  gsub(match[i],pattern = "c",replacement = "")
          match[i] <-  gsub(match[i],pattern = "\\(",replacement = "")
          match[i] <-  gsub(match[i],pattern = ")",replacement = "")
          
          block_names <- unlist(strsplit(match[i], split = ","))
        } else {
          block_names <- match[i]
        }
        match_list[[i]] <- block_names
        }
        
      

      block_dat$Block_2 <- NA
      
      for(j in 1:length(match_list)){
        
        if(length(match_list[[j]]) > 1){
          
          block_code <- ""
          for(k in 1:length(match_list[[j]])){
            block_code <- paste(block_code, match_list[[j]][k], sep = "")
          }
        } else {
          block_code <- match_list[[j]][1]
        }
        
      
        for(k in 1:length(match_list[[j]])){
          block_dat$Block_2 <- ifelse(grepl(block_dat$Block, pattern = match_list[[j]][k]), block_code, block_dat$Block_2)
        }
          
        
      }
      
      block_dat$Block <- block_dat$Block_2
      block_dat <- block_dat %>% dplyr::select(-Block_2)
      block_dat$num_combine <- NA
      block_dat$Code <- NA
      
      blocks <- unique(block_dat$Block)
      
      for(blo in blocks){
        
        land <- block_dat %>% dplyr::filter(Block == blo) %>% distinct(LUI) %>% pull()
        
        for(LU in land){
          
          num <- block_dat %>% dplyr::filter(Block == blo & LUI == LU) %>% distinct(SSBS) %>% nrow()
          
          cat(sprintf("%s block in %s study has %i sites of %s LUI that are to be merged \n", blo, study, num, LU))
          
          block_dat$Code <- ifelse(grepl(block_dat$Block, pattern = blo) & grepl(block_dat$LUI, pattern = LU,fixed = TRUE), paste(blo,LU,sep = "_"),block_dat$Code)
          block_dat$num_combine <- ifelse(grepl(block_dat$Block, pattern = blo) & grepl(block_dat$LUI, pattern = LU,fixed = TRUE), num,block_dat$num_combine)
        }
        
      }
      
      block_com <- merge_by(block_dat, Code)
    
      
      for(blo in blocks){
        
        codes <- unique(block_com$Code)[grep(unique(block_com$Code), pattern = blo)]
        
        block_com <- block_com %>% dplyr::mutate(SSB = ifelse(Block == blo, paste(study, blo), paste(SSB)))
        
        
        j <- 1  
        for(code in codes){
          block_com <- block_com %>% dplyr::mutate(SSBS = ifelse(Code == code, paste(study ,blo, j), paste(SSBS)))
          j <- j + 1
        }
        
      }
      
      combine_dat <- block_com %>% dplyr::mutate(SSBS = factor(SSBS),
                                                 SS = factor(SS),
                                                 SSB = factor(SSB),
                                                 SSS = SSBS) %>% 
        dplyr::select(-Code) %>% dplyr::select(-num_combine) %>% dplyr::select(-LUI)
      
      
      
    } ### block names merge ends here
   
    if(merge_code == "habitat"){
      
      habitat_dat <- study_dat %>% dplyr::mutate(LUI = paste(Predominant_habitat,Use_intensity, sep = "_"))
      habitat_dat$Code <- NA
      habitat_dat$num_combine <- NA
      
      habitat <- as.character(unique(habitat_dat$Habitat_as_described))
      
      
      for(hab in habitat){
        
        land <- habitat_dat %>% dplyr::filter(Habitat_as_described == hab) %>% distinct(LUI) %>% pull()
        
        for(LU in land){
          
          num <- habitat_dat %>% dplyr::filter(Habitat_as_described == hab & LUI == LU) %>% distinct(SSBS) %>% nrow()
          
          cat(sprintf("%s block in %s study has %i sites of %s LUI that are to be merged \n", hab, study, num, LU))
          
          habitat_dat$Code <- ifelse(grepl(habitat_dat$Habitat_as_described, pattern = hab, fixed = TRUE) & grepl(habitat_dat$LUI, pattern = LU, fixed = TRUE),
                                     paste(hab,LU,sep = "_"),habitat_dat$Code)
          habitat_dat$num_combine <- ifelse(grepl(habitat_dat$Habitat_as_described, pattern = hab, fixed = TRUE) & grepl(habitat_dat$LUI, pattern = LU,fixed = TRUE),
                                            num,habitat_dat$num_combine)
        }
        
      }
      
      habitat_com <- merge_by(habitat_dat, Code)
      
      
      j <- 1  
      for(code in unique(habitat_com$Code)){
        habitat_com <- habitat_com %>% dplyr::mutate(SSBS = ifelse(Code == code, paste(study, j), paste(SSBS)))
        j <- j + 1
      }
      
      
      combine_dat <- habitat_com %>% dplyr::mutate(SSBS = factor(SSBS),
                                                 SS = factor(SS),
                                                 SSB = factor(SSB),
                                                 SSS = SSBS) %>% 
        dplyr::select(-Code) %>% dplyr::select(-num_combine) %>% dplyr::select(-LUI)
      
      
      
    } ### end of habitat merge
    
    combine_dat <- coordinate_resolution(combine_dat)
    
    
    fun_data <- fun_data %>% dplyr::filter(SS != study) %>% rbind(combine_dat)
  
  }
  
  
  
  }
  
  return(fun_data)
  
}


