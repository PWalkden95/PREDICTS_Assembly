require(tidyverse)
require(gtools)

combine_study <- function(data, stud_combine) {
  
  study <- data %>% filter(grepl(SS, pattern = stud_combine)) %>% droplevels()
  
  cat(sprintf("Number of blocks in study: %i\n", length(unique(study$SS))))
  
  cross <- as.data.frame(unclass(table(study$SS, study$Site_name)))
  
  for(block in rownames(cross)){
    cat(sprintf("Number of sites in block %s: %i\n", block, length(which(cross[block,] > 0))))
  }
  
  if(nrow(cross) == 1){
    cat("Study only has a single block")
    return(data)
  }  else {
    combinations <- gtools::combinations(n = length(rownames(cross)), r = 2, v = rownames(cross))
    
    for(i in 1:nrow(combinations)){
      same <- which(c(cross[combinations[1],] > 0) == c(cross[combinations[2],] > 0))
      diff <- as.numeric(which(c(cross[combinations[1],] > 0) != c(cross[combinations[2],] > 0)))
      
      cat(sprintf("Block %s surveyed %i of the same sites as block %s\n", combinations[1], length(same), combinations[2]))
      
      if(is_empty(diff)){
        cat("No difference in sites")
      } else {
        cat(sprintf("sites different: %i \n", length(diff)))}
    }
    
    cat("Combining the same sites....\n")
    
    combine_stud <- study %>% group_by(Site_name, Jetz_Name) %>%
      dplyr::mutate(Effort_Corrected_Measurement = sum(Effort_Corrected_Measurement)) %>% ungroup() %>%
      dplyr::distinct(Site_name, Jetz_Name, .keep_all = TRUE) 
    
    cat("Resolving SSBS names...\n")
    
    site_names <- dplyr::distinct(combine_stud, Site_name) %>% pull() %>% as.character()
    
    j <- 1
    for(site in site_names){
      combine_stud <- combine_stud %>% dplyr::mutate(SSBS = ifelse(Site_name == site, paste(stud_combine," 1 ",j, sep = ""), paste(SSBS)))
      j <- j + 1
    }
    
    combine_stud <- combine_stud %>% dplyr::mutate(SSBS = factor(SSBS),
                                                   SS = paste(stud_combine," 1"),
                                                   SS = factor(SS),
                                                   SSB = paste(stud_combine," 1 ", "1",sep = ""),
                                                   SSB = factor(SSB),
                                                   SSS = SSBS)
    
    
    combine_stud <- data %>% filter(!(grepl(SS, pattern = stud_combine))) %>% rbind(combine_stud)
    
    return(combine_stud)
  }
}