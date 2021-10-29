rm(list = ls())

require(tidyverse)




 PREDICTS_tpds <- readRDS("Outputs/PREDICTS_sites_tpds.rds")
 PREDICTS_tpds_for <- readRDS("Outputs/PREDICTS_sites_for_tpds.rds")
 PREDICTS_randomisations <- readRDS("Outputs/randomisations_TPD_morpho.rds")
 PREDICTS <- readRDS("Outputs/refined_predicts.rds") %>% dplyr::distinct(SSBS, Predominant_habitat, Use_intensity, Biome, UN_subregion, Realm) %>% 
   dplyr::mutate(Predominant_habitat = ifelse(grepl(Predominant_habitat, pattern = "secondary", ignore.case = TRUE), "Secondary vegetation",
                                              paste(Predominant_habitat))) %>%  data.frame()
 species_TPD <- readRDS("Outputs/species_tpds_morpho.rds")

 
 TPD_LU <- data.frame(SSBS = names(PREDICTS_tpds)) %>% dplyr::left_join(PREDICTS) %>% dplyr::filter(Predominant_habitat != "Cannot decice",
                                                                                                    Use_intensity != "Cannot decide")
 
 primary <- TPD_LU %>% dplyr::filter(Predominant_habitat == "Primary forest", Realm == "Afrotropic") %>% distinct(SSBS) %>% pull()
 secondary <- TPD_LU %>% dplyr::filter(Predominant_habitat == "Secondary vegetation", Realm == "Afrotropic") %>% distinct(SSBS) %>% pull()
 
 
 
 source("Functions/TPD_3D_plots.R") 


 land_use <- PREDICTS %>% dplyr::distinct(Predominant_habitat, .keep_all = FALSE) %>% dplyr::filter(Predominant_habitat != "Cannot decide") %>% 
   pull() %>% as.character()
 

 table(TPD_LU$Predominant_habitat, TPD_LU$UN_subregion)
 table(TPD_LU$Predominant_habitat, TPD_LU$Biome)
 table(TPD_LU$Predominant_habitat, TPD_LU$Realm)
 
 
 realms <- TPD_LU %>% distinct(Realm) %>% pull() %>% as.character()
  
 
 
  ##### check whether the Biogeographic-Land-use type combinations are sufficiently sampled 

  
  combinations <- expand.grid(land_use,realms) %>% set_colnames(c("land_use","realm"))
  
  
  
  Check_list <- list()
  for(i in 1:nrow(combinations)){
  
    Land_use <- as.character(combinations[i,"land_use"])
    Realm <- as.character(combinations[i,"realm"])
    
    list <-   TPD_site_check(LU = Land_use , realm = Realm)
    
    Check_list[[i]] <- list
    names(Check_list)[i] <- paste(Land_use, Realm, sep = "/")
  }
###########################################################
  ##########################################################
  
  
  for(i in 1:length(Check_list)){
    if(is.null(Check_list[[i]]$ggplots)){
      next()
    }
    
    figure <- Check_list[[i]]$ggplots
    figure <- annotate_figure(figure,
                    top = text_grob(paste(names(Check_list)[i]), color = "red", face = "bold", size = 14))
    
    plot(figure)
    }
  ######################################################################################################  
  ######################################################################################################
  ######################################################################################################
  
  
  # evaluating the FD metric plots some realm/land_use type combinTIONS HAVE suffiecients sites to compare
  
  # Afrotropic: Cropland/Primary/Secondary
  # Australasia; Secondary and primary
  # Indo-Malay: Plantation forest/primary forest/ Secondary veg
  # Nearctic: Primary forest/non-forest maybe pasture/cropland
  # Neotropic: Pasture, Plantation, Primary forest/Secondary
  # Palearctic: Pasture/plantation forest/primary forest/secondary vegetation
  
  
  ## so lets just have a look at the Neotropics as an example
  
  dir.create("Outputs/TPD_3D_Plots")
  
  for(r in realms){
  dir.create(paste("Outputs/TPD_3D_Plots",r, sep = "/"))  
  }
  
  
  dissim_list <- list()
  for(r in realms){
    
    if(r == "Afrotropic"){
      land_uses <- c("Cropland","Primary forest","Primary non-forest", "Secondary vegetation","Plantation forest")
    }
    if(r == "Australasia"){
      land_uses <- c("Primary forest","Primary non-forest", "Secondary vegetation")
    }
    if(r == "Indo-Malay"){
      land_uses <- c("Primary forest","Plantation forest")
    }
    if(r == "Nearctic"){
      land_uses <- c("Cropland","Primary forest","Primary non-forest","Pasture")
    }
    if(r == "Neotropic"){
      land_uses <- c("Cropland","Primary forest","Pasture","Plantation forest","Secondary vegetation")
    }
    if(r == "Palearctic"){
      land_uses <- c("Primary forest","Pasture","Plantation forest","Secondary vegetation")
    }
    
  LU_combo <- matrix(gtools::combinations(n = length(land_uses), r = 2, v = land_uses,set = TRUE), ncol = 2) %>% data.frame() %>% set_colnames(c("LU1","LU2"))
  
  
  array_dimnames <- list(c(land_uses),c(land_uses),c("Dissimilarity","Shared","Not_shared"))
  LU_array <- array(rep(NA,(length(land_uses)^2)*3),dim = c(length(land_uses),length(land_uses),3), dimnames = array_dimnames)
  
  for(i in 1:nrow(LU_combo)){
    
    
    sites_1 <- TPD_LU %>% dplyr::filter(Realm == r, Predominant_habitat == LU_combo[i,"LU1"]) %>% dplyr::distinct(SSBS) %>% pull()
    sites_2 <- TPD_LU %>% dplyr::filter(Realm == r, Predominant_habitat == LU_combo[i,"LU2"]) %>% dplyr::distinct(SSBS) %>% pull()
    
    dissimilarity <- Calc_dissim(data = PREDICTS_tpds,sites1 = sites_1, sites2 = sites_2)
    
    LU_array[LU_combo[i,"LU1"],LU_combo[i,"LU2"],"Dissimilarity"] <- dissimilarity[["dissim"]][["dissimilarity"]]
    LU_array[LU_combo[i,"LU1"],LU_combo[i,"LU2"],"Shared"] <- dissimilarity[["dissim"]][["P_shared"]]
    LU_array[LU_combo[i,"LU1"],LU_combo[i,"LU2"],"Not_shared"] <- dissimilarity[["dissim"]][["P_non_shared"]]
    
    LU_array[LU_combo[i,"LU2"],LU_combo[i,"LU1"],"Dissimilarity"] <- dissimilarity[["dissim"]][["dissimilarity"]]
    LU_array[LU_combo[i,"LU2"],LU_combo[i,"LU1"],"Shared"] <- dissimilarity[["dissim"]][["P_shared"]]
    LU_array[LU_combo[i,"LU2"],LU_combo[i,"LU1"],"Not_shared"] <- dissimilarity[["dissim"]][["P_non_shared"]]
  
    # TPD_Diff_plot(data = PREDICTS_tpds, sites1 = sites_1, sites2 =sites_2, T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body",
    #               method = "prob", title = paste("3D_TPD_Plot",LU_combo[i,"LU1"],LU_combo[i,"LU2"], sep = "_"))
    # rgl.snapshot(paste("Outputs/TPD_3D_Plots/",r,"/",LU_combo[i,"LU1"],"_",LU_combo[i,"LU2"],"_difference_", "3dplot.png",sep = ""), fmt = 'png')
    
  }
    
  dissim_list[[r]][["dissim_array"]] <- LU_array
  
  
  }

  ######################################
    ###################################
  ## calculate dissimilarity between sites and the randomised communitites 
  
  
  ran_dissim_list <- list()
  i <- 1
  for(r in realms){
   
    if(r == "Afrotropic"){
      land_uses <- c("Cropland","Primary forest","Primary non-forest", "Secondary vegetation","Plantation forest")
    }
    if(r == "Australasia"){
      land_uses <- c("Primary forest","Primary non-forest", "Secondary vegetation")
    }
    if(r == "Indo-Malay"){
      land_uses <- c("Primary forest","Plantation forest")
    }
    if(r == "Nearctic"){
      land_uses <- c("Cropland","Primary forest","Primary non-forest","Pasture")
    }
    if(r == "Neotropic"){
      land_uses <- c("Cropland","Primary forest","Pasture","Plantation forest","Secondary vegetation")
    }
    if(r == "Palearctic"){
      land_uses <- c("Primary forest","Pasture","Plantation forest","Secondary vegetation")
    }
   
    
    for(l in land_uses){
      
      ransites <- TPD_LU %>% dplyr::filter(Predominant_habitat == l, Realm == r) %>% dplyr::distinct(SSBS) %>% pull() 
      
      random_dissim <- Calc_dissim_random(data = PREDICTS_tpds, randata = PREDICTS_randomisations, sites = ransites, threshold = 0.95)
      
      ran_dissim_list[i] <- random_dissim
      names(ran_dissim_list)[i] <- paste(r,l,sep = "_")
      
      i <- i +1
    
      }
    
     
      
  }
  
  ransites <- TPD_LU %>% dplyr::filter(Predominant_habitat == "Primary non-forest", Realm == "Afrotropic") %>% dplyr::distinct(SSBS) %>% pull()
  
  TPD_ranDiff_plot(data = PREDICTS_tpds, randata = PREDICTS_randomisations, site = ransites, threshold = 0.95, T1lab = "locomotion", T2lab = "foraging",
                   T3lab = "body", title = "hehe", method = "prob")
  
  
  
  
  
  
  