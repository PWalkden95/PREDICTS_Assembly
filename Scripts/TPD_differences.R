rm(list = ls())

require(tidyverse)




 PREDICTS_tpds <- readRDS("Outputs/PREDICTS_sites_tpds.rds")
 PREDICTS_tpds_for <- readRDS("Outputs/PREDICTS_sites_for_tpds.rds")
 PREDICTS <- readRDS("Outputs/refined_predicts.rds") %>% dplyr::distinct(SSBS, Predominant_habitat, Use_intensity, Biome, UN_subregion, Realm) %>% 
   dplyr::mutate(Predominant_habitat = ifelse(grepl(Predominant_habitat, pattern = "secondary", ignore.case = TRUE), "Secondary vegetation",
                                              paste(Predominant_habitat))) 

 
 TPD_LU <- data.frame(SSBS = names(PREDICTS_tpds)) %>% dplyr::left_join(PREDICTS)
 
 source("Functions/TPD_3D_plots.R") 


 land_use <- PREDICTS %>% dplyr::distinct(Predominant_habitat, .keep_all = FALSE) %>% dplyr::filter(Predominant_habitat != "Cannot decide") %>% 
   pull() %>% as.character()
 

 table(TPD_LU$Predominant_habitat, TPD_LU$UN_subregion)
 table(TPD_LU$Predominant_habitat, TPD_LU$Biome)
 table(TPD_LU$Predominant_habitat, TPD_LU$Realm)
 
 primary_for <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[5], Realm == "Afrotropic") %>% pull(SSBS) %>% as.character()
 primary_non_for <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[6]) %>% pull(SSBS) %>% as.character()
 secondary <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[2]) %>% pull(SSBS) %>% as.character()
 plantation <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[1], Realm == "Afrotropic" ) %>% pull(SSBS) %>% as.character()
 urban <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[7] ) %>% pull(SSBS) %>% as.character()
 pasture <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[4] ) %>% pull(SSBS) %>% as.character()
 cropland <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[3] ) %>% pull(SSBS) %>% as.character()
 


    
  TPD_3d_plot(PREDICTS_tpds, sites = primary_for,  T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body", method = "prob")
  TPD_3d_plot(PREDICTS_tpds, sites = primary_for,  T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body", method = "percentile")
  
  TPD_3d_plot(PREDICTS_tpds, sites = primary_non_for,  T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body")
  TPD_3d_plot(PREDICTS_tpds, sites = secondary,  T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body")

  
  

  
  TPD_Diff_plot(data = PREDICTS_tpds, sites1 = primary_for, sites2 = plantation, T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body",
                method = "percentile")
  TPD_Diff_plot(data = PREDICTS_tpds, sites1 = primary_for, sites2 = plantation, T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body",
                method = "prob")

  
