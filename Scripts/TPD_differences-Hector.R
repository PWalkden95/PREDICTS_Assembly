rm(list = ls())

require(tidyverse)




 PREDICTS_tpds <- readRDS("Outputs/PREDICTS_sites_tpds.rds")
 randomisations_tpds <- readRDS("Outputs/randomisations_TPD_morpho.rds")
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
 
 primary_for <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[5], Realm == "Australasia") %>% pull(SSBS) %>% as.character()
 primary_non_for <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[6], Realm == "Australasia") %>% pull(SSBS) %>% as.character()
 secondary <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[2], Realm == "Australasia") %>% pull(SSBS) %>% as.character()
 plantation <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[1]) %>% pull(SSBS) %>% as.character()
 urban <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[7] ) %>% pull(SSBS) %>% as.character()
 pasture <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[4] ) %>% pull(SSBS) %>% as.character()
 cropland <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[3] ) %>% pull(SSBS) %>% as.character()
 


    
  TPD_3d_plot(PREDICTS_tpds, sites = TPD_LU$SSBS[-1],  T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body", method = "prob")
  TPD_3d_plot(PREDICTS_tpds, sites = primary_for,  T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body", method = "percentile")
  
  TPD_3d_plot(PREDICTS_tpds, sites = primary_non_for,  T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body", method = "prob")
  TPD_3d_plot(PREDICTS_tpds, sites = primary_for,  T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body", method = "prob")

  
  
  
  TPD_Diff_plot(data = PREDICTS_tpds, sites1 = primary_for, sites2 = primary_non_for, T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body",
                method = "prob")
  TPD_Diff_plot(data = PREDICTS_tpds, sites1 = primary_for, sites2 = secondary, T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body",
                method = "prob")

  site_dissim <- Calc_dissim(PREDICTS_tpds, sites1 = primary_for, sites2 = primary_non_for)
  
  
  TPD_ranDiff_plot(data = PREDICTS_tpds,randata = randomisations_tpds,site = TPD_LU$SSBS[-1],T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body",
                   method = "prob")
  

  
  