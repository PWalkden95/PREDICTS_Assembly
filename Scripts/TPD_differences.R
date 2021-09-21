rm(list = ls())


 PREDICTS_tpds <- readRDS("Outputs/PREDICTS_sites_tpds.rds")
 PREDICTS_tpds_for <- readRDS("Outputs/PREDICTS_sites_for_tpds.rds")
 PREDICTS <- readRDS("Outputs/refined_predicts.rds") %>% dplyr::distinct(SSBS, Predominant_habitat, Use_intensity) %>% 
   dplyr::mutate(Predominant_habitat = ifelse(grepl(Predominant_habitat, pattern = "Primary"), "Primary", 
                                              paste(Predominant_habitat)),
                 Predominant_habitat = ifelse(grepl(Predominant_habitat, pattern = "secondary", ignore.case = TRUE), "Secondary vegetation",
                                              paste(Predominant_habitat))) 

 
 TPD_LU <- data.frame(SSBS = names(PREDICTS_tpds)) %>% dplyr::left_join(PREDICTS)
 
 source("Functions/TPD_3D_plots.R") 


 land_use <- PREDICTS %>% dplyr::distinct(Predominant_habitat, .keep_all = FALSE) %>% dplyr::filter(Predominant_habitat != "Cannot decide") %>% 
   pull() %>% as.character()
 

 LU <- land_use[2]
 
  dat <- TPD_LU %>% dplyr::filter(Predominant_habitat == LU ) %>% pull(SSBS) %>% as.character() 

    
  TPD_3d_plot(PREDICTS_tpds, sites = plantation,  T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body")

  
  primary <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[5] ) %>% pull(SSBS) %>% as.character()
  secondary <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[2] ) %>% pull(SSBS) %>% as.character()
  plantation <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[1] ) %>% pull(SSBS) %>% as.character()
  urban <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[6] ) %>% pull(SSBS) %>% as.character()

  
  TPD_Diff_plot(data = PREDICTS_tpds, sites1 = primary, sites2 = plantation, T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body")
 