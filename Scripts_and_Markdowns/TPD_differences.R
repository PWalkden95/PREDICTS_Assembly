rm(list = ls())

require(tidyverse)




 PREDICTS_tpds <- readRDS("Outputs/PREDICTS_sites_tpds.rds")
 
 PREDICTS_randomisations <- readRDS("Outputs/randomisations_TPD_morpho.rds")
 
 PREDICTS_full <- readRDS("Outputs/refined_predicts.rds")
 
 PREDICTS <- PREDICTS_full %>%  ## PREDICTS data 
   dplyr::distinct(SSBS, Predominant_habitat,Use_intensity, Realm, SS) %>% ## pull out land_use type, Subregion, realm etc
   dplyr::mutate(Predominant_habitat = ifelse(grepl(Predominant_habitat, pattern = "secondary", ignore.case = TRUE), "Secondary vegetation",paste(Predominant_habitat)), 
                 Predominant_habitat = ifelse(grepl(Predominant_habitat, pattern = "primary", ignore.case = TRUE), "Primary vegetation",paste(Predominant_habitat)),
                 agri_land_use = ifelse(Predominant_habitat %in% c("Pasture","Cropland") & Use_intensity %in% c("Minimal use"), "Minimal agriculture", paste(Predominant_habitat)),
                 agri_land_use = ifelse(Predominant_habitat %in% c("Pasture","Cropland") & Use_intensity %in% c("Light use","Intense use"), "Intensive agriculture", paste(agri_land_use))) %>% 
   data.frame() 
 
 species_TPD <- readRDS("Outputs/species_tpds_morpho.rds")

 species_pools <- readRDS("Outputs/predicts_sites_species_pools.rds")

 Forage <- readRDS(file = "../PREDICTS_Taxonomy/PREDICTS_imputed_BL_traits_forage.rds") 


 TPD_LU <- data.frame(SSBS = names(PREDICTS_tpds)) %>% dplyr::left_join(PREDICTS) %>% dplyr::filter(Predominant_habitat != "Cannot decide")
 ## create a list with the land use classifications for each realm
 
 
 table(TPD_LU$Predominant_habitat, TPD_LU$Realm)
 
 
 realm_land_uses <- list()
 realm_land_uses[["Afrotropic"]] <- c("Primary vegetation", "Secondary vegetation","Plantation forest","Cropland","Pasture","Urban")
 realm_land_uses[["Australasia"]] <- c("Primary vegetation", "Secondary vegetation", "Pasture","Urban","Intensive agriculture")
 realm_land_uses[["Indo-Malay"]] <- c("Primary vegetation","Secondary vegetation","Plantation forest","Cropland")
 realm_land_uses[["Nearctic"]] <- c("Primary vegetation", "Pasture","Cropland", "Urban")
 realm_land_uses[["Neotropic"]] <- c("Primary vegetation","Secondary vegetation","Plantation forest","Pasture","Cropland","Urban")
 realm_land_uses[["Palearctic"]] <- c("Primary vegetation","Secondary vegetation", "Plantation forest","Cropland","Urban") 
 
 source("Functions/TPD_3D_plots.R") 


 land_use <- PREDICTS %>% dplyr::distinct(Predominant_habitat, .keep_all = FALSE) %>% dplyr::filter(Predominant_habitat != "Cannot decide") %>% 
   pull() %>% as.character() %>% c("Intensive agriculture","Minimal agriculture")
 
 
 realms <- TPD_LU %>% distinct(Realm) %>% pull() %>% as.character()
  
 



     dir.create("Outputs/TPD_3D_Plots")
  
  TPD_for_mapping_data <- list()
  TPD_for_mapping_data_random <- list()
 
  cumulative_legend_data <- c()
  
  for(r in realms){
    dir.create(paste("Outputs/TPD_3D_Plots",r, sep = "/"))
    dir.create(paste("Outputs/TPD_3D_Plots",r,"land_uses" ,sep = "/"))
    dir.create(paste("Outputs/TPD_3D_Plots",r,"land_use_differences" ,sep = "/"))
    dir.create(paste("Outputs/TPD_3D_Plots",r,"randomised_holes" ,sep = "/"))
  
  for(LU in land_use){

    
    if(LU %in% c("Minimal agriculture","Intensive agriculture")){
      sites_lu <- TPD_LU %>% dplyr::filter(agri_land_use == LU, Realm == r) %>% dplyr::distinct(SSBS) %>% pull()  
    } else {
      sites_lu <- TPD_LU %>% dplyr::filter(Predominant_habitat == LU, Realm == r) %>% dplyr::distinct(SSBS) %>% pull()
    }
    
  print(paste(r,LU))
    
  if(length(sites_lu) > 0){
    
    site_data <- TPD_plot_data(data = PREDICTS_tpds,sites_lu)
    legend_data <- site_data[["pl_dat"]][["prob"]]
    legend_data <- legend_data[legend_data > 0]
    cumulative_legend_data <- c(cumulative_legend_data,legend_data)
    
    
   TPD_3d_plot(data = PREDICTS_tpds, sites = sites_lu, T1lab = "Locomotion", T2lab = "Foraging", T3lab = "Body", method = "prob",
               save = TRUE, file = paste("Outputs/TPD_3D_Plots/",r,"/land_uses/",LU,"_TPD_plot.png",sep = ""),title = "",  grid = FALSE, free_limits = FALSE)
    
    # tpd_for_dat <- TPD_forage_mapping_data(data = PREDICTS_tpds, randata = PREDICTS_randomisations, for_data = Forage, sites = sites_lu)
     
   # TPD_for_mapping_data[[r]][[LU]] <- tpd_for_dat$observed
  #  TPD_for_mapping_data_random[[r]][[LU]] <- tpd_for_dat$random
    
    
}
  }    
  }
  
  
  # colours for TPD_plots
  legend_col <- c()
for(i in seq(0,1,0.0002)){
  colour <- as.numeric(quantile(cumulative_legend_data,i))
  legend_col <- c(legend_col,colour)

  }
    
  write_rds(legend_col, file = "Functions/TPD_colours.rds")
  
    find_position <- function(x,y){
  
    value <- which(x > y)
    value <- value[length(value)]
    
    return(value)   
  }
  

  
  write_rds(file = "Outputs/TPD_forage_mapping.rds", TPD_for_mapping_data)
  write_rds(file = "Outputs/TPD_forage_mapping_random.rds", TPD_for_mapping_data_random)
  
  TPD_for_mapping_data_2 <- readRDS("Outputs/TPD_forage_mapping.rds")
  TPD_for_mapping_data_random_2 <- readRDS("Outputs/TPD_forage_mapping_random.rds")

  
  
  LU <- "Primary vegetation"  
for(r in names(TPD_for_mapping_data)){
  
  dir.create(paste("Outputs/TPD_3D_Plots",r,"foraging_guilds_mapped" ,sep = "/"))
  dir.create(paste("Outputs/TPD_3D_Plots",r,"foraging_guilds_mapped" ,"null",sep = "/"))
  dir.create(paste("Outputs/TPD_3D_Plots",r,"foraging_guilds_mapped" ,"observed",sep = "/"))
  
  for(LU in names(TPD_for_mapping_data[[r]])){
    
    tpd_for_dat <- TPD_for_mapping_data[[r]][[LU]]
    
    TPD_forage_mapping_plot(data = PREDICTS_tpds,tpd_for_dat,T1lab = "", T2lab = "", T3lab = "",
                            title = "", save = TRUE,
                            s_file = paste("Outputs/TPD_3D_Plots/",r,"/foraging_guilds_mapped/","observed/",LU,"_TPD_plot.png",sep = ""),
                            animation = TRUE,
                            a_file = paste("Outputs/TPD_3D_Plots/",r,"/foraging_guilds_mapped/","observed/",LU,"_TPD_animation",sep = ""))
    
    
    tpd_for_dat_random <- TPD_for_mapping_data_random[[r]][[LU]]
    
    
    TPD_forage_mapping_plot(data = PREDICTS_tpds, tpd_for_dat_random,T1lab = "Locomotion", T2lab = "Foraging", T3lab = "Body",
                            title = paste(r,LU,"null_trophic_niche_map",sep = "_"), save = TRUE,
                            s_file = paste("Outputs/TPD_3D_Plots/",r,"/foraging_guilds_mapped/","null/",LU,"_TPD_plot.png",sep = ""),
                            animation = TRUE,
                            a_file = paste("Outputs/TPD_3D_Plots/",r,"/foraging_guilds_mapped/","null/",LU,"_TPD_animation",sep = ""))
    
  }
}
  
  
  lu_combos <- matrix(gtools::combinations(n = length(land_use), r = 2, v = land_use,set = TRUE), ncol = 2) %>%
    data.frame() %>% set_colnames(c("LU1","LU2"))
  
  
  for(r in realms){
    
  
  for(i in 1:nrow(lu_combos)){
    
    
    sites_1 <- TPD_LU %>% dplyr::filter(Realm == r, Predominant_habitat == lu_combos[i,"LU1"]) %>% dplyr::distinct(SSBS) %>% pull()
    sites_2 <- TPD_LU %>% dplyr::filter(Realm == r, Predominant_habitat == lu_combos[i,"LU2"]) %>% dplyr::distinct(SSBS) %>% pull()
    
    
    if(length(sites_1) > 0 & length(sites_2) > 0){ 
  
     TPD_Diff_plot(data = PREDICTS_tpds, sites1 = sites_1, sites2 =sites_2, T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body",
                   method = "prob", title = paste(r,"Difference",lu_combos[i,"LU1"],lu_combos[i,"LU2"], sep = ""),
                   save = TRUE, file = paste("Outputs/TPD_3D_Plots/",r,"/land_use_differences/",lu_combos[i,"LU1"],"_",lu_combos[i,"LU2"],"_difference_TPD_plot.png",sep = ""))
    
  }
  }
}
  ######################################
    ###################################
  ## calculate dissimilarity between sites and the randomised communitites 
  
    for(r in realms){
   
   
    
    for(LU in land_use){
      
      ransites <- TPD_LU %>% dplyr::filter(Predominant_habitat == LU, Realm == r) %>% dplyr::distinct(SSBS) %>% pull() 
      
      if(length(ransites) > 0){
      
      TPD_ranDiff_plot(data = PREDICTS_tpds, randata = PREDICTS_randomisations, sites = ransites, threshold = 1, T1lab = "locomotion", T2lab = "foraging",
                       T3lab = "body", title = "", method = "prob",save = TRUE, 
                       file = paste("Outputs/TPD_3D_Plots/",r,"/randomised_holes/",LU,"_random_holes_TPD_plot.png",sep = ""),observed = FALSE, grid = FALSE,free_limits = TRUE )
      
      
      }
    }
     
      
  }
  

  
  
  proportional_occupancy_df <- c()
  
  for(r in realms){
  
    land_uses <- realm_land_uses[[r]]
    
    
    for(LU in land_uses){
      
      lu_sites <- TPD_LU %>% dplyr::filter(Predominant_habitat == LU, Realm == r) %>% dplyr::distinct(SSBS) %>% dplyr::pull()
      
      
      for_mapping_data <- TPD_for_mapping_data[[r]][[LU]]
      for_mapping_data_random <- TPD_for_mapping_data_random[[r]][[LU]]
      
      
      for(g in c("In","Gr","Fr","Ne","Om")){
      prop_occ <- proportional_occupancy(data = PREDICTS_tpds, randata = PREDICTS_randomisations,
                             fordata = for_mapping_data, ranfordata = for_mapping_data_random,
                             sites = lu_sites, guild = g)
      
      
      df <- data.frame(Realm = r, Land_use = LU, Guild = g, 
                       proportional_occupancy = prop_occ$proportional_occupancy, Similarity = prop_occ$similarity)
      
      proportional_occupancy_df <- rbind(proportional_occupancy_df,df)
      
          }
    }
      
  }
  
  


  