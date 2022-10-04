{r species difference, fig.width= 10}


diff_species <- function(realm){
  
  dif_sp <- list()
  
  land_uses <- realm_land_uses[[realm]]
  
  ## get each combination of land uses within the realm
  
  LU_combo <- matrix(gtools::combinations(n = length(land_uses), r = 2, v = land_uses,set = TRUE), ncol = 2) %>% data.frame() %>% set_colnames(c("LU1","LU2")) %>% reorder_combinations()
  
  
  for(i in 1:nrow(LU_combo)){
    
    
    sites_1 <- TPD_LU %>% dplyr::filter(Realm == realm, Predominant_habitat == LU_combo[i,"LU1"]) %>% dplyr::distinct(SSBS) %>% pull()
    sites_2 <- TPD_LU %>% dplyr::filter(Realm == realm, Predominant_habitat == LU_combo[i,"LU2"]) %>% dplyr::distinct(SSBS) %>% pull()
    
    
    loss_gain <- TPD_species_occupancy(data = PREDICTS_tpds, randata = PREDICTS_randomisations,sites1 = sites_1, sites2 = sites_2)
    
    gain_sp <- split_sp_func(pool_sp(loss_gain$gain[["gain_species"]], sites = c(sites_2), pool = "0.9"))
    lost_sp <- split_sp_func(pool_sp(loss_gain$lost[["lost_species"]], sites = c(sites_1), pool = "0.9"))
    
    
    dif_sp[[paste(LU_combo[i,"LU1"],LU_combo[i,"LU2"], sep = "_")]]$gain <- gain_sp  
    dif_sp[[paste(LU_combo[i,"LU1"],LU_combo[i,"LU2"], sep = "_")]]$gain_vol <- loss_gain$gain$gain_volume
    
    dif_sp[[paste(LU_combo[i,"LU1"],LU_combo[i,"LU2"], sep = "_")]]$lost <- lost_sp
    dif_sp[[paste(LU_combo[i,"LU1"],LU_combo[i,"LU2"], sep = "_")]]$lost_vol <- loss_gain$lost$lost_volume
    
  }
  
  return(dif_sp)
  
}

options(future.globals.maxSize = 8000 * 1024^2)

if(any(grepl(list.files("../Outputs"), pattern = "pairwise_sp_difs.rds"))){
  pairwise_sp_difs <- readRDS(grep(list.files("../Outputs",full.names = TRUE),pattern = "pairwise_sp_difs.rds", value = TRUE))
} else{
  plan(multicore(workers = 8))
  
  realm_list <- lapply(realms, paste)
  
  pairwise_sp_difs <- future_lapply(realm_list, diff_species, future.label = TRUE)
  
  closeAllConnections()
  
  write_rds("Outputs/pairwise_sp_difs.rds", pairwise_sp_difs)
}


pairwise_diff_plot_sp <- c()

for( r in realms){
  
  land_uses <- realm_land_uses[[r]]
  
  
  ## get each combination of land uses within the realm
  
  LU_combo <- matrix(gtools::combinations(n = length(land_uses), r = 2, v = land_uses,set = TRUE), ncol = 2) %>% data.frame() %>% set_colnames(c("LU1","LU2")) %>% reorder_combinations()
  
  
  for(i in 1:nrow(LU_combo)){
    
    LU_sp_difs <- pairwise_sp_difs[[r]][[i]]
    
    gain_sp <- data.frame(Birdlife_Name = LU_sp_difs[["gain"]], sp_differences = "gain", Land_use = paste(LU_combo[i,1], LU_combo[i,2],sep = "/"), Realm = r, vol = LU_sp_difs[["gain_vol"]])
    loss_sp <- data.frame(Birdlife_Name = LU_sp_difs[["lost"]], sp_differences = "lost", Land_use = paste(LU_combo[i,1], LU_combo[i,2],sep = "/"), Realm = r,vol = LU_sp_difs[["lost_vol"]])
    
    sp_differences <- rbind(gain_sp,loss_sp) %>% dplyr::left_join(Forage[,c("Birdlife_Name", "Trophic_Level", "Trophic_Niche", "Foraging_Niche")], by = "Birdlife_Name")
    
    
    
    sp_differences$Trophic_Niche <- factor(sp_differences$Trophic_Niche)
    
    
    pairwise_diff_plot_sp <- rbind(pairwise_diff_plot_sp,sp_differences)
    
    
  }
}

priority_combos <- c( "Primary vegetation/Secondary vegetation",
                      "Primary vegetation/Plantation forest",
                      "Primary vegetation/Pasture",
                      "Primary vegetation/Cropland",
                      "Primary vegetation/Urban",
                      "Secondary vegetation/Plantation forest",
                      "Secondary vegetation/Pasture",
                      "Secondary vegetation/Cropland",
                      "Secondary vegetation/Urban",
                      "Plantation forest/Pasture",
                      "Plantation forest/Cropland",
                      "Plantation forest/Urban",
                      "Pasture/Cropland",
                      "Pasture/Urban",
                      "Cropland/Urban")

pairwise_diff_plot_sp$Land_use <- factor(pairwise_diff_plot_sp$Land_use,
                                         levels = priority_combos)


pairwise_diff_plot_sp$Trophic_Niche <- factor(pairwise_diff_plot_sp$Trophic_Niche, levels = c("In", "Fr","Gr", "Ne", "Hb.T", "Hb.A", "Vt", "Aq.p","Om" ))

pairwise_diff_plot_sp <- pairwise_diff_plot_sp %>% dplyr::mutate(Foraging_Niche = ifelse(Foraging_Niche == "For.Gen", "For.Gen","Spec"))




for(r in realms){
  
  data <- pairwise_diff_plot_sp%>% dplyr::filter(Realm == r)
  
  
  
  barwidth <- 0.45
  
  for(d_g in c("In","Gr","Fr","Ne","Om")){
    
    p_d <- data %>% filter(Trophic_Niche == d_g) %>% group_by(Trophic_Niche, Land_use, sp_differences,vol, Foraging_Niche) %>% dplyr::summarise(count = n(), .groups = "drop") %>% dplyr::mutate(relative_sp_occupancy = count/vol) 
    
    if(nrow(p_d) == 0){
      next()
    }
    
    priority_order <- priority_combos[priority_combos %in% unique(p_d$Land_use)]
    
    p_d$Land_use <- factor(p_d$Land_use, levels = priority_order)
    
    
    df_1 <- p_d %>% dplyr::filter(sp_differences == "gain")
    df_2 <- p_d %>% dplyr::filter(sp_differences == "lost")   
    
    
    hole_sp_plot <- ggplot() +
      geom_bar(data = df_1, aes(x = as.numeric(Land_use), y = count, fill = factor(ifelse(Foraging_Niche == "For.Gen","Gain Generalist","Gain Specialist"), levels = c("Gain Specialist","Gain Generalist")) ), stat = "identity", width = barwidth) +
      geom_bar(data = df_2, aes(x = as.numeric(Land_use) + barwidth , y = count, fill = factor(ifelse(Foraging_Niche == "For.Gen","Lost Generalist","Lost Specialist"),levels = c("Lost Specialist","Lost Generalist"))), stat = "identity", width = barwidth) +
      facet_grid(. ~Trophic_Niche) +
      scale_fill_manual(name = "Foraging_Niche", values=c("turquoise4","tomato3","seagreen4", "orange2")) +
      ggtitle(paste(r,sep = "_")) +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90)) +
      scale_x_continuous(name = "Station",
                         limits = c(0.5, length(levels(p_d$Land_use))+1),
                         breaks = unique(sort(as.numeric(p_d$Land_use))),
                         labels = levels(p_d$Land_use)) 
    
    
    
    
    
    plot(hole_sp_plot)
    
  }
}
