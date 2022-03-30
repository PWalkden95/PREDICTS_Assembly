rm(list = ls())


require(tidyverse)

source("Functions/TPD_3D_plots.R")

PREDICTS_tpds <- readRDS("Outputs/PREDICTS_sites_tpds.rds")
PREDICTS_randomisations <- readRDS("Outputs/randomisations_TPD_morpho.rds")

PREDICTS_full <- readRDS("Outputs/refined_predicts.rds")

PREDICTS <- PREDICTS_full %>%  ## PREDICTS data 
  dplyr::distinct(SSBS, Predominant_habitat, Realm, SS) %>% ## pull out land_use type, Subregion, realm etc
  dplyr::mutate(Predominant_habitat = ifelse(grepl(Predominant_habitat, pattern = "secondary", ignore.case = TRUE), "Secondary vegetation",paste(Predominant_habitat)), 
                Predominant_habitat = ifelse(grepl(Predominant_habitat, pattern = "primary", ignore.case = TRUE), "Primary vegetation",paste(Predominant_habitat))) %>% data.frame() ## merge all secondary sites together


site_info <- PREDICTS %>% dplyr::filter(SSBS %in% c(names(PREDICTS_tpds))) %>% 
  distinct(SSBS,Realm,Predominant_habitat) %>% dplyr:: filter(Predominant_habitat != "Cannot decide")


Forage <- readRDS("../PREDICTS_Taxonomy/PREDICTS_imputed_BL_traits_forage.rds")

drop_sp <- readRDS("Outputs/assembly_drop_spp.rds")

species_pool <- readRDS("Outputs/predicts_sites_species_pools.rds")

table(site_info$Predominant_habitat, site_info$Realm)


realms <- as.character(unique(site_info$Realm))

land_uses <- c("Urban","Cropland","Pasture","Plantation forest","Secondary vegetation","Primary vegetation")


limits <- max(site_info %>% group_by(Predominant_habitat,Realm) %>% dplyr::summarise(count = n()) %>% pull(count))


for(r in realms){
  data <- site_info %>% dplyr::filter(Realm == r) %>% dplyr::group_by(Predominant_habitat) %>% dplyr::summarise(count = n()) %>%
    dplyr::mutate(Predominant_habitat = factor(Predominant_habitat,
                                                                                   levels = land_uses)) %>% data.frame()
  
  rownames(data) <- data$Predominant_habitat
  
  if(any(!(land_uses %in% rownames(data)))){
    
    for(lu in land_uses[which(!land_uses %in% rownames(data))]){
    
      data[lu,"Predominant_habitat"] <- lu
      data[lu,"count"] <- 0
      
    }
    
    
  }
  
  
  
  data <- data[land_uses,]
  
    colours <- c("ivory4","gold1","khaki","lightgreen","olivedrab3","green4")
    
  plot <- ggplot(data = data, aes (x = Predominant_habitat,y=count)) +
    ylim(0,limits) +
    geom_bar(stat = "identity", show.legend = FALSE, fill = colours) +
    coord_flip() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background = element_blank())

  
    png(filename = paste("C:/Users/patri/Desktop/sampling",".png", sep = r), width = 200, height = 150, units = "mm", res = 500)
    plot(plot)
    dev.off()
  
}

#########################################################################
#########################################################################
#########################################################################
#### Dietary guilds


guilds <- c("In","Om","Gr","Fr","Ne")
guild_colours <- c("mediumblue","white","orange2","red2","olivedrab")


r <- "Palearctic"

for(r in realms){
    realm_sites <- site_info %>% dplyr::filter(Realm == r) %>% distinct(SSBS) %>% pull() %>% as.character()
  
    obvs_species <- PREDICTS_full %>% dplyr::filter(!(Birdlife_Name %in% drop_sp), SSBS %in% realm_sites, Effort_Corrected_Measurement > 0) %>%
      distinct(Birdlife_Name) %>% pull()
    
    obvs_df <- data.frame(Birdlife_Name = obvs_species) %>% 
      dplyr::left_join(Forage[,c("Birdlife_Name", "Trophic_Niche")], by = "Birdlife_Name") %>%
      dplyr::group_by(Trophic_Niche) %>% dplyr::summarise(count = n()) %>% dplyr::filter(Trophic_Niche %in% guilds) %>%
      data.frame()
    
    rownames(obvs_df) <- obvs_df$Trophic_Niche
    
    if(any(!(guilds %in% rownames(obvs_df)))){
      
      for(g in guilds[which(!guilds %in% rownames(obvs_df))]){
        
        obvs_df[g,"Trophic_Niche"] <- g
        obvs_df[g,"count"] <- 0
        
      }
      
      
    }
    
    obvs_df$Trophic_Niche <- factor(obvs_df$Trophic_Niche, levels = guilds)
    
    
    
    obvs_df <- obvs_df[guilds,]
    
    
    random_sp <- c()
    for(rs in realm_sites){
      random_sp <- unique(c(random_sp,species_pool[[rs]][["0.9"]]))
    }
    
    random_df <- data.frame(Birdlife_Name = random_sp) %>% 
      dplyr::left_join(Forage[,c("Birdlife_Name", "Trophic_Niche")], by = "Birdlife_Name") %>%
      dplyr::group_by(Trophic_Niche) %>% dplyr::summarise(count = n()) %>% dplyr::filter(Trophic_Niche %in% guilds) %>%
      data.frame()
    
    rownames(random_df) <- random_df$Trophic_Niche
    
    if(any(!(guilds %in% rownames(obvs_df)))){
      
      for(g in guilds[which(!guilds %in% rownames(obvs_df))]){
        
        obvs_df[g,"Trophic_Niche"] <- g
        obvs_df[g,"count"] <- 0
        
      }
      
      
    }
    
    
    random_df$Trophic_Niche <- factor(random_df$Trophic_Niche, levels = guilds)
    
    random_df <- random_df[guilds,]
    
    
    plot <- ggplot() +
      geom_bar(data = random_df, aes (x = c(0.3,0.65,1,1.35,1.70),y=count), stat = "identity", show.legend = FALSE, fill = guild_colours, colour = c("white","black","white","white","white"), width = 0.3) +
      geom_bar(data = obvs_df, aes (x = c(2.25,2.60,2.95,3.3,3.65),y=count), stat = "identity", show.legend = FALSE, fill = guild_colours, colour = c("white","black","white","white","white"), width = 0.3) +
      geom_vline(xintercept = 1.95, lwd = 2) + 
      xlim(0,5) +
      ylim(0,1261) +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.background = element_blank())
    
    plot(plot)
    
    png(filename = paste("C:/Users/patri/Desktop/species_sampling",".png", sep = r), width = 200, height = 150, units = "mm", res = 500)
    plot(plot)
    dev.off()

    }


##############################################################################################
##############################################################################################
##############################################################################################


obvs_species <- PREDICTS_full %>% dplyr::filter(!(Birdlife_Name %in% drop_sp), SSBS %in% site_info$SSBS, Effort_Corrected_Measurement > 0) %>%
  distinct(Birdlife_Name) %>% pull()

obvs_df <- data.frame(Birdlife_Name = obvs_species) %>% 
  dplyr::left_join(Forage[,c("Birdlife_Name", "Trophic_Niche")], by = "Birdlife_Name") %>%
  dplyr::group_by(Trophic_Niche) %>% dplyr::summarise(count = n()) %>% dplyr::filter(Trophic_Niche %in% guilds) %>%
  data.frame()

rownames(obvs_df) <- obvs_df$Trophic_Niche

if(any(!(guilds %in% rownames(obvs_df)))){
  
  for(g in guilds[which(!guilds %in% rownames(obvs_df))]){
    
    obvs_df[g,"Trophic_Niche"] <- g
    obvs_df[g,"count"] <- 0
    
  }
  
  
}

obvs_df$Trophic_Niche <- factor(obvs_df$Trophic_Niche, levels = guilds)



obvs_df <- obvs_df[guilds,]


random_sp <- c()
for(rs in site_info$SSBS){
  random_sp <- unique(c(random_sp,species_pool[[rs]][["0.9"]]))
}

random_df <- data.frame(Birdlife_Name = random_sp) %>% 
  dplyr::left_join(Forage[,c("Birdlife_Name", "Trophic_Niche")], by = "Birdlife_Name") %>%
  dplyr::group_by(Trophic_Niche) %>% dplyr::summarise(count = n()) %>% dplyr::filter(Trophic_Niche %in% guilds) %>%
  data.frame()

rownames(random_df) <- random_df$Trophic_Niche

if(any(!(guilds %in% rownames(obvs_df)))){
  
  for(g in guilds[which(!guilds %in% rownames(obvs_df))]){
    
    obvs_df[g,"Trophic_Niche"] <- g
    obvs_df[g,"count"] <- 0
    
  }
  
  
}


random_df$Trophic_Niche <- factor(random_df$Trophic_Niche, levels = guilds)

random_df <- random_df[guilds,]


plot <- ggplot() +
  geom_bar(data = random_df, aes (x = c(0.3,0.65,1,1.35,1.70),y=count), stat = "identity", show.legend = FALSE, fill = guild_colours, colour = c("white","black","white","white","white"), width = 0.3) +
  geom_bar(data = obvs_df, aes (x = c(2.25,2.60,2.95,3.3,3.65),y=count), stat = "identity", show.legend = FALSE, fill = guild_colours, colour = c("white","black","white","white","white"), width = 0.3) +
  geom_vline(xintercept = 1.95, lwd = 2) + 
  xlim(0,5) +
  ylim(0,2972) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank())

plot(plot)

png(filename = paste("C:/Users/patri/Desktop/species_sampling",".png", sep = "Total"), width = 200, height = 150, units = "mm", res = 500)
plot(plot)
dev.off()




