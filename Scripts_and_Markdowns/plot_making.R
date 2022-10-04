rm(list = ls())


require(tidyverse)

source("Functions/TPD_3D_plots.R")

PREDICTS_tpds <- readRDS("Outputs/PREDICTS_sites_tpds.rds")
PREDICTS_randomisations <-
  readRDS("Outputs/randomisations_TPD_morpho.rds")

PREDICTS_full <- readRDS("Outputs/refined_predicts.rds")

PREDICTS <- PREDICTS_full %>%  ## PREDICTS data
  dplyr::distinct(SSBS, Predominant_habitat, Realm, SS,Longitude,Latitude) %>% ## pull out land_use type, Subregion, realm etc
  dplyr::mutate(
    Predominant_habitat = ifelse(
      grepl(
        Predominant_habitat,
        pattern = "secondary",
        ignore.case = TRUE
      ),
      "Secondary vegetation",
      paste(Predominant_habitat)
    ),
    Predominant_habitat = ifelse(
      grepl(
        Predominant_habitat,
        pattern = "Primary forest",
        ignore.case = TRUE
      ),
      "Primary vegetation",
      paste(Predominant_habitat)
    )
  ) %>% dplyr::filter(Predominant_habitat != "Primary non-forest") %>% data.frame()  ## merge all secondary sites together


site_info <-
  PREDICTS %>% dplyr::filter(SSBS %in% c(names(PREDICTS_tpds))) %>%
  distinct(SSBS, Realm, Predominant_habitat,Longitude,Latitude) %>% dplyr::filter(Predominant_habitat != "Cannot decide")


Forage <-
  readRDS("../PREDICTS_Taxonomy/PREDICTS_imputed_BL_traits_forage.rds")

drop_sp <- readRDS("Outputs/assembly_drop_spp.rds")

species_pool <- readRDS("Outputs/predicts_sites_species_pools.rds")

table(site_info$Predominant_habitat, site_info$Realm)


realms <- as.character(unique(site_info$Realm))

land_uses <-
  c(
    "Urban",
    "Cropland",
    "Pasture",
    "Plantation forest",
    "Secondary vegetation",
    "Primary vegetation"
  )


limits <-
  max(
    site_info %>% group_by(Predominant_habitat, Realm) %>% dplyr::summarise(count = n()) %>% pull(count)
  )

r <- "Palearctic"
for (r in realms) {
  data <-
    site_info %>% dplyr::filter(Realm == r) %>% dplyr::group_by(Predominant_habitat) %>% dplyr::summarise(count = n()) %>%
    dplyr::mutate(Predominant_habitat = factor(Predominant_habitat,
                                               levels = land_uses)) %>% data.frame()
  
  rownames(data) <- data$Predominant_habitat
  
  if (any(!(land_uses %in% rownames(data)))) {
    for (lu in land_uses[which(!land_uses %in% rownames(data))]) {
      data[lu, "Predominant_habitat"] <- lu
      data[lu, "count"] <- 0
      
    }
    
    
  }
  
  
  
  data <- data[land_uses,]
  
  colours <-
    c(
      "#718879",
      "#E3D438",
      "#EBF787",
      "springgreen2",
      "olivedrab2",
      "chartreuse4"
    )
  
  plot <-
    ggplot(data = data, aes (x = Predominant_habitat, y = count)) +
    ylim(0, limits) +
    geom_bar(stat = "identity",
             show.legend = FALSE,
             fill = colours) +
    coord_flip() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_blank()
    )
  
  
    png(
      filename = paste("Figures/plots_for_figures/sampling", ".png", sep = r),
      width = 200,
      height = 150,
      units = "mm",
      res = 500
    )
  plot(plot)
  dev.off()
  
}

#########################################################################
#########################################################################
#########################################################################
#### Dietary guilds


guilds <- c("In", "Om", "Gr", "Fr", "Ne")
guild_colours <-
  c("mediumblue", "orchid", "orange2", "red2", "olivedrab")




for (r in realms) {

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
    geom_bar(
      data = random_df,
      aes (x = c(0.3, 0.65, 1, 1.35, 1.70), y = count),
      stat = "identity",
      show.legend = FALSE,
      fill = guild_colours,
      colour = c("black", "black", "black", "black", "black"),
      width = 0.3
    ) +
    geom_bar(
      data = obvs_df,
      aes (x = c(2.25, 2.60, 2.95, 3.3, 3.65), y = count),
      stat = "identity",
      show.legend = FALSE,
      fill = guild_colours,
      colour = c("black", "black", "black", "black", "black"),
      width = 0.3
    ) +
    geom_vline(xintercept = 1.95, lwd = 2) +
    xlim(0, 5) +
    ylim(0, 1261) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_blank()
    )
  
  plot(plot)
  
  png(
    filename = paste("Figures/plots_for_figures/species_sampling", ".png", sep = r),
    width = 200,
    height = 150,
    units = "mm",
    res = 500
  )
  plot(plot)
  dev.off()
  
}


##############################################################################################
##############################################################################################
##############################################################################################


obvs_species <-
  PREDICTS_full %>% dplyr::filter(
    !(Birdlife_Name %in% drop_sp),
    SSBS %in% site_info$SSBS,
    Effort_Corrected_Measurement > 0
  ) %>%
  distinct(Birdlife_Name) %>% pull()

obvs_df <- data.frame(Birdlife_Name = obvs_species) %>%
  dplyr::left_join(Forage[, c("Birdlife_Name", "Trophic_Niche")], by = "Birdlife_Name") %>%
  dplyr::group_by(Trophic_Niche) %>% dplyr::summarise(count = n()) %>% dplyr::filter(Trophic_Niche %in% guilds) %>%
  data.frame()

rownames(obvs_df) <- obvs_df$Trophic_Niche

if (any(!(guilds %in% rownames(obvs_df)))) {
  for (g in guilds[which(!guilds %in% rownames(obvs_df))]) {
    obvs_df[g, "Trophic_Niche"] <- g
    obvs_df[g, "count"] <- 0
    
  }
  
  
}

obvs_df$Trophic_Niche <-
  factor(obvs_df$Trophic_Niche, levels = guilds)



obvs_df <- obvs_df[guilds,]


random_sp <- c()
for (rs in site_info$SSBS) {
  random_sp <- unique(c(random_sp, species_pool[[rs]][["0.9"]]))
}

random_df <- data.frame(Birdlife_Name = random_sp) %>%
  dplyr::left_join(Forage[, c("Birdlife_Name", "Trophic_Niche")], by = "Birdlife_Name") %>%
  dplyr::group_by(Trophic_Niche) %>% dplyr::summarise(count = n()) %>% dplyr::filter(Trophic_Niche %in% guilds) %>%
  data.frame()

rownames(random_df) <- random_df$Trophic_Niche

if (any(!(guilds %in% rownames(obvs_df)))) {
  for (g in guilds[which(!guilds %in% rownames(obvs_df))]) {
    obvs_df[g, "Trophic_Niche"] <- g
    obvs_df[g, "count"] <- 0
    
  }
  
  
}


random_df$Trophic_Niche <-
  factor(random_df$Trophic_Niche, levels = guilds)

random_df <- random_df[guilds,]


plot <- ggplot() +
  geom_bar(
    data = random_df,
    aes (x = c(0.3, 0.65, 1, 1.35, 1.70), y = count),
    stat = "identity",
    show.legend = FALSE,
    fill = guild_colours,
    colour = c("black", "black", "black", "black", "black"),
    width = 0.3
  ) +
  geom_bar(
    data = obvs_df,
    aes (x = c(2.25, 2.60, 2.95, 3.3, 3.65), y = count),
    stat = "identity",
    show.legend = FALSE,
    fill = guild_colours,
    colour = c("black", "black", "black", "black", "black"),
    width = 0.3
  ) +
  geom_vline(xintercept = 1.95, lwd = 2) +
  xlim(0, 5) +
  ylim(0, 2972) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank()
  )

plot(plot)

png(
  filename = paste("Figures/plots_for_figures/species_sampling", ".png", sep = "Total"),
  width = 200,
  height = 150,
  units = "mm",
  res = 500
)
plot(plot)
dev.off()


#### going to need to make a map showing the site points 
#### 
#### 
#### 



library(tidyverse) ## for wrangling and general data handling
library(terra) # mapping
library(stars)
library(sp) ## mapping
library(raster) ## mapping
require(sf) ## mapping
require(doParallel) ## parallelisation
require(foreach) ## parallelisation
require(fasterize)


biome_polygons <- st_read("../../Datasets/WWF_Terrestrial_Biomes/official/wwf_terr_ecos.shp")


biome_polygons <- biome_polygons %>% dplyr::filter(REALM %in% c("NT","PA","NA","IM","AT","AA"))


realm_colours <- c("yellowgreen","royalblue4","darkturquoise","orangered4","goldenrod1","darkorchid4")


wm<-map_data("world") %>% filter(region != "Antartica") %>% fortify() 


boundry_plot <-ggplot()+ coord_fixed()+
  geom_map(data =wm, map = wm,
           aes(group = group, map_id= region),
           fill = "darkgrey") 



i <- 1
for(realm in unique(biome_polygons$REALM)){

  realm_poly <- biome_polygons %>% dplyr::filter(REALM == realm)
  realm_poly <- st_combine(realm_poly$geometry)
  
  boundry_plot <- boundry_plot +
    geom_sf(data = st_as_sf(realm_poly), fill = realm_colours[i], colour = realm_colours[i])

  
    i <- i + 1
}


boundry_plot <- boundry_plot +
  geom_point(data = site_info, aes(x = Longitude, y = Latitude), size = 7, alpha = 0.2) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), panel.background = element_blank()) +
  xlim(-155,171) +
  ylim(-50,90)


ggsave(filename = "Figures/plots_for_figures/Rplot_map.png",plot = boundry_plot, device = "png", dpi = 600)


plot(boundry_plot)






