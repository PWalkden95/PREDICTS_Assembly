rm(list = ls())

require(gamm4)
require(mgcv)
require(tidyverse)
require(doParallel)






source("Functions/TPD_3D_plots.R")


PREDICTS_tpds <-
  readRDS("Outputs/PREDICTS_sites_tpds.rds") # morphometric TPDs of observed sites

PREDICTS_randomisations <-
  readRDS("Outputs/randomisations_TPD_morpho.rds") ## morphometric TPD of randomisised sites




PREDICTS_full <- readRDS("Outputs/refined_predicts.rds")

PREDICTS <- PREDICTS_full %>%  ## PREDICTS data
  dplyr::distinct(SSBS, Predominant_habitat, Realm, SS) %>% ## pull out land_use type, Subregion, realm etc
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
        pattern = "primary",
        ignore.case = TRUE
      ),
      "Primary vegetation",
      paste(Predominant_habitat)
    )
  ) %>% data.frame() ## merge all secondary sites together



TPD_land_uses <-
  data.frame(SSBS = names(PREDICTS_tpds)) %>% dplyr::left_join(PREDICTS, by = "SSBS") %>% dplyr::filter(Predominant_habitat != "Cannot decide")
### check a table
table(TPD_land_uses$Predominant_habitat, TPD_land_uses$Realm)


land_uses <-
  c(
    "Primary vegetation",
    "Secondary vegetation",
    "Pasture",
    "Cropland",
    "Plantation forest",
    "Urban"
  )
realms <-
  c("Neotropic",
    "Afrotropic",
    "Palearctic",
    "Nearctic",
    "Indo-Malay",
    "Australasia")


land_use_colours <-
  data.frame(
    land_use = land_uses,
    colours = c(
      "chartreuse4",
      "olivedrab2",
      "#EBF787",
      "#E3D438",
      "springgreen2",
      "#718879"
    )
  )
rownames(land_use_colours) <- land_uses


################################
################################
####### Again at the site level there isn't anything really that clear let's try again at the Realm level




realm_level_hole_frame <- c()

for (land_use in land_uses) {
  for (realm in realms) {
    sites <-
      TPD_land_uses %>% dplyr::filter(Predominant_habitat == land_use, Realm == realm) %>% dplyr::distinct(SSBS) %>% pull()
    
    if (is_empty(sites)) {
      next()
    }
    
    LU_holes <-
      TPD_holes(
        data = PREDICTS_tpds,
        randata = PREDICTS_randomisations,
        sites = sites,
        threshold = 0.95,
        minimum_points = 10
      )
    
    hole_data <- LU_holes$internal$metrics_frame
    
    hole_data$land_use <- land_use
    hole_data$realm <- realm
    
    
    realm_level_hole_frame <-
      rbind(realm_level_hole_frame, hole_data)
  }
}

write_rds(file = "Outputs/realm_level_hole_dbscan_ninety_five.rds", realm_level_hole_frame)

realm_level_hole_frame <-
  readRDS("Outputs/realm_level_hole_dbscan.rds")


realm_level_hole_frame$land_use <-
  factor(
    realm_level_hole_frame$land_use,
    levels = c(
      "Primary vegetation",
      "Secondary vegetation",
      "Plantation forest",
      "Pasture",
      "Cropland",
      "Urban"
    )
  )


realm_level_hole_frame_2 <-
  realm_level_hole_frame[-which(
    realm_level_hole_frame$land_use == "Urban" &
      realm_level_hole_frame$realm == "Neotropic" |
      realm_level_hole_frame$land_use == "Secondary vegetation" &
      realm_level_hole_frame$realm == "Nearctic" |
      realm_level_hole_frame$land_use == "Plantation forest" &
      realm_level_hole_frame$realm == "Nearctic" |
      realm_level_hole_frame$land_use == "Urban" &
      realm_level_hole_frame$realm == "Australasia" |
      realm_level_hole_frame$land_use == "Plantation forest" &
      realm_level_hole_frame$realm == "Australasia" |
      realm_level_hole_frame$land_use == "Cropland" &
      realm_level_hole_frame$realm == "Australasia" |
      realm_level_hole_frame$land_use == "Pasture" &
      realm_level_hole_frame$realm == "Indo-Malay" |
      realm_level_hole_frame$land_use == "Urban" &
      realm_level_hole_frame$realm == "Indo-Malay"
  ), ]

realm_level_hole_frame_2$hypervolume_occupancy <-
  scale(
    log(realm_level_hole_frame_2$hypervolume_occupancy),
    center = TRUE,
    scale = TRUE
  )[, 1]

realm_level_hole_frame_2$total_hole_volume <-
  scale(
    sqrt(realm_level_hole_frame_2$total_hole_volume),
    center = TRUE,
    scale = TRUE
  )[, 1]




realm_level_occupancy_data <-
  realm_level_hole_frame_2 %>% dplyr::group_by(land_use, realm) %>%
  dplyr::summarise(hypervolume_occupancy = hypervolume_occupancy,
                   total_volume = max(total_hole_volume),
                   mean_total_volume = mean(total_hole_volume))



realm_level_radius_data <- realm_level_hole_frame_2 %>% group_by(land_use,radius) %>%
  dplyr::summarise(mean_hole_volume = mean(total_hole_volume))



radius_volume <-
  ggplot(data = realm_level_hole_frame_2, aes(x = radius, y = total_hole_volume, colour = land_use)) +
  geom_point()

plot(radius_volume)


mean_volume_radius <- 
  ggplot(data = realm_level_radius_data, aes(x = radius, y = mean_hole_volume, colour = land_use)) +
  geom_point()


plot(mean_volume_radius)  
  
occupancy_volume <-
  ggplot(data = realm_level_occupancy_data,
         aes(x = hypervolume_occupancy, y = total_volume, colour = land_use)) +
  geom_point()

plot(occupancy_volume)


radius_number_of_holes <-
  ggplot(data = realm_level_hole_frame, aes(x = radius, y = number_of_holes, colour = land_use)) +
  geom_point()

plot(radius_number_of_holes)

occupancy_holes <-
  ggplot(data = realm_level_hole_frame,
         aes(x = hypervolume_occupancy, y = number_of_holes, colour = land_use)) +
  geom_point()

plot(occupancy_holes)


hist(scale(
  sqrt(realm_level_hole_frame$total_hole_volume),
  center = TRUE,
  scale = TRUE
))


realm_level_hole_frame_2 <-
  realm_level_hole_frame[-which(
    realm_level_hole_frame$land_use == "Urban" &
      realm_level_hole_frame$realm == "Neotropic" |
      realm_level_hole_frame$land_use == "Secondary vegetation" &
      realm_level_hole_frame$realm == "Nearctic" |
      realm_level_hole_frame$land_use == "Plantation forest" &
      realm_level_hole_frame$realm == "Nearctic" |
      realm_level_hole_frame$land_use == "Urban" &
      realm_level_hole_frame$realm == "Australasia" |
      realm_level_hole_frame$land_use == "Plantation forest" &
      realm_level_hole_frame$realm == "Australasia" |
      realm_level_hole_frame$land_use == "Cropland" &
      realm_level_hole_frame$realm == "Australasia" |
      realm_level_hole_frame$land_use == "Pasture" &
      realm_level_hole_frame$realm == "Indo-Malay" |
      realm_level_hole_frame$land_use == "Urban" &
      realm_level_hole_frame$realm == "Indo-Malay"
  ), ]

realm_level_hole_frame_2$hypervolume_occupancy <-
  scale(
    log(realm_level_hole_frame_2$hypervolume_occupancy),
    center = TRUE,
    scale = TRUE
  )[, 1]

realm_level_hole_frame_2$total_hole_volume <-
  scale(
    sqrt(realm_level_hole_frame_2$total_hole_volume),
    center = TRUE,
    scale = TRUE
  )[, 1]



realm_gam <- mgcv::gam(
  total_hole_volume ~ s(radius) +
    s(radius, by = land_use) +
    hypervolume_occupancy +
    #te(hypervolume_occupancy,radius)+
    hypervolume_occupancy:land_use  +
    land_use,
  data = realm_level_hole_frame_2,
  method = "REML"
)



summary(realm_gam)
anova(realm_gam)

radius_sequence <- seq(min(realm_level_hole_frame_2$radius), max(realm_level_hole_frame_2$radius), 0.01)




hole_occupancy <- seq(min(realm_level_hole_frame_2$hypervolume_occupancy),
    max(realm_level_hole_frame_2$hypervolume_occupancy),length.out = 50)


k <- 1

for(occupancy in hole_occupancy){


gamm_plot_radius <- c()



for (land_use in land_uses) {
  land_use_data <-
    realm_level_hole_frame_2[realm_level_hole_frame_2$land_use == land_use, "radius"]
  
  #hyper_size <- range(unique(realm_level_hole_frame_2[realm_level_hole_frame_2$land_use == land_use, "hypervolume_occupancy"]))
  #hyper_size <- seq(hyper_size[1],hyper_size[2], length.out = 100)
  
  
  for (i in land_use_data) {
  #    mean_volume <- median(unique(realm_level_hole_frame_2[realm_level_hole_frame_2$land_use == land_use,"hypervolume_occupancy"]))
  mean_volume <-
      median(unique(realm_level_hole_frame_2$hypervolume_occupancy))
    
    
    data <- data.frame(
      radius = i,
      land_use = land_use,
      predict(
        realm_gam,
        newdata = data.frame(
          radius = i,
          land_use = land_use,
          hypervolume_occupancy = occupancy
        ),
        se.fit = TRUE
      )
    )
    
    gamm_plot_radius <- rbind(gamm_plot_radius, data)
  }
}


gamm_plot_radius$land_use <-
  factor(
    gamm_plot_radius$land_use,
    levels = c(
      "Primary vegetation",
      "Secondary vegetation",
      "Pasture",
      "Cropland",
      "Plantation forest",
      "Urban"
    )
  )


gamm_plot_radius$upper <-
  gamm_plot_radius$fit + gamm_plot_radius$se.fit
gamm_plot_radius$lower <-
  gamm_plot_radius$fit - gamm_plot_radius$se.fit


occupancy_land_uses_high <- as.character(unique(realm_level_hole_frame_2[realm_level_hole_frame_2$hypervolume_occupancy >= occupancy,"land_use"]))
occupancy_land_uses_low <- as.character(unique(realm_level_hole_frame_2[realm_level_hole_frame_2$hypervolume_occupancy <= occupancy,"land_use"]))


gamm_plot_radius <- gamm_plot_radius %>% dplyr::filter(land_use %in% occupancy_land_uses_high & land_use %in% occupancy_land_uses_low)



test <-
  ggplot(data = gamm_plot_radius, aes(x = radius, y = fit, group = land_use)) +
  #geom_point(data = realm_level_radius_data, aes( x = radius, y = mean_hole_volume, colour = land_use), alpha = 0.7, size = 3) +
  geom_line(size = 1.5, show.legend = FALSE,aes(colour = land_use)) +
  geom_ribbon(aes(
    ymin = lower,
    ymax = upper,
    fill = land_use
  ),
  alpha = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = land_use_colours[as.character(unique(gamm_plot_radius$land_use)), "colours"]) +
  scale_fill_manual(values = land_use_colours[as.character(unique(gamm_plot_radius$land_use)), "colours"]) +
  theme(panel.background = element_rect(fill = 'grey', color = 'grey'),
panel.grid.major = element_line(color = 'grey'),
panel.grid.minor = element_line(color = 'grey')) +
  xlab("Radius") +
  ylab("Total Hole Volume") +
  ylim(-4,3) +
  xlim(min(radius_sequence), max(radius_sequence)) +
  ggtitle(paste(round(occupancy,digits = 2)))


plot(test)


ggsave(filename = paste("Figures/Animations/radius_hole_volume/radius_hole",k,".png"), plot = test, device = "png", height = 10, width = 15, dpi = 300)

k <- k +1
}



###################################################
###################################################

k <- 1

for(rad in radius_sequence) {


gamm_plot_hypervolume_occupancy <- c()


for (land_use in land_uses) {
  #land_use_data <- realm_level_hole_frame_2[realm_level_hole_frame_2$land_use == land_use,"radius"]
  
  hyper_size <-
    range(unique(realm_level_hole_frame_2[realm_level_hole_frame_2$land_use == land_use, "hypervolume_occupancy"]))
  hyper_size <- seq(hyper_size[1], hyper_size[2], length.out = 100)
  
  
  for (i in hyper_size) {
    
    #mean_radius <-
      round(median(unique(realm_level_hole_frame_2[realm_level_hole_frame_2$land_use == land_use, "radius"])), digits = 2)
    #  mean_volume <- median(unique(realm_level_hole_frame_2$hypervolume_occupancy))
    #mean_radius <- round(median(unique(realm_level_hole_frame_2$radius)),digits = 2)
    
    data <- data.frame(
      radius = i,
      land_use = land_use,
      predict(
        realm_gam,
        newdata = data.frame(
          radius = rad,
          land_use = land_use,
          hypervolume_occupancy = i
        ),
        se.fit = TRUE
      )
    )
    
    gamm_plot_hypervolume_occupancy <-
      rbind(gamm_plot_hypervolume_occupancy, data)
  }
}


gamm_plot_hypervolume_occupancy$land_use <-
  factor(
    gamm_plot_hypervolume_occupancy$land_use,
    levels = c(
      "Primary vegetation",
      "Secondary vegetation",
      "Pasture",
      "Cropland",
      "Plantation forest",
      "Urban"
    )
  )

gamm_plot_hypervolume_occupancy$upper <-
  gamm_plot_hypervolume_occupancy$fit + gamm_plot_hypervolume_occupancy$se.fit
gamm_plot_hypervolume_occupancy$lower <-
  gamm_plot_hypervolume_occupancy$fit - gamm_plot_hypervolume_occupancy$se.fit


radius_points <- c()

for (land_use in land_uses) {
  mean_radius <-
    round(median(unique(realm_level_hole_frame_2[realm_level_hole_frame_2$land_use == land_use, "radius"])), digits = 2)
  
  radius_data <-
    realm_level_hole_frame_2 %>% dplyr::mutate(radius = as.character(radius)) %>%
    dplyr::filter(radius == as.character(rad))
  
  radius_data <- radius_data[radius_data$land_use == land_use,]
  
  radius_points <- rbind(radius_points, radius_data)
  
}


gamm_plot_hypervolume_occupancy <- gamm_plot_hypervolume_occupancy %>% dplyr::filter(land_use %in% as.character(unique(radius_points$land_use)))




test <-
  ggplot(data = gamm_plot_hypervolume_occupancy, aes(x = radius, y = fit, group = land_use)) +
  geom_line(size = 1.5, show.legend = FALSE,aes(colour = land_use)) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper,
      fill = land_use
    ),
    alpha = 0.6,
    show.legend = FALSE
  ) +
    geom_point(
    data = radius_points,
    aes(x = hypervolume_occupancy, y = total_hole_volume, colour = land_use)
    ,
    size = 8,
    alpha = 0.6,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = land_use_colours[as.character(unique(gamm_plot_hypervolume_occupancy$land_use)), "colours"]) +
  scale_fill_manual(values = land_use_colours[as.character(unique(gamm_plot_hypervolume_occupancy$land_use)), "colours"]) +
  ylab("Total Hole Volume") +
  xlab("Hypervolume occupancy") +
  theme(panel.background = element_rect(fill = 'grey', color = 'grey'),
        panel.grid.major = element_line(color = 'grey'),
        panel.grid.minor = element_line(color = 'grey')) +
  ylim(-4,3)+
  xlim(-3,2)


plot(test)

ggsave(filename = paste("Figures/Animations/occupancy_hole/occupancy_hole",k,".png"), plot = test, device = "png", height = 10, width = 15, dpi = 300)

k <- k +1

}




#######################################################
#######################################################
#######################################################

realm_gam_holes <- mgcv::gam(
  number_of_holes ~ s(radius) +
    s(radius, by = land_use) +
    hypervolume_occupancy +
    hypervolume_occupancy:land_use +
    #te(hypervolume_occupancy,radius)+
    # s(hypervolume_occupancy, by = land_use)  +
    land_use,
  data = realm_level_hole_frame_2,
  family = "poisson",
  method = "REML"
)


summary(realm_gam_holes)
anova(realm_gam_holes)

gamm_plot_holes <- c()


for (land_use in land_uses) {
  land_use_data <-
    realm_level_hole_frame_2[realm_level_hole_frame_2$land_use == land_use, "radius"]
  
  #hyper_size <- range(unique(realm_level_hole_frame_2[realm_level_hole_frame_2$land_use == land_use, "hypervolume_occupancy"]))
  #hyper_size <- seq(hyper_size[1],hyper_size[2], length.out = 100)
  
  
  for (i in land_use_data) {
    #mean_volume <- median(unique(realm_level_hole_frame_2[realm_level_hole_frame_2$land_use == land_use,"hypervolume_occupancy"]))
    mean_volume <-
      mean(unique(realm_level_hole_frame_2$hypervolume_occupancy))
    
    
    data <- data.frame(
      radius = i,
      land_use = land_use,
      predict(
        realm_gam_holes,
        newdata = data.frame(
          radius = i,
          land_use = land_use,
          hypervolume_occupancy = mean_volume
        ),
        se.fit = TRUE
      )
    )
    
    gamm_plot_holes <- rbind(gamm_plot_holes, data)
  }
}


gamm_plot_holes$land_use <-
  factor(
    gamm_plot_holes$land_use,
    levels = c(
      "Primary vegetation",
      "Secondary vegetation",
      "Plantation forest",
      "Pasture",
      "Cropland",
      "Urban"
    )
  )

gamm_plot_holes$upper <-
  gamm_plot_holes$fit + gamm_plot_holes$se.fit






test <-
  ggplot(data = gamm_plot_holes, aes(x = radius, y = fit, colour = land_use)) +
  geom_point() +
  #  geom_point(data = radius_one, aes(x = hypervolume_occupancy, y = number_of_holes)
  #            ,size = 7, alpha = 0.5) +
  ylab("number of holes") +
  xlab("Radius")


plot(test)

######################################################
######################################################
######################################################



gamm_plot_holes_hypervolume_occupancy <- c()


for (land_use in land_uses) {
  #land_use_data <- realm_level_hole_frame_2[realm_level_hole_frame_2$land_use == land_use,"radius"]
  
  hyper_size <-
    range(unique(realm_level_hole_frame_2[realm_level_hole_frame_2$land_use == land_use, "hypervolume_occupancy"]))
  hyper_size <- seq(hyper_size[1], hyper_size[2], length.out = 100)
  
  
  for (i in hyper_size) {
    #mean_volume <- mean(unique(realm_level_hole_frame_2[realm_level_hole_frame_2$land_use == land_use,"hypervolume_occupancy"]))
    mean_radius <-
      round(median(unique(realm_level_hole_frame_2$radius)), digits = 2)
    
    
    data <- data.frame(
      radius = i,
      land_use = land_use,
      mean_proportion = predict(
        realm_gam_holes,
        newdata = data.frame(
          radius = mean_radius,
          land_use = land_use,
          hypervolume_occupancy = i
        )
      )
    )
    
    gamm_plot_holes_hypervolume_occupancy <-
      rbind(gamm_plot_holes_hypervolume_occupancy, data)
  }
}



gamm_plot_holes_hypervolume_occupancy$land_use <-
  factor(
    gamm_plot_holes_hypervolume_occupancy$land_use,
    levels = c(
      "Primary vegetation",
      "Secondary vegetation",
      "Plantation forest",
      "Pasture",
      "Cropland",
      "Urban"
    )
  )





radius_one <-
  realm_level_hole_frame_2 %>% dplyr::mutate(radius = as.character(radius)) %>%
  dplyr::filter(radius == "0.75")

plot(radius_one$number_of_holes ~ radius_one$land_use)


gamm_plot_2 <-
  gamm_plot_holes_hypervolume_occupancy %>% dplyr::filter(land_use %in% unique(radius_one$land_use))

test <-
  ggplot(data = gamm_plot_2, aes(x = radius, y = mean_proportion, colour = land_use)) +
  geom_line(size = 1.5,
            show.legend = FALSE,
            fill =) +
  scale_colour_manual(values = land_use_colours[order(land_uses), "colours"]) +
  #geom_point(data = radius_one, aes(x = hypervolume_occupancy, y = total_hole_volume)
  #          ,size = 7, alpha = 0.5) +
  ylab("Number of Holes") +
  xlab("Hypervolume Occupancy") +
  theme(
    panel.background = element_rect(fill = 'lightgrey', color = 'lightgrey'),
    panel.grid.major = element_line(color = 'lightgrey'),
    panel.grid.minor = element_line(color = 'lightgrey')
  )


plot(test)
