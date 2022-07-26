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
    )
  ) %>% data.frame() ## merge all secondary sites together



TPD_land_uses <-
  data.frame(SSBS = names(PREDICTS_tpds)) %>% dplyr::left_join(PREDICTS, by = "SSBS") %>% dplyr::filter(Predominant_habitat != "Cannot decide")
### check a table
table(TPD_land_uses$Predominant_habitat, TPD_land_uses$Realm)


land_uses <-
  c(
    "Primary forest",
    "Primary non-forest",
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
      "chartreuse3",
      "olivedrab2",
      "#EBF787",
      "#E3D438",
      "springgreen2",
      "#718879"
    )
  )
rownames(land_use_colours) <- land_uses


registerDoParallel(cores = 12)


hole_frame <- foreach(sites = unique(TPD_land_uses$SSBS),.combine = "rbind",.inorder = FALSE,.packages = c("tidyverse","magrittr","geometry")) %dopar% {
  
  
  
  LU_holes <- TPD_holes(data = PREDICTS_tpds, randata = PREDICTS_randomisations, sites = sites, threshold = 1, minimum_points = 10)
  
  h_dat <- LU_holes$internal$metrics_frame
  
  if(is_null(h_dat)){
    hole_frame <- c()
  } else {
    
    h_dat$land_use <- TPD_land_uses %>% dplyr::filter(SSBS == sites) %>% dplyr::distinct(Predominant_habitat) %>% pull()
    h_dat$SSBS <- sites 
    h_dat$SS <- TPD_land_uses %>% dplyr::filter(SSBS == sites) %>% dplyr::distinct(SS) %>% pull()
    
    
    
    hole_frame <- h_dat
  }
  
  return(hole_frame)
}

registerDoSEQ()
closeAllConnections()

length(unique(hole_frame$SSBS)) ## 20 sites that couldn't compute the holes metric for -- usually that there are too few cells
##  occupied


## write


hole_frame <- hole_frame %>% dplyr::filter(!is.na(hypervolume_occupancy))
hole_frame <- hole_frame %>% dplyr::mutate(SS = factor(SS),
                                           SSBS = factor(SSBS),
                                           land_use = factor(land_use, levels = c("Primary forest",
                                                                                  "Primary non-forest",
                                                                                  "Secondary vegetation",
                                                                                  "Plantation forest",
                                                                                  "Pasture",
                                                                                  "Cropland",
                                                                                  "Urban")))


write_rds(x = hole_frame, file = "Outputs/site_holes_dbscan.rds")

hole_frame <- readRDS("Outputs/site_holes_dbscan.rds")




hole_frame$land_use <- factor(hole_frame$land_use, levels = c("Primary forest","Primary non-forest", "Secondary vegetation", "Plantation forest",
                                                              "Pasture", "Cropland", "Urban"))


try(memory.limit(120000000000))

gam <- gamm4::gamm4(total_hole_volume ~ #s(radius, by = land_use) + 
                      s(radius) + 
                      hypervolume_occupancy +
                      hypervolume_occupancy:land_use +
                      land_use,
                    random = ~(1|SS), data = hole_frame)


summary(gam$gam)




land_use_plot <- c()

median_volume <- median(unique(hole_frame$hypervolume_occupancy))

median_radius <- median(unique(hole_frame$radius))



for(land_use in land_uses){
  
  
  
  data <- data.frame(land_use = land_use, predict(gam$gam, newdata = data.frame(radius = median_radius,
                                                                                land_use = land_use,
                                                                                hypervolume_occupancy = median_volume), se.fit = TRUE, re.form = NA))
  
  
  land_use_plot <- rbind(land_use_plot,data)
  
}


land_use_plot$fit <- land_use_plot$fit - land_use_plot[land_use_plot$land_use == "Primary forest", "fit"]
land_use_plot$lower <- land_use_plot$fit - land_use_plot$se.fit
land_use_plot$upper <- land_use_plot$fit + land_use_plot$se.fit

land_use_plot$land_use <- factor(land_use_plot$land_use, levels = c("Primary forest","Primary non-forest", "Secondary vegetation", "Plantation forest", "Pasture", "Cropland", "Urban"))


land_use_effect_plot <- ggplot(data = land_use_plot, aes(x = land_use, y = fit, colour = land_use )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin= lower , ymax= upper), colour="black", width=.1, position=position_dodge(0.5), linetype = 1) +
  geom_point(size = 5)  +
  theme_classic()


plot(land_use_effect_plot)












































linear_model_data <- hole_frame %>% dplyr::group_by(SSBS,SS) %>% 
  dplyr::summarise(total_hole_volume = max(total_hole_volume),
                   hypervolume_occupancy = hypervolume_occupancy,
                   land_use = land_use) %>%
  dplyr::distinct() %>% dplyr::left_join(TPD_land_uses[,c("SSBS","Realm")], by = "SSBS")



lmm <- lmer(total_hole_volume ~ 
              hypervolume_occupancy:Realm +
              hypervolume_occupancy +
              hypervolume_occupancy:land_use +
              land_use:Realm + land_use + (1|SS),
            data = linear_model_data)

summary(lmm)
car::Anova(lmm)


summary(gam)


gam <- read_rds("Outputs/site_level_gamm.rds")

summary(gam$mer)

unique(hole_frame$hypervolume_occupancy)

realm_plot_list <- list()
realm <- "Nearctic"
for(realm in realms){
  
  gamm_plot <-c() 
  
  for(land_use in land_uses){
    
    if(realm == "Nearctic" & land_use == "Secondary vegetation"){
      next()
    }
    
    hyper_size <- range(linear_model_data[linear_model_data$land_use == land_use &
                                            linear_model_data$Realm == realm, "hypervolume_occupancy"])
    hyper_size <- seq(hyper_size[1],hyper_size[2], length.out = 100)
    
    for(i in hyper_size){
      
      
      
      
      data <- data.frame(radius = i, land_use = land_use, 
                         mean_proportion = predict(lmm, newdata = data.frame(land_use = land_use, 
                                                                             hypervolume_occupancy = i,
                                                                             Realm = realm)
                                                   ,re.form = NA))
      
      gamm_plot <- rbind(gamm_plot, data)
    }
  }
  
  
  gamm_plot$land_use <- factor(gamm_plot$land_use, levels = c("Primary vegetation", "Secondary vegetation", "Plantation forest",
                                                              "Pasture", "Cropland", "Urban"))
  
  
  
  
  test <- ggplot(data = gamm_plot, aes(x = radius, y = mean_proportion, colour = land_use)) +
    geom_point()
  
  plot(test)
  
  
}





gam_number_of_holes <- mgcv::gamm(number_of_holes ~ s(radius, by = land_use),random = list(SS = ~1), data = hole_frame)


gamm_plot <- c() 
for(i in seq(0,1.5, 0.01)){
  for(land_use in land_uses){
    
    data <- data.frame(radius = i, land_use = land_use, 
                       mean_proportion = predict(gam_number_of_holes$gam, newdata = data.frame(radius = i, land_use = land_use)))
    
    gamm_plot <- rbind(gamm_plot, data)
  }
}

gamm_plot$land_use <- factor(gamm_plot$land_use, levels = c("Primary vegetation", "Secondary vegetation", "Plantation forest",
                                                            "Pasture", "Cropland", "Urban"))

test <- ggplot(data = gamm_plot, aes(x = radius, y = mean_proportion, colour = land_use)) +
  geom_point()

plot(test)

