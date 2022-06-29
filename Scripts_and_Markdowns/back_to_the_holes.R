land_uses <- c("Primary vegetation", "Secondary vegetation", "Plantation forest", "Pasture", "Cropland", "Urban")
realms <- c("Neotropic", "Afrotropic","Palearctic","Nearctic","Indo-Malay", "Australasia")


require(gamm4)
require(mgcv)

hole_frame <- c()


afro_sites <- TPD_LU %>% dplyr::filter(Realm == "Afrotropic")


for(sites in unique(afro_sites$SSBS)){
print(sites)  
  
  LU_holes <- TPD_holes(data = PREDICTS_tpds, randata = PREDICTS_randomisations, sites = sites, threshold = 1, minimum_points = 10)
  
  h_dat <- LU_holes$internal$metrics_frame
  if(is_null(h_dat)){
    next()
  }
  
  h_dat$land_use <- TPD_LU %>% dplyr::filter(SSBS == sites) %>% dplyr::distinct(Predominant_habitat) %>% pull()
  

  
  
  hole_frame <- rbind(hole_frame, h_dat)
    
}



hole_frame_plot <- hole_frame %>% group_by(radius, land_use) %>% dplyr::summarise(mean_holes = mean(number_of_holes),
                                                                                  mean_proportion = mean(mean_proportion),
                                                                                  mean_t_proportion = mean(total_proportion)) 

hole_frame_plot$land_use <- factor(hole_frame_plot$land_use, levels = c("Primary vegetation", "Secondary vegetation", "Plantation forest",
                                                              "Pasture", "Cropland", "Urban"))

gam <- mgcv::gam(mean_proportion ~ s(radius, by = land_use) , data = hole_frame_plot)

summary(gam)


plot(gam)


test <- ggplot(data = hole_frame_plot, aes(x = radius, y = mean_holes, colour = land_use)) +
  geom_point()

plot(test)



test_2 <- ggplot(data = hole_frame_plot, aes(x = radius, y = mean_proportion, colour = land_use)) +
  geom_point()

plot(test_2)
