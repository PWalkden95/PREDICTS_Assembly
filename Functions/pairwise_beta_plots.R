beta_data <- beta_boot
realm <- "Neotropic"


pairwise_beta_plot <- function(beta_data, realm){

plot_colours <- c(`Primary forest_Primary forest` = "green4",
                  `Primary forest_Primary_non forest` = "charteuse4",
                  `Primary forest_Secondary vegetation` = "chartreuse3",
                  `Primary forest_Plantation forest` = "seagreen3",
                  `Primary forest_Pasture` = "darkseagreen1",
                  `Primary forest_Cropland` = "darkolivegreen1",
                  `Primary forest_Minimal agriculture` = "darkseagreen1",
                  `Primary forest_Intensive agriculture` = "darkolivegreen1",
                  `Primary forest_Urban` = "palegreen",
                  `Primary non-forest_Primary non-forest` = "green4",
                  `Primary non-forest_Secondary vegetation` = "chartreuse3",
                  `Primary non-forest_Plantation forest` = "seagreen3",
                  `Primary non-forest_Pasture` = "darkseagreen1",
                  `Primary non-forest_Cropland` = "darkolivegreen1",
                  `Primary non-forest_Minimal agriculture` = "darkseagreen1",
                  `Primary non-forest_Intensive agriculture` = "darkolivegreen1",
                  `Primary non-forest_Urban` = "palegreen",
                  `Primary vegetation_Primary vegetation` = "green4",
                  `Primary vegetation_Secondary vegetation` = "chartreuse3",
                  `Primary vegetation_Plantation forest` = "seagreen3",
                  `Primary vegetation_Pasture` = "darkseagreen1",
                  `Primary vegetation_Cropland` = "darkolivegreen1",
                  `Primary vegetation_Minimal agriculture` = "darkseagreen1",
                  `Primary vegetation_Intensive agriculture` = "darkolivegreen1",
                  `Primary vegetation_Urban` = "palegreen",
                  `Secondary vegetation_Secondary vegetation` = "olivedrab3",
                  `Secondary vegetation_Plantation forest` = "lawngreen",
                  `Secondary vegetation_Pasture` = "darkseagreen2",
                  `Secondary vegetation_Cropland` = "olivedrab1",
                  `Secondary vegetation_Minimal agriculture` = "darkseagreen2",
                  `Secondary vegetation_Intensive agriculture` = "olivedrab1",
                  `Secondary vegetation_Urban` = "lightskyblue",
                  `Plantation forest_Plantation forest` ="lightgreen",
                  `Plantation forest_Pasture` = "darkolivegreen2",
                  `Plantation forest_Cropland` = "seagreen2",
                  `Plantation forest_Minimal agriculture` = "darkolivegreen2",
                  `Plantation forest_Intensive agriculture` = "seagreen2",
                  `Plantation forest_Urban` = "mediumseagreen",
                  Pasture_Pasture = "khaki",
                  `Pasture_Cropland` = "goldenrod3",
                  `Minimal agriculture_Intensive agriculture` = "goldenrod3",
                  `Minimal agriculture_Urban` = "darkgoldenrod4",
                  `Pasture_Urban` = "darkgoldenrod4",
                  Cropland_Cropland ="gold1",
                  Cropland_Urban = "brown",
                  `Intensive agriculture_Urban` = "brown",
                  Urban_Urban = "ivory4")

shorthand <- c(`Primary forest` = "PriFor",
               `Primary non-forest` = "PriNFor",
  `Primary vegetation` = "PriVeg",
               `Secondary vegetation` = "SecVeg",
               `Plantation forest` = "PlnFor",
               `Pasture` = "Pas",
               `Cropland` = "Crp",
               `Minimal agriculture` = "MinAgr",
               `Intensive agriculture` = "InAgr",
               `Urban` = "Urb")



p_order <- list("1" = c(1),
                "2" = c(1,2,3),
                "3" = c(1,2,4,3,5,6),
                "4" = c(1,2,5,3,6,8,4,7,9,10),
                "5" = c(1,2,6,3,7,10,4,8,11,13,5,9,12,14,15),
                "6" = c(1,2,7,3,8,12,4,9,13,16,5,10,14,17,19,6,11,15,18,20,21),
                "7" = c(1,2,8,3,9,14,4,10,15,19,5,11,16,20,23,6,12,17,21,24,26,7,13,18,22,25,27,28))

s_order <- list("2" = c(2),
                "3" = c(2,3,5),
                "4" = c(2,3,4,6,7,9),
                "5" = c(2,3,4,5,7,8,9,11,12,14),
                "6" = c(2,3,4,5,6,8,9,10,11,13,14,15,17,18,20),
                "7" = c(2,3,4,5,6,7,9,10,11,12,13,15,16,17,18,20,21,22,24,25,27))


p_dat <- beta_data %>% dplyr::filter(Realm == realm)



plot_LU <- c("Primary forest","Primary non-forest", "Secondary vegetation", "Plantation forest", "Pasture", "Cropland", "Urban")

data_LU <- unique(unlist(str_split(p_dat$Land_use_comp,pattern =  "_")))


plot_LU <- plot_LU[plot_LU %in% data_LU]

combos <- c()
for(i in 1:length(plot_LU)){
  combos <- rbind(combos, matrix(rep(plot_LU[i], 2), ncol = 2))
}


combos <- rbind(combos,gtools::combinations(n = length(plot_LU), r = 2, v = plot_LU) %>% reorder_combinations())

plot_combo <- paste(combos[,1], combos[,2], sep = "_")




lvls <- beta_levels[beta_levels %in% plot_combo]

plot_lvls <- lvls[p_order[[as.character(length(plot_LU))]]]



level_plots <- list()

pri_median <- median(p_dat %>% dplyr::filter(Land_use_comp == "Primary forest_Primary forest") %>% pull(value))


for(l in plot_lvls){
  
  pd <- p_dat %>% dplyr::filter(Land_use_comp == l)

  if(nrow(pd) == 0){
     beta_plot <- ggplot() +
       theme(axis.title.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank(),
             axis.title.x=element_blank(),
             axis.text.x=element_blank(),
             axis.ticks.x=element_blank(),
             panel.background = element_rect(fill = "black", colour = "black", size = 0.5))
     
    } else {
  
  col <- plot_colours[l]
  
  beta_plot <- ggplot(data = pd, aes(x = 0.1, y = value)) +
    geom_flat_violin( fill = col, alpha = 0.6, colour = "white") +
    #geom_boxplot(width = 0.1, color="black", alpha=0.2, show.legend = FALSE, outlier.alpha = 0) +
    ylim(min(p_dat$value),max(p_dat$value)) +
    geom_hline(yintercept = pri_median, linetype = "dashed", lwd = 1.5) +
    coord_flip()  + 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background = element_blank())
  
  }
  level_plots[[l]] <- beta_plot
  
  }

############################
############################
# Get significances

## sig_plots


none <- ggplot() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank())



####################################################



contrast_1 <- unlist(str_split(pairwise_sig$contrast, pattern = " - "))[seq(1,nrow(pairwise_sig)*2,2)]
contrast_2 <- unlist(str_split(pairwise_sig$contrast, pattern = " - "))[seq(2,nrow(pairwise_sig)*2,2)]


realm_rows <- which(grepl(contrast_1, pattern = realm) &
                      grepl(contrast_1, pattern = "Primary forest_Primary forest") &
                      grepl(contrast_2, pattern = realm))

contrast_1 <- contrast_1[realm_rows]
contrast_2 <- contrast_2[realm_rows]

realm_conts <- pairwise_sig %>% dplyr::slice(realm_rows) 

diag_lvls <- paste(plot_LU, plot_LU, sep = "_")[-1]

for(d in diag_lvls){
  sig <- realm_conts[which(grepl(contrast_2, pattern = d)),"p.value"]
  est <- realm_conts[which(grepl(contrast_2, pattern = d)),"estimate"]

  
  lvl_median <- median(p_dat %>% dplyr::filter(Land_use_comp == d) %>% pull(value))
  
  lvl_sig <- level_plots[[d]]
  

  if(is.na(sig)){
    next()
  }
  
      if(sig >= 0.05){
    } else {
      if(sig < 0.05 & sig > 0.01 ){
        lvl_sig <- lvl_sig +
          annotate("text", label = "*", x = 0.45, y = lvl_median+0.05, size = 10, colour = ifelse(est < 0, "blue", "red"),fontface = "bold")
        level_plots[[d]] <- lvl_sig
      } else {
        if(sig < 0.01 & sig > 0.001){
          lvl_sig <- lvl_sig +
            annotate("text", label = "**", x = 0.45, y = lvl_median+0.05, size = 10, colour = ifelse(est < 0, "blue", "red"),fontface = "bold")
          level_plots[[d]] <- lvl_sig
        } else {
          if(sig < 0.001){
            lvl_sig <- lvl_sig +
              annotate("text", label = "***", x = 0.45, y = lvl_median+0.05, size = 10, colour = ifelse(est < 0, "blue", "red"), fontface = "bold")
            level_plots[[d]] <- lvl_sig
          }
        }
      }
    }
    
  

  }



sig_lvls <- lvls[s_order[[as.character(length(plot_LU))]]]


significance_plots <- list()

for(l in sig_lvls){
  
  sig <- realm_conts[which(grepl(contrast_2, pattern = l)),"p.value"]
  est <- realm_conts[which(grepl(contrast_2, pattern = l)),"estimate"]
  
  if(is_empty(sig)){
    blank <- none +
      theme(panel.background = element_rect(fill = "black", colour = "black", size = 0.5))
    significance_plots <- c(significance_plots, list(blank))
  } else {
  
  
  if(is.na(sig)){
    blank <- none +
      theme(panel.background = element_rect(fill = "black", colour = "black", size = 0.5))
    significance_plots <- c(significance_plots, list(blank))
  } else {
  
  if(sig >= 0.05){
    significance_plots <- c(significance_plots,list(none))
  } else {
    if(sig < 0.05 & sig > 0.01 ){
      one <- none +
        xlim(0,1) +
        ylim(0,1) +
        annotate("text", label = "*", x = 0.5, y = 0.3, size = 15, colour = ifelse(est < 0, "blue", "red"),fontface = "bold")
      significance_plots <- c(significance_plots, list(one))
    } else {
      if(sig < 0.01 & sig > 0.001){
        two <- none +
          xlim(0,1) +
          ylim(0,1) +
          annotate("text", label = "**", x = 0.5, y = 0.3, size = 15, colour = ifelse(est < 0, "blue", "red"),fontface = "bold")
        significance_plots <- c(significance_plots, list(two))
      } else {
        if(sig < 0.001){
          three <- none +
            xlim(0,1) +
            ylim(0,1) +
            annotate("text", label = "***", x = 0.5, y = 0.3, size = 15, colour = ifelse(est < 0, "blue", "red"),fontface = "bold")
          significance_plots <- c(significance_plots, list(three))
        }
      }
    }
  }
    
  }
  }
} 

############################################
############################################




final_plots <- list()

pstart <- 1  
pend <- 0


sstart <- 1

s_length <- length(plot_LU) -1


for(i in 1:length(plot_LU)){
  pend <- pend + i
  
  s_leng <- s_length - i
  
send <- sstart + s_leng
  
  

if(i == length(plot_LU)){
  final_plots <- c(final_plots, level_plots[pstart:pend])
} else {
  final_plots <- c(final_plots, level_plots[pstart:pend],significance_plots[sstart:send])
}

  pstart <- pend + 1
  sstart <- send + 1
  
   
  
  }


figure <- ggarrange(plotlist = final_plots, ncol = length(plot_LU), nrow = length(plot_LU))


# 
# 
# ypos <- seq(0,1,by = 1/length(plot_LU))
# ypos <- ypos[-length(ypos)]
# 
# xpos <- ypos
# xpos <- xpos + (1/length(plot_LU)/2)
# 
# ypos <- ypos[order(ypos, decreasing = TRUE)]
# ypos <- ypos + (1/length(plot_LU)/2)
# 
# 
# plot_short <- shorthand[plot_LU]
# 
# figure <-  annotate_figure(figure,
#                             left = text_grob(plot_short[1], face = "bold", y = ypos[1], size = 20))
# 
# 
# for(i in 2:length(plot_LU)){
# figure <- figure +
#   annotate(geom = "text", y = ypos[i], x = 0.08, label = plot_short[i], fontface = "bold", size = 8)
# 
# }
# 
# 
# 
# xstart <- 0.1 + (0.9/length(plot_LU))/2
# 
# figure <- annotate_figure(figure,
#                 top = text_grob(plot_short[1], face = "bold", rot = 90, size = 20, x = xstart))
# 
# for(i in 2:length(plot_LU)){
#   
#   xstart <- xstart + (0.9/length(plot_LU))
#   
#   figure <- figure +
#     annotate(geom = "text", y = 0.87, x = xstart, label = plot_short[i], fontface = "bold", size = 8, angle = 90)
# 
# }



return(figure)

}