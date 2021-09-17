
########################################
###### TPD 3D PLOT SCRIPT ##############
#######################################


###  Infographic plot to exemplify the characterisation and occupancy of trait space 

rm(list = ls())

require(gstat)
require(sf)
require(ggpubr) ## for multiple plots
require(magrittr) ## piping 
require(tidyverse) ## data manipulations
require(rgl) ## 3D plotting

options(rgl.printRglwidget = TRUE)

## loading in our TPDs and PREDICTS database -- PREDICTS will be used to get the site land-use classifications for sites that can be combined within study

PREDICTS_tpds <- readRDS("Outputs/PREDICTS_sites_tpds.rds")
PREDICTS_tpds_for <- readRDS("Outputs/PREDICTS_sites_for_tpds.rds")
PREDICTS <- readRDS("Outputs/refined_predicts.rds") %>% dplyr::distinct(SSBS, Predominant_habitat, Use_intensity)


### function to plot a 2D plot of pairwise trait axes input is a matrix with three columns the first two being trait axes values labelled T1 and T2,
### the third column being the TPD probability of occupancy


TPD_plot_fun <- function(mat){
  
  plot <- with(mat, ggplot2::ggplot(mat, ggplot2::aes(x = T1, 
                                                      y = T2, fill = prob), interpolate = TRUE)) + ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient(na.value = "grey80") + 
    ggplot2::theme(axis.line = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_text(colour = "black", 
                                                       size = 8), axis.text.y = ggplot2::element_text(colour = "black", 
                                                                                                      size = 8),  
                   legend.position = "none", 
                   plot.title = ggplot2::element_text(size = 10), 
                   plot.background = ggplot2::element_blank(), 
                   panel.background = ggplot2::element_rect(size = 1, 
                                                            linetype = "solid", color = "black"), 
                   panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank())
  
  return(plot)
  
}

#########################
######################## function to gather data in the required format for the figure making can handle multiple sites at once so when it comes to 
####################### getting sites of the same land-use type you can just input that as a single vector


TPD_plot_data <- function(data,site){
  
  pl_dat <- matrix(rep(0,(50^3)*4), ncol = 4)
  pl_dat[,1] <- data[[1]][["data"]][["evaluation_grid"]][[1]]
  pl_dat[,2] <- data[[1]][["data"]][["evaluation_grid"]][[2]]
  pl_dat[,3] <- data[[1]][["data"]][["evaluation_grid"]][[3]]
  
  for(sit in site){
    pl_dat[,4] <- pl_dat[,4] + data[[sit]][["TPDc"]][["RelativeAbundance"]]
  }
  
  
  pl_dat[,4] <- pl_dat[,4]/length(site)  
  colnames(pl_dat) <- c("T1","T2","T3","prob")  
  pl_dat <- data.frame(pl_dat)
  
  T21_dat <- pl_dat %>% dplyr::group_by(T2,T1) %>% dplyr::summarise(prob = mean(prob, na.rm = TRUE)) %>%
    dplyr::mutate(prob = ifelse(prob == 0,NA,prob*50))%>% set_colnames(c("T1","T2","prob")) %>% data.frame()
  T31_dat <- pl_dat %>% dplyr::group_by(T3,T1) %>% dplyr::summarise(prob = mean(prob, na.rm = TRUE)) %>%
    dplyr::mutate(prob = ifelse(prob == 0,NA,prob*50)) %>% set_colnames(c("T1","T2","prob")) %>% data.frame()
  T23_dat <- pl_dat %>% dplyr::group_by(T2,T3) %>% dplyr::summarise(prob = mean(prob, na.rm = TRUE)) %>%
    dplyr::mutate(prob = ifelse(prob == 0,NA,prob*50)) %>% set_colnames(c("T1","T2","prob")) %>% data.frame()
  
  
  trait_plot_data <- list()
  trait_plot_data$pl_dat <- pl_dat
  trait_plot_data$T21_dat <- T21_dat
  trait_plot_data$T31_dat <- T31_dat
  trait_plot_data$T23_dat <- T23_dat
  
  return(trait_plot_data)
  
  
}




################ examle of how to use the data and how to represent the three trait axes in 3 2D plots

# plot_data <- TPD_plot_data(data = PREDICTS_tpds,site)
# 
# ##
# 
# T21_plo <- TPD_plot_fun(plot_data[["T21_dat"]]) + ylab("Locomotion") + theme(axis.title.x = element_blank())
# T31_plo <- TPD_plot_fun(plot_data[["T31_dat"]]) + xlab("Body") + theme(axis.title.y = element_blank())
# T23_plo <- TPD_plot_fun(plot_data[["T23_dat"]]) + xlab("Foraging") + ylab("Body")
# 
# TPD_grid <-  ggarrange(T21_plo,T31_plo,T23_plo, ncol = 2, nrow = 2 , align = "hv")
# 
# plot(TPD_grid)


#########################################
#########################################
# Now try to represent the TPD in three dimensions
# constructing this plot comes in a number of stages and requires the use of the rgl package 


TPD_3d_plot <- function(data,sites,T1lab,T2lab,T3lab){
  
  data_3d <- TPD_plot_data(data,sites)
  
  # first identify which cells in 3D space are occupied 
  
  filled_cells <- data_3d[["pl_dat"]] %>% dplyr::group_by(T2,T1,T3) %>% dplyr::summarise(prob = mean(prob, na.rm = TRUE)) %>%
    dplyr::mutate(prob = ifelse(prob == 0,NA,prob)) %>% filter(!is.na(prob))
  
  x <- filled_cells$T2
  y <- filled_cells$T1
  z <- filled_cells$T3
  
  
  my.colors<-colorRampPalette(c("white", "orange", "green", "darkgreen")) #creates a function my.colors which interpolates n colors between blue, white and red
  color.df<-data.frame(prob=filled_cells$prob[order(filled_cells$prob)], color.name=my.colors(length(filled_cells$prob))) %>% distinct(prob,.keep_all = TRUE)#generates 2001 colors from the color ramp
  filled_cells_col <- filled_cells %>% dplyr::left_join(color.df)
  
  
  my_colour_2 <- colorRampPalette(c("yellow","orange","red"))
  
  #####################################
  
  
  #######################################
  ###########################################
  
  
  T21_col_mat <- matrix(rep("#FFFFFF",50*50),
                        nrow = 50, ncol = 50)
  
  
  T21_col <- matrix(rep(min(data_3d[["pl_dat"]]$T3),50*50),
                    nrow = 50, ncol = 50)
  rownames(T21_col) <- unique(data_3d[["T21_dat"]]$T1)
  colnames(T21_col) <- unique(data_3d[["T21_dat"]]$T2)
  
  for(i in 1:nrow(data_3d[["T21_dat"]])){
    T21_col[as.character(data_3d[["T21_dat"]][i,1]),as.character(data_3d[["T21_dat"]][i,2])] <- data_3d[["T21_dat"]][i,3]
  }
  
  T21_vals <- T21_col[which(!is.na(T21_col))]
  
  
  ###rows are x
  ### columns are y
  
  colours_T21 <- data.frame(value=T21_vals[order(T21_vals)], color.name=my_colour_2(length(T21_vals))) 
  
  T21_col_mat[which(!is.na(T21_col))[order(T21_vals)]] <- colours_T21$color.name
  
  
  T21_col_mat <- as.character(T21_col_mat)
  
  
  ##########################################################
  ############################################################
  
  
  
  T31_col_mat <- matrix(rep("#FFFFFF",50*50),
                        nrow = 50, ncol = 50)
  
  
  T31_col <- matrix(rep(unique(data_3d[["pl_dat"]]$T3),50),
                    nrow = 50, ncol = 50)
  rownames(T31_col) <- unique(data_3d[["T31_dat"]]$T1)
  colnames(T31_col) <- unique(data_3d[["T31_dat"]]$T2)
  
  for(i in 1:nrow(data_3d[["T31_dat"]])){
    T31_col[as.character(data_3d[["T31_dat"]][i,1]),as.character(data_3d[["T31_dat"]][i,2])] <- data_3d[["T31_dat"]][i,3]
  }
  
  
  
  T31_vals <- T31_col[which(!is.na(T31_col))]
  
  
  ###rows are x which is 
  ### columns are y
  
  colours_T31 <- data.frame(value=T31_vals[order(T31_vals)], color.name=my_colour_2(length(T31_vals)))
  
  
  T31_col_mat[which(!is.na(T31_col))[order(T31_vals)]] <- colours_T31$color.name
  
  
  T31_col_mat <- as.character(T31_col_mat)
  
  
  ###################################################
  ##################################################
  
  
  ##########################################################
  ############################################################
  
  
  
  T23_col_mat <- matrix(rep("#FFFFFF",50*50),
                        nrow = 50, ncol = 50)
  
  
  T23_col <- matrix(rep(unique(data_3d[["pl_dat"]]$T3),50),
                    nrow = 50, ncol = 50)
  rownames(T23_col) <- unique(data_3d[["T23_dat"]]$T1)
  colnames(T23_col) <- unique(data_3d[["T23_dat"]]$T2)
  
  for(i in 1:nrow(data_3d[["T23_dat"]])){
    T23_col[as.character(data_3d[["T23_dat"]][i,1]),as.character(data_3d[["T23_dat"]][i,2])] <- data_3d[["T23_dat"]][i,3]
  }
  
  
  
  T23_vals <- T23_col[which(!is.na(T23_col))]
  
  
  ###rows are x which is 
  ### columns are y
  
  colours_T23 <- data.frame(value=T23_vals[order(T23_vals)], color.name=my_colour_2(length(T23_vals)))
  
  
  T23_col_mat[which(!is.na(T23_col))[order(T23_vals)]] <- colours_T23$color.name
  
  
  T23_col_mat <- as.character(T23_col_mat)
  
  #######################################
  ###################################
  
  xmax <- max(data_3d[["pl_dat"]]$T2)
  xmin <- min(data_3d[["pl_dat"]]$T2)
  ymax <- max(data_3d[["pl_dat"]]$T1)
  ymin <- min(data_3d[["pl_dat"]]$T1)
  zmax <- max(data_3d[["pl_dat"]]$T3)
  zmin <- min(data_3d[["pl_dat"]]$T3)
  
  scale <- mean(c(dist(c(xmin,xmax))[1],
       dist(c(ymin,ymax))[1],
       dist(c(zmin,zmax))[1])) / 100
  
  
  plot3d(x, y, z, box = FALSE,xlab = "",ylab = "",zlab = "",
         type ="s", radius = scale,alpha = 0.8, xlim = c(xmin,xmax),
         ylim = c(ymin,ymax),
         zlim = c(zmin,zmax), col = filled_cells_col$color.name, lwd = 0.1)
  surface3d(x = unique(data_3d[["T21_dat"]]$T1),
            y = unique(data_3d[["T21_dat"]]$T2),
            z = matrix(rep(min(data_3d[["pl_dat"]]$T3),50*50),
                       nrow = 50, ncol = 50),
            alpha = 0.5, lit = FALSE, front = "lines", back = "lines")
  surface3d(x = unique(data_3d[["T21_dat"]]$T1),
            y = unique(data_3d[["T21_dat"]]$T2),
            z = matrix(rep(min(data_3d[["pl_dat"]]$T3),50*50),
                       nrow = 50, ncol = 50),
            color = T21_col_mat, smooth = FALSE, lit = FALSE)
  surface3d(x = rep(min(data_3d[["T21_dat"]]$T1),50),
            y = unique(data_3d[["T21_dat"]]$T2),
            z = matrix(rep(unique(data_3d[["pl_dat"]]$T3),50),
                       nrow = 50, ncol = 50),
            alpha = 0.5, lit = FALSE, front = "lines", back = "lines")
  surface3d(x = rep(min(data_3d[["T21_dat"]]$T1),50),
            y = unique(data_3d[["T21_dat"]]$T2),
            z = matrix(rep(unique(data_3d[["pl_dat"]]$T3),50),
                       nrow = 50, ncol = 50),
            color = T31_col_mat,lit = FALSE, smooth = FALSE)
  surface3d(x = unique(data_3d[["T21_dat"]]$T1),
            y = rep(min(data_3d[["T21_dat"]]$T2),50),
            z = t(matrix(rep(unique(data_3d[["pl_dat"]]$T3),50),
                         nrow = 50, ncol = 50)),
            alpha = 0.5, lit = FALSE, front = "lines", back = "lines")
  surface3d(x = unique(data_3d[["T21_dat"]]$T1),
            y = rep(min(data_3d[["T21_dat"]]$T2),50),
            z = t(matrix(rep(unique(data_3d[["pl_dat"]]$T3),50),
                         nrow = 50, ncol = 50)),
            color = T23_col_mat, lit = FALSE, smooth = FALSE)
  title3d(main = "3D_TPD_Plot", xlab = T2lab, ylab = T1lab, zlab = T3lab)
  
  
  
  
}




sites_lu <- data.frame(SSBS= names(PREDICTS_tpds)) %>% dplyr::left_join(PREDICTS[,c("SSBS","Predominant_habitat","Use_intensity")])

primary <- sites_lu %>% dplyr::filter(grepl(Predominant_habitat, pattern = "Primary")) %>% pull(SSBS)



TPD_3d_plot(sites = primary,data = PREDICTS_tpds, T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body")

#####################################################
#####################################################
## SD plot showing the difference between two set of sites 

secondary <- sites_lu %>% dplyr::filter(grepl(Predominant_habitat, pattern = "secondary",ignore.case = TRUE)) %>% pull(SSBS)

data <- PREDICTS_tpds
sites1 <- primary
sites2 <- secondary

TPD_Diff_Func <- function(data,sites1,sites2,T1lab,T2lab,T3lab){
  
  sites1_data <- TPD_plot_data(data,sites1)
  sites2_data <- TPD_plot_data(data,sites2)
 
  filled_cells_1 <- sites1_data[["pl_dat"]] %>% dplyr::group_by(T2,T1,T3) %>% dplyr::summarise(prob = mean(prob, na.rm = TRUE)) %>%
    dplyr::mutate(prob = ifelse(prob == 0,NA,prob)) %>% filter(!is.na(prob))
  
  filled_cells_2 <- sites2_data[["pl_dat"]] %>% dplyr::group_by(T2,T1,T3) %>% dplyr::summarise(prob_2 = mean(prob, na.rm = TRUE)) %>%
    dplyr::mutate(prob_2 = ifelse(prob_2 == 0,NA,prob_2)) %>% filter(!is.na(prob_2))
  
  diff_cells <- filled_cells_1 %>% dplyr::left_join(filled_cells_2, by = c("T1","T2","T3")) %>% data.frame()
  
  cells_frame <- filled_cells_2 %>% dplyr::left_join(filled_cells_1, by = c("T1","T2","T3")) %>% dplyr::filter(is.na(prob)) %>%
    dplyr::relocate(prob, .before = prob_2) %>% rbind(diff_cells) %>% dplyr::mutate(prob = ifelse(is.na(prob), 0, prob),
                                                                                    prob_2 = ifelse(is.na(prob_2), 0, prob_2),
                                                                                    diff = prob_2 - prob) 
    
  
  
  x <- cells_frame$T2
  y <- cells_frame$T1
  z <- cells_frame$T3
  
  
  ##### HIGH CELLS 
  
  pos_cells <- cells_frame %>% dplyr::filter(diff > 0) %>% pull(diff)
  
  my.colors_high <-colorRampPalette(c("tomato", "orangered", "orangered4", "red")) #creates a function my.colors which interpolates n colors between blue, white and red
  color.df_high<-data.frame(diff=pos_cells[order(pos_cells)], color.name=my.colors_high(length(pos_cells))) %>% distinct(diff,.keep_all = TRUE)#generates 2001 colors from the color ramp
  
  
  ##### LOW CELLS 
  
  neg_cells <- cells_frame %>% dplyr::filter(diff < 0) %>% pull(diff)
  
  my.colors_low <-colorRampPalette(c("blue", "steelblue4", "steelblue1", "lightskyblue1")) #creates a function my.colors which interpolates n colors between blue, white and red
  color.df_low<-data.frame(diff=neg_cells[order(neg_cells)], color.name=my.colors_low(length(neg_cells))) %>% distinct(diff,.keep_all = TRUE)#generates 2001 colors from the color ramp
  
  color.df <- rbind(color.df_high,color.df_low)
  
  filled_cells_col <- cells_frame %>% dplyr::left_join(color.df, by = "diff") %>% 
    dplyr::mutate(color.name = ifelse(is.na(color.name), "#FFFFFF",color.name))
  
  

 
  ############################################################
  ############################################################
  
  
  my_colour_2p <- colorRampPalette(c("yellow","orange","red"))
  my_colour_2n <- colorRampPalette(c("blue","steelblue1","lightskyblue1"))
  
  #######################################
  ###########################################
  
  
  T21_col_mat <- matrix(rep("#FFFFFF",50*50),
                        nrow = 50, ncol = 50)
  
  
  T21_col <- matrix(rep(min(sites1_data[["pl_dat"]]$T3),50*50),
                    nrow = 50, ncol = 50)
  rownames(T21_col) <- unique(sites1_data[["T21_dat"]]$T1)
  colnames(T21_col) <- unique(sites1_data[["T21_dat"]]$T2)
  
  s1d <- sites1_data[["T21_dat"]] %>% dplyr::mutate(prob = ifelse(is.na(prob), 0, prob))
  s2d <- sites2_data[["T21_dat"]] %>% dplyr::mutate(prob = ifelse(is.na(prob), 0, prob))
  
  for(i in 1:nrow(sites1_data[["T21_dat"]])){
    T21_col[as.character(sites1_data[["T21_dat"]][i,1]),as.character(sites1_data[["T21_dat"]][i,2])] <- 
      ifelse((s2d[i,3] - s1d[i,3]) == 0, NA,(s2d[i,3] - s1d[i,3]))
  }
  
  T21_vals <- T21_col[which(!is.na(T21_col))]
  
  pos_T21_vals <- T21_vals[which(T21_vals > 0)]
  pos_T21_vals <- pos_T21_vals[order(pos_T21_vals)]
  pos_T21_cols <- my_colour_2p(length(pos_T21_vals))

  neg_T21_vals <- T21_vals[which(T21_vals < 0)]
  neg_T21_vals <- neg_T21_vals[order(neg_T21_vals)]
  neg_T21_cols <- my_colour_2n(length(neg_T21_vals))
  
  colours_T21 <- data.frame(value = c(neg_T21_vals,pos_T21_vals), color.name = c(neg_T21_cols,pos_T21_cols)) 
  
  ###rows are x
  ### columns are y
  
  
  T21_col_mat[which(!is.na(T21_col))[order(T21_vals)]] <- colours_T21$color.name
  
  
  T21_col_mat <- as.character(T21_col_mat)
  
  
  #######################################
  ###########################################
  
  
  T31_col_mat <- matrix(rep("#FFFFFF",50*50),
                        nrow = 50, ncol = 50)
  
  
  T31_col <- matrix(rep(min(sites1_data[["pl_dat"]]$T3),50*50),
                    nrow = 50, ncol = 50)
  rownames(T31_col) <- unique(sites1_data[["T31_dat"]]$T1)
  colnames(T31_col) <- unique(sites1_data[["T31_dat"]]$T2)
  
  s1d <- sites1_data[["T31_dat"]] %>% dplyr::mutate(prob = ifelse(is.na(prob), 0, prob))
  s2d <- sites2_data[["T31_dat"]] %>% dplyr::mutate(prob = ifelse(is.na(prob), 0, prob))
  
  for(i in 1:nrow(sites1_data[["T31_dat"]])){
    T31_col[as.character(sites1_data[["T31_dat"]][i,1]),as.character(sites1_data[["T31_dat"]][i,2])] <- 
      ifelse((s2d[i,3] - s1d[i,3]) == 0, NA,(s2d[i,3] - s1d[i,3]))
  }
  
  T31_vals <- T31_col[which(!is.na(T31_col))]
  
  pos_T31_vals <- T31_vals[which(T31_vals > 0)]
  pos_T31_vals <- pos_T31_vals[order(pos_T31_vals)]
  pos_T31_cols <- my_colour_2p(length(pos_T31_vals))
  
  neg_T31_vals <- T31_vals[which(T31_vals < 0)]
  neg_T31_vals <- neg_T31_vals[order(neg_T31_vals)]
  neg_T31_cols <- my_colour_2n(length(neg_T31_vals))
  
  colours_T31 <- data.frame(value = c(neg_T31_vals,pos_T31_vals), color.name = c(neg_T31_cols,pos_T31_cols)) 
  
  ###rows are x
  ### columns are y
  
  
  T31_col_mat[which(!is.na(T31_col))[order(T31_vals)]] <- colours_T31$color.name
  
  
  T31_col_mat <- as.character(T31_col_mat)
  
  #######################################
  ###########################################
  
  
  T23_col_mat <- matrix(rep("#FFFFFF",50*50),
                        nrow = 50, ncol = 50)
  
  
  T23_col <- matrix(rep(min(sites1_data[["pl_dat"]]$T3),50*50),
                    nrow = 50, ncol = 50)
  rownames(T23_col) <- unique(sites1_data[["T23_dat"]]$T1)
  colnames(T23_col) <- unique(sites1_data[["T23_dat"]]$T2)
  
  s1d <- sites1_data[["T23_dat"]] %>% dplyr::mutate(prob = ifelse(is.na(prob), 0, prob))
  s2d <- sites2_data[["T23_dat"]] %>% dplyr::mutate(prob = ifelse(is.na(prob), 0, prob))
  
  for(i in 1:nrow(sites1_data[["T23_dat"]])){
    T23_col[as.character(sites1_data[["T23_dat"]][i,1]),as.character(sites1_data[["T23_dat"]][i,2])] <- 
      ifelse((s2d[i,3] - s1d[i,3]) == 0, NA,(s2d[i,3] - s1d[i,3]))
  }
  
  
  T23_vals <- T23_col[which(!is.na(T23_col))]
  
  pos_T23_vals <- T23_vals[which(T23_vals > 0)]
  pos_T23_vals <- pos_T23_vals[order(pos_T23_vals)]
  pos_T23_cols <- my_colour_2p(length(pos_T23_vals))
  
  neg_T23_vals <- T23_vals[which(T23_vals < 0)]
  neg_T23_vals <- neg_T23_vals[order(neg_T23_vals)]
  neg_T23_cols <- my_colour_2n(length(neg_T23_vals))
  
  colours_T23 <- data.frame(value = c(neg_T23_vals,pos_T23_vals), color.name = c(neg_T23_cols,pos_T23_cols)) 
  
  ###rows are x
  ### columns are y
  
  
  T23_col_mat[which(!is.na(T23_col))[order(T23_vals)]] <- colours_T23$color.name
  
  
  T23_col_mat <- as.character(T23_col_mat)
  
  
  #######################################
  ###################################
  
  xmax <- max(sites1_data[["pl_dat"]]$T2)
  xmin <- min(sites1_data[["pl_dat"]]$T2)
  ymax <- max(sites1_data[["pl_dat"]]$T1)
  ymin <- min(sites1_data[["pl_dat"]]$T1)
  zmax <- max(sites1_data[["pl_dat"]]$T3)
  zmin <- min(sites1_data[["pl_dat"]]$T3)
  
  scale <- mean(c(dist(c(xmin,xmax))[1],
                  dist(c(ymin,ymax))[1],
                  dist(c(zmin,zmax))[1])) / 100
  
  
  
  plot3d(x, y, z, box = FALSE,xlab = "",ylab = "",zlab = "",
         type ="s", radius = scale,alpha = 0.8, xlim = c(xmin,xmax),
         ylim = c(ymin,ymax),
         zlim = c(zmin,zmax), col = filled_cells_col$color.name, lwd = 0.1)
  surface3d(x = unique(sites1_data[["T21_dat"]]$T1),
            y = unique(sites1_data[["T21_dat"]]$T2),
            z = matrix(rep(min(sites1_data[["pl_dat"]]$T3),50*50),
                       nrow = 50, ncol = 50),
            alpha = 0.5, lit = FALSE, front = "lines", back = "lines")
  surface3d(x = unique(sites1_data[["T21_dat"]]$T1),
            y = unique(sites1_data[["T21_dat"]]$T2),
            z = matrix(rep(min(sites1_data[["pl_dat"]]$T3),50*50),
                       nrow = 50, ncol = 50),
            color = T21_col_mat, smooth = FALSE, lit = FALSE)
  surface3d(x = rep(min(sites1_data[["T21_dat"]]$T1),50),
            y = unique(sites1_data[["T21_dat"]]$T2),
            z = matrix(rep(unique(sites1_data[["pl_dat"]]$T3),50),
                       nrow = 50, ncol = 50),
            alpha = 0.5, lit = FALSE, front = "lines", back = "lines")
  surface3d(x = rep(min(sites1_data[["T21_dat"]]$T1),50),
            y = unique(sites1_data[["T21_dat"]]$T2),
            z = matrix(rep(unique(sites1_data[["pl_dat"]]$T3),50),
                       nrow = 50, ncol = 50),
            color = T31_col_mat,lit = FALSE, smooth = FALSE)
  surface3d(x = unique(sites1_data[["T21_dat"]]$T1),
            y = rep(min(sites1_data[["T21_dat"]]$T2),50),
            z = t(matrix(rep(unique(sites1_data[["pl_dat"]]$T3),50),
                         nrow = 50, ncol = 50)),
            alpha = 0.5, lit = FALSE, front = "lines", back = "lines")
  surface3d(x = unique(sites1_data[["T21_dat"]]$T1),
            y = rep(min(sites1_data[["T21_dat"]]$T2),50),
            z = t(matrix(rep(unique(sites1_data[["pl_dat"]]$T3),50),
                         nrow = 50, ncol = 50)),
            color = T23_col_mat, lit = FALSE, smooth = FALSE)
  title3d(main = "3D_TPD_Plot", xlab = T2lab, ylab = T1lab, zlab = T3lab)
  
  
}

TPD_Diff_Func(sites1 = primary,sites2 = secondary,data = PREDICTS_tpds, T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body")







