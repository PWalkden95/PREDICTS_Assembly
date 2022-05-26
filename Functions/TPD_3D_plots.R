
markdown_rds_open <- function(path) {
  open_file <- try(readRDS(path))
  
  if (class(open_file) == "try-error") {
    open_file <- readRDS(paste("../", path, sep = ""))
  }
  
  return(open_file)
  
}


########################################
###### TPD 3D PLOT SCRIPT ##############
#######################################


###  Infographic plot to exemplify the characterisation and occupancy of trait space
###


require(gstat)
require(sf)
require(ggpubr) ## for multiple plots
require(magrittr) ## piping
require(tidyverse) ## data manipulations
require(rgl) ## 3D plotting
require(geometry)
require(fastcluster)
require(Hmisc)
require(magick)
library(webshot2)

TPD_colours <- markdown_rds_open("Functions/TPD_colours.rds")


find_position <- function(x,y){
  
  value <- which(x > y)
  
  if(is_empty(value)){
    value <- 1
  } else {
  value <- value[length(value)]
  }
  
  return(value)   
}

options(rgl.printRglwidget = TRUE)

## loading in our TPDs and PREDICTS database -- PREDICTS will be used to get the site land-use classifications for sites that can be combined within study
#
# PREDICTS_tpds <- readRDS("Outputs/PREDICTS_sites_tpds.rds")
# PREDICTS_tpds_for <- readRDS("Outputs/PREDICTS_sites_for_tpds.rds")
# PREDICTS <- readRDS("Outputs/refined_predicts.rds") %>% dplyr::distinct(SSBS, Predominant_habitat, Use_intensity)
#


### percentileTPD: function to transform probabilities from TPDc object into quantiles
### x is a TPDc object
percentile_cells <- function(data) {
  x <- data
  x$percentile <- NA
  x$original_rows <- as.numeric(rownames(x))
  x <- x %>% dplyr::arrange(desc(prob))
  x[, "percentile"] <- cumsum(x[, "prob"])
  x <- x %>% dplyr::arrange(original_rows)
  x <- x %>% dplyr::select(-original_rows)
  return(x)
}
### function to plot a 2D plot of pairwise trait axes input is a matrix with three columns the first two being trait axes values labelled T1 and T2,
### the third column being the TPD probability of occupancy


TPD_plot_fun <- function(mat) {
  plot <- with(mat, ggplot2::ggplot(mat, ggplot2::aes(
    x = T1,
    y = T2, fill = prob
  ), interpolate = TRUE)) + ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient(na.value = "grey80") +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(colour = "black",
                                          size = 8),
      axis.text.y = ggplot2::element_text(colour = "black",
                                          size = 8),
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 10),
      plot.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        size = 1,
        linetype = "solid",
        color = "black"
      ),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank()
    )
  
  return(plot)
  
}

#########################
######################## function to gather data in the required format for the figure making can handle multiple sites at once so when it comes to
####################### getting sites of the same land-use type you can just input that as a single vector


TPD_plot_data <- function(data, site) {
  pl_dat <- matrix(rep(0, (50 ^ 3) * 4), ncol = 4)
  pl_dat[, 1] <- data[["data"]][["evaluation_grid"]][[1]]
  pl_dat[, 2] <- data[["data"]][["evaluation_grid"]][[2]]
  pl_dat[, 3] <- data[["data"]][["evaluation_grid"]][[3]]
  
  for (sit in site) {
    pl_dat[, 4] <-
      pl_dat[, 4] + data[[sit]][["TPDc"]][["RelativeAbundance"]]
  }
  
  
  pl_dat[, 4] <- pl_dat[, 4] / length(site)
  colnames(pl_dat) <- c("T1", "T2", "T3", "prob")
  pl_dat <- data.frame(pl_dat)
  
  T21_dat <-
    pl_dat %>% dplyr::group_by(T2, T1) %>% dplyr::summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% set_colnames(c("T1", "T2", "prob")) %>% data.frame() %>% percentile_cells() %>%
    dplyr::filter(!is.na(T1))
  T31_dat <-
    pl_dat %>% dplyr::group_by(T3, T1) %>% dplyr::summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% set_colnames(c("T1", "T2", "prob")) %>% data.frame() %>% percentile_cells() %>%
    dplyr::filter(!is.na(T1))
  T23_dat <-
    pl_dat %>% dplyr::group_by(T2, T3) %>% dplyr::summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% set_colnames(c("T1", "T2", "prob")) %>% data.frame() %>% percentile_cells() %>%
    dplyr::filter(!is.na(T1))
  
  
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

hypervolume_occupied_cells <- function(data){
  occupied_cells <- data %>% dplyr::group_by(T1,T2,T3) %>%
    dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% dplyr::filter(!is.na(prob)) %>% data.frame
  
  return(occupied_cells)
}


TPD_3d_plot <-
  function(data,
           sites = NULL,
           T1lab = "",
           T2lab = "",
           T3lab = "",
           method = "prob",
           title = "",
           save = FALSE,
           file,
           limits = limits,
           scale = scale,
           grid = FALSE,
           free_limits = FALSE) {
    
    
      if(!all(sites %in% names(data))) {
        sites <- sites[which(sites %in% names(data))]
      }
      
      
      data_3d <- TPD_plot_data(data, sites)
      
      # first identify which cells in 3D space are occupied
      
      filled_cells <- hypervolume_occupied_cells(data_3d[["pl_dat"]]) %>% percentile_cells()
      
      
      
      if (method == "prob") {
        filled_cells <- filled_cells %>% dplyr::rename(value = prob)
        my.colors <-
          colorRampPalette(c("white", "orange", "green", "darkgreen"))
        my_colour_2 <- colorRampPalette(c("yellow", "orange", "red"))
      } else {
        filled_cells <- filled_cells %>% dplyr::rename(value = percentile)
        my.colors <-
          colorRampPalette(c("darkgreen", "green", "orange", "white"))
        my_colour_2 <- colorRampPalette(c("red", "orange", "yellow"))
      }
      
      
      x <- filled_cells$T2
      y <- filled_cells$T1
      z <- filled_cells$T3
      
      
      
      
      
      #creates a function my.colors which interpolates n colors between blue, white and red
      color.df <-
        data.frame(value = TPD_colours[order(TPD_colours)], color.name =
                     my.colors(length(TPD_colours))) %>% distinct(value, .keep_all = TRUE)#generates 2001 colors from the color ramp
     
      
      filled_cells_col <- filled_cells %>% 
        dplyr::mutate(color.name = color.df$color.name[unlist(lapply(filled_cells$value,find_position, y = color.df$value))])
      

  
      
      #####################################
      
      
      #######################################
      ###########################################
      if (grid) {
        T21_col_mat <- matrix(rep("#FFFFFF", 50 * 50),
                              nrow = 50, ncol = 50)
        
        
        T21_col <- matrix(rep(min(data_3d[["pl_dat"]]$T3), 50 * 50),
                          nrow = 50, ncol = 50)
        rownames(T21_col) <- unique(data_3d[["T21_dat"]]$T1)
        colnames(T21_col) <- unique(data_3d[["T21_dat"]]$T2)
        
        for (i in 1:nrow(data_3d[["T21_dat"]])) {
          T21_col[as.character(data_3d[["T21_dat"]][i, 1]), as.character(data_3d[["T21_dat"]][i, 2])] <-
            ifelse(method == "prob",
                   data_3d[["T21_dat"]][i, 3],
                   data_3d[["T21_dat"]][i, 4])
          
        }
        
        T21_vals <- T21_col[which(!is.na(T21_col))]
        
        
        ###rows are x
        ### columns are y
        
        colours_T21 <-
          data.frame(value = T21_vals[order(T21_vals)], color.name = my_colour_2(length(T21_vals)))
        
        T21_col_mat[which(!is.na(T21_col))[order(T21_vals)]] <-
          colours_T21$color.name
        
        
        T21_col_mat <- as.character(T21_col_mat)
        
        
        ##########################################################
        ############################################################
        
        
        
        T31_col_mat <- matrix(rep("#FFFFFF", 50 * 50),
                              nrow = 50, ncol = 50)
        
        
        T31_col <- matrix(rep(unique(data_3d[["pl_dat"]]$T3), 50),
                          nrow = 50, ncol = 50)
        rownames(T31_col) <- unique(data_3d[["T31_dat"]]$T1)
        colnames(T31_col) <- unique(data_3d[["T31_dat"]]$T2)
        
        for (i in 1:nrow(data_3d[["T31_dat"]])) {
          T31_col[as.character(data_3d[["T31_dat"]][i, 1]), as.character(data_3d[["T31_dat"]][i, 2])] <-
            ifelse(method == "prob",
                   data_3d[["T31_dat"]][i, 3],
                   data_3d[["T31_dat"]][i, 4])
        }
        
        
        
        
        T31_vals <- T31_col[which(!is.na(T31_col))]
        
        
        ###rows are x which is
        ### columns are y
        
        colours_T31 <-
          data.frame(value = T31_vals[order(T31_vals)], color.name = my_colour_2(length(T31_vals)))
        
        
        T31_col_mat[which(!is.na(T31_col))[order(T31_vals)]] <-
          colours_T31$color.name
        
        
        T31_col_mat <- as.character(T31_col_mat)
        
        
        ###################################################
        ##################################################
        
        
        ##########################################################
        ############################################################
        
        
        
        T23_col_mat <- matrix(rep("#FFFFFF", 50 * 50),
                              nrow = 50, ncol = 50)
        
        
        T23_col <- matrix(rep(unique(data_3d[["pl_dat"]]$T3), 50),
                          nrow = 50, ncol = 50)
        rownames(T23_col) <- unique(data_3d[["T23_dat"]]$T1)
        colnames(T23_col) <- unique(data_3d[["T23_dat"]]$T2)
        
        for (i in 1:nrow(data_3d[["T23_dat"]])) {
          T23_col[as.character(data_3d[["T23_dat"]][i, 1]), as.character(data_3d[["T23_dat"]][i, 2])] <-
            ifelse(method == "prob",
                   data_3d[["T23_dat"]][i, 3],
                   data_3d[["T23_dat"]][i, 4])
        }
        
        
        
        T23_vals <- T23_col[which(!is.na(T23_col))]
        
        
        ###rows are x which is
        ### columns are y
        
        colours_T23 <-
          data.frame(value = T23_vals[order(T23_vals)], color.name = my_colour_2(length(T23_vals)))
        
        
        T23_col_mat[which(!is.na(T23_col))[order(T23_vals)]] <-
          colours_T23$color.name
        
        
        T23_col_mat <- as.character(T23_col_mat)
      }
    #######################################
    ###################################
    
    
    if (free_limits) {
      limits <- c(
        xmax = max(filled_cells$T2),
        xmin = min(filled_cells$T2),
        ymax = max(filled_cells$T1),
        ymin = min(filled_cells$T1),
        zmax = max(filled_cells$T3),
        zmin = min(filled_cells$T3)
      )
      
      scale <- mean(c(dist(c(
        limits["xmin"], limits["xmax"]
      ))[1],
      dist(c(
        limits["ymin"], limits["xmax"]
      ))[1],
      dist(c(
        limits["zmin"], limits["zmax"]
      ))[1])) / 75
    } else {
      limits <- c(
        xmax = max(data[["data"]][["evaluation_grid"]][[2]]),
        xmin = min(data[["data"]][["evaluation_grid"]][[2]]),
        ymax = max(data[["data"]][["evaluation_grid"]][[1]]),
        ymin = min(data[["data"]][["evaluation_grid"]][[1]]),
        zmax = max(data[["data"]][["evaluation_grid"]][[3]]),
        zmin = min(data[["data"]][["evaluation_grid"]][[3]])
      )
      
      scale <- mean(c(dist(c(
        limits["xmin"], limits["xmax"]
      ))[1],
      dist(c(
        limits["ymin"], limits["xmax"]
      ))[1],
      dist(c(
        limits["zmin"], limits["zmax"]
      ))[1])) / 100
    }
    
    
    
    
      tpd_col <- filled_cells_col[, "color.name"]
    
    
    
    clear3d()
    rgl.viewpoint(theta = 50,
                  phi = 25,
                  zoom = 13 / 16)
    plot3d(
      x,
      y,
      z,
      box = FALSE,
      xlab = "",
      ylab = "",
      zlab = "",
      type = "s",
      radius = scale,
      alpha = 0.8,
      xlim = c(limits["xmin"], limits["xmax"]),
      ylim = c(limits["ymin"], limits["ymax"]),
      zlim = c(limits["zmin"], limits["zmax"]),
      col = tpd_col,
      lwd = 0.1
    )
    
    if (grid) {
      surface3d(
        x = unique(data_3d[["T21_dat"]]$T1),
        y = unique(data_3d[["T21_dat"]]$T2),
        z = matrix(rep(min(data_3d[["pl_dat"]]$T3), 50 * 50),
                   nrow = 50, ncol = 50),
        alpha = 0.5,
        lit = FALSE,
        front = "lines",
        back = "lines"
      )
      surface3d(
        x = unique(data_3d[["T21_dat"]]$T1),
        y = unique(data_3d[["T21_dat"]]$T2),
        z = matrix(rep(min(data_3d[["pl_dat"]]$T3), 50 * 50),
                   nrow = 50, ncol = 50),
        color = T21_col_mat,
        smooth = FALSE,
        lit = FALSE
      )
      surface3d(
        x = rep(min(data_3d[["T21_dat"]]$T1), 50),
        y = unique(data_3d[["T21_dat"]]$T2),
        z = matrix(rep(unique(data_3d[["pl_dat"]]$T3), 50),
                   nrow = 50, ncol = 50),
        alpha = 0.5,
        lit = FALSE,
        front = "lines",
        back = "lines"
      )
      surface3d(
        x = rep(min(data_3d[["T21_dat"]]$T1), 50),
        y = unique(data_3d[["T21_dat"]]$T2),
        z = matrix(rep(unique(data_3d[["pl_dat"]]$T3), 50),
                   nrow = 50, ncol = 50),
        color = T31_col_mat,
        lit = FALSE,
        smooth = FALSE
      )
      surface3d(
        x = unique(data_3d[["T21_dat"]]$T1),
        y = rep(min(data_3d[["T21_dat"]]$T2), 50),
        z = t(matrix(
          rep(unique(data_3d[["pl_dat"]]$T3), 50),
          nrow = 50, ncol = 50
        )),
        alpha = 0.5,
        lit = FALSE,
        front = "lines",
        back = "lines"
      )
      surface3d(
        x = unique(data_3d[["T21_dat"]]$T1),
        y = rep(min(data_3d[["T21_dat"]]$T2), 50),
        z = t(matrix(
          rep(unique(data_3d[["pl_dat"]]$T3), 50),
          nrow = 50, ncol = 50
        )),
        color = T23_col_mat,
        lit = FALSE,
        smooth = FALSE
      )
      
    }
    title3d(
      main = title,
      xlab = T2lab,
      ylab = T1lab,
      zlab = T3lab
    )
    
    if (save) {
      rgl.snapshot(file, fmt = 'png')
    }
    
    
    
  }


# sites_lu <- data.frame(SSBS= names(PREDICTS_tpds)[-1]) %>% dplyr::left_join(PREDICTS[,c("SSBS","Predominant_habitat","Use_intensity")])
#
# primary <- sites_lu %>% dplyr::filter(grepl(Predominant_habitat, pattern = "Primary")) %>% pull(SSBS)
#
#
#
# TPD_3d_plot(sites = sites_lu$SSBS,data = PREDICTS_tpds, T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body", method = "prob")
#
# #####################################################
# #####################################################
# ## SD plot showing the difference between two set of sites
#
# secondary <- sites_lu %>% dplyr::filter(grepl(Predominant_habitat, pattern = "secondary",ignore.case = TRUE)) %>% pull(SSBS)
#
# data <- PREDICTS_tpds
# sites1 <- primary
# sites2 <- secondary
# method <- "prob"

TPD_Diff_plot <-
  function(data,
           sites1,
           sites2,
           T1lab,
           T2lab,
           T3lab,
           method,
           title,
           save = FALSE,
           file) {
    if (!all(sites1 %in% names(data))) {
      sites1 <- sites1[which(sites1 %in% names(data))]
    }
    
    if (!all(sites2 %in% names(data))) {
      sites2 <- sites2[which(sites2 %in% names(data))]
    }
    
    ### extract the TPD data for teh two sets of sites
    
    sites1_data <- TPD_plot_data(data, sites1)
    sites2_data <- TPD_plot_data(data, sites2)
    
    ## for the 3D plot need to just get the cells which are functionally occupied
    
    filled_cells_1 <-
      sites1_data[["pl_dat"]] %>% dplyr::group_by(T2, T1, T3) %>%
      dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% filter(!is.na(prob)) %>% data.frame()
    
    filled_cells_1 <- percentile_cells(filled_cells_1)
    
    filled_cells_2 <-
      sites2_data[["pl_dat"]] %>% dplyr::group_by(T2, T1, T3) %>%
      dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% filter(!is.na(prob)) %>% data.frame()
    
    filled_cells_2 <-
      percentile_cells(filled_cells_2) %>% dplyr::rename(prob_2 = prob, percentile_2 = percentile)
    
    
    if (method == "prob") {
      filled_cells_1 <- filled_cells_1 %>% dplyr::rename(value = prob)
      filled_cells_2 <-
        filled_cells_2 %>% dplyr::rename(value_2 = prob_2)
      my.colors_high <-
        colorRampPalette(c("blue", "steelblue4", "steelblue1", "lightskyblue1"))
      my.colors_low <-
        colorRampPalette(c("tomato", "orangered", "orangered4", "red"))
      my_colour_2p <-
        colorRampPalette(c("blue", "steelblue1", "lightskyblue1"))
      my_colour_2n <- colorRampPalette(c("yellow", "orange", "red"))
    } else {
      filled_cells_1 <-
        filled_cells_1 %>% dplyr::rename(value = percentile)
      filled_cells_2 <-
        filled_cells_2 %>% dplyr::rename(value_2 = percentile_2)
      my.colors_high <-
        colorRampPalette(c("lightskyblue1", "steelblue1", "steelblue4", "blue"))
      my.colors_low <-
        colorRampPalette(c("red", "orangered4", "orangered", "tomato"))
      my_colour_2p <-
        colorRampPalette(c("lightskyblue1", "steelblue1", "blue"))
      my_colour_2n <- colorRampPalette(c("red", "orange", "yellow"))
    }
    
    
    
    ####################################
    # WORK OUT THE DIFFERENCE ##########
    ####################################
    
    ## join the two data frames together with the secondary sites prob and percentile renamed
    
    diff_cells <-
      filled_cells_1 %>% dplyr::left_join(filled_cells_2, by = c("T1", "T2", "T3")) %>% data.frame()
    
    ##### work out with difference in prob/percentile of cells and whether the occupancy of the cell is new (functional gain) or lost (functional loss)
    
    cells_frame <-
      filled_cells_2 %>% dplyr::left_join(filled_cells_1, by = c("T1", "T2", "T3")) %>% dplyr::filter(is.na(value)) %>%
      dplyr::relocate(value, .before = value_2) %>% rbind(diff_cells) %>% dplyr::mutate(
        lost_cell = ifelse(is.na(value_2), TRUE, FALSE),
        gain_cell = ifelse(is.na(value), TRUE, FALSE),
        value = ifelse(is.na(value), 0, value),
        value_2 = ifelse(is.na(value_2), 0, value_2),
        diff = value_2 - value
      )
    
    which(cells_frame$diff > as.numeric(quantile(cells_frame$diff, 0.975)))
    which(cells_frame$diff < as.numeric(quantile(cells_frame$diff, 0.025)))
    #### for percentile the difference is a bit complicated if the difference between sites 1 and two is positive that means that functional occupancy has gone down and vice versa.
    # with increaseing high values means increasing less occupancy and increasing low values means increasing functional occupancy gain
    
    
    ##### HIGH CELLS --- looking at percentil if the difference is less then that means the occupancy goes up ( cells in which occupancy goes up)
    
    
    ## identify cells which have a negative diff
    if (method == "prob") {
      pos_cells <-
        cells_frame %>% dplyr::filter(diff > 0,!gain_cell,!lost_cell) %>% pull(diff)
    } else {
      pos_cells <-
        cells_frame %>% dplyr::filter(diff < 0,!gain_cell,!lost_cell) %>% pull(diff)
    }
    ### blue is going to represent gain here ### Ramp goes from high to low so for negative(positive cells) the scale should be
    
    color.df_high <-
      data.frame(diff = pos_cells[order(pos_cells)], color.name = my.colors_high(length(pos_cells))) %>% distinct(diff, .keep_all = TRUE)#generates 2001 colors from the color ramp
    
    
    ##### LOW CELLS
    
    ## identify cells which have a positive diff(lose occupancy)
    if (method == "prob") {
      neg_cells <-
        cells_frame %>% dplyr::filter(diff < 0,!gain_cell,!lost_cell) %>% pull(diff)
    } else {
      neg_cells <-
        cells_frame %>% dplyr::filter(diff > 0,!gain_cell,!lost_cell) %>% pull(diff)
    }
    ## red is going to represent loss ## ramp is high to low
    
    color.df_low <-
      data.frame(diff = neg_cells[order(neg_cells)], color.name = my.colors_low(length(neg_cells))) %>% distinct(diff, .keep_all = TRUE)#generates 2001 colors from the color ramp
    
    color.df <- rbind(color.df_high, color.df_low)
    
    filled_cells_col <-
      cells_frame %>% dplyr::left_join(color.df, by = "diff") %>%
      dplyr::mutate(color.name = ifelse(is.na(color.name), "#FFFFFF", color.name))
    
    ### If the cells are gained between the two sets of sites then the cell will colour green
    ### if the cells are lost between the two sets of sites then the cell will colour black
    
    filled_cells_col <-
      filled_cells_col %>% dplyr::mutate(color.name = ifelse(gain_cell, "#CCFF99", paste(color.name)))
    filled_cells_col <-
      filled_cells_col %>% dplyr::mutate(color.name = ifelse(lost_cell, "#000000", paste(color.name)))
    
    
    
    ############################################################
    ############################################################
    
    x <- filled_cells_col$T2
    y <- filled_cells_col$T1
    z <- filled_cells_col$T3
    
    
    
    
    #######################################
    ###########################################
    
    
    T21_col_mat <- matrix(rep("#FFFFFF", 50 * 50),
                          nrow = 50, ncol = 50)
    
    T21_col <- matrix(rep(min(sites1_data[["pl_dat"]]$T3), 50 * 50),
                      nrow = 50, ncol = 50)
    rownames(T21_col) <- unique(sites1_data[["T21_dat"]]$T1)
    colnames(T21_col) <- unique(sites1_data[["T21_dat"]]$T2)
    
    s1d <- sites1_data[["T21_dat"]]
    s2d <- sites2_data[["T21_dat"]]
    if (method == "prob") {
      colnames(s1d)[3] <- "value"
      colnames(s2d)[3] <- "value"
      s1d <-
        s1d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
      s2d <-
        s2d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
    } else {
      colnames(s1d)[4] <- "value"
      colnames(s2d)[4] <- "value"
      s1d <-
        s1d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
      s2d <-
        s2d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
    }
    
    
    
    
    
    for (i in 1:nrow(sites1_data[["T21_dat"]])) {
      T21_col[as.character(sites1_data[["T21_dat"]][i, 1]), as.character(sites1_data[["T21_dat"]][i, 2])] <-
        ifelse(method == "prob", ifelse((s2d[i, 3] - s1d[i, 3]) == 0, NA, (s2d[i, 3] - s1d[i, 3])),
               ifelse((s2d[i, 4] - s1d[i, 4]) == 0, NA, (s2d[i, 4] - s1d[i, 4])))
    }
    
    T21_vals <- T21_col[which(!is.na(T21_col))]
    
    
    if (method == "prob") {
      pos_T21_vals <- T21_vals[which(T21_vals > 0)]
    } else {
      pos_T21_vals <- T21_vals[which(T21_vals < 0)]
    }
    
    pos_T21_vals <- pos_T21_vals[order(pos_T21_vals)]
    pos_T21_cols <- my_colour_2p(length(pos_T21_vals))
    
    if (method == "prob") {
      neg_T21_vals <- T21_vals[which(T21_vals < 0)]
    } else {
      neg_T21_vals <- T21_vals[which(T21_vals > 0)]
    }
    
    neg_T21_vals <- neg_T21_vals[order(neg_T21_vals)]
    neg_T21_cols <- my_colour_2n(length(neg_T21_vals))
    
    colours_T21 <-
      data.frame(
        value = c(neg_T21_vals, pos_T21_vals),
        color.name = c(neg_T21_cols, pos_T21_cols)
      )
    
    ###rows are x
    ### columns are y
    
    
    T21_col_mat[which(!is.na(T21_col))[order(T21_vals)]] <-
      colours_T21$color.name
    
    
    T21_col_mat <- as.character(T21_col_mat)
    
    
    #######################################
    ###########################################
    
    
    T31_col_mat <- matrix(rep("#FFFFFF", 50 * 50),
                          nrow = 50, ncol = 50)
    
    
    T31_col <- matrix(rep(min(sites1_data[["pl_dat"]]$T3), 50 * 50),
                      nrow = 50, ncol = 50)
    rownames(T31_col) <- unique(sites1_data[["T31_dat"]]$T1)
    colnames(T31_col) <- unique(sites1_data[["T31_dat"]]$T2)
    
    s1d <- sites1_data[["T31_dat"]]
    s2d <- sites2_data[["T31_dat"]]
    if (method == "prob") {
      colnames(s1d)[3] <- "value"
      colnames(s2d)[3] <- "value"
      s1d <-
        s1d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
      s2d <-
        s2d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
    } else {
      colnames(s1d)[4] <- "value"
      colnames(s2d)[4] <- "value"
      s1d <-
        s1d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
      s2d <-
        s2d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
    }
    
    
    for (i in 1:nrow(sites1_data[["T31_dat"]])) {
      T31_col[as.character(sites1_data[["T31_dat"]][i, 1]), as.character(sites1_data[["T31_dat"]][i, 2])] <-
        ifelse(method == "prob", ifelse((s2d[i, 3] - s1d[i, 3]) == 0, NA, (s2d[i, 3] - s1d[i, 3])),
               ifelse((s2d[i, 4] - s1d[i, 4]) == 0, NA, (s2d[i, 4] - s1d[i, 4])))
    }
    
    T31_vals <- T31_col[which(!is.na(T31_col))]
    
    if (method == "prob") {
      pos_T31_vals <- T31_vals[which(T31_vals > 0)]
    } else {
      pos_T31_vals <- T31_vals[which(T31_vals < 0)]
    }
    
    pos_T31_vals <- pos_T31_vals[order(pos_T31_vals)]
    pos_T31_cols <- my_colour_2p(length(pos_T31_vals))
    
    if (method == "prob") {
      neg_T31_vals <- T31_vals[which(T31_vals < 0)]
    } else {
      neg_T31_vals <- T31_vals[which(T31_vals > 0)]
    }
    
    neg_T31_vals <- neg_T31_vals[order(neg_T31_vals)]
    neg_T31_cols <- my_colour_2n(length(neg_T31_vals))
    
    colours_T31 <-
      data.frame(
        value = c(neg_T31_vals, pos_T31_vals),
        color.name = c(neg_T31_cols, pos_T31_cols)
      )
    
    ###rows are x
    ### columns are y
    
    
    T31_col_mat[which(!is.na(T31_col))[order(T31_vals)]] <-
      colours_T31$color.name
    
    
    T31_col_mat <- as.character(T31_col_mat)
    
    #######################################
    ###########################################
    
    
    T23_col_mat <- matrix(rep("#FFFFFF", 50 * 50),
                          nrow = 50, ncol = 50)
    
    
    T23_col <- matrix(rep(min(sites1_data[["pl_dat"]]$T3), 50 * 50),
                      nrow = 50, ncol = 50)
    rownames(T23_col) <- unique(sites1_data[["T23_dat"]]$T1)
    colnames(T23_col) <- unique(sites1_data[["T23_dat"]]$T2)
    
    s1d <- sites1_data[["T23_dat"]]
    s2d <- sites2_data[["T23_dat"]]
    if (method == "prob") {
      colnames(s1d)[3] <- "value"
      colnames(s2d)[3] <- "value"
      s1d <-
        s1d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
      s2d <-
        s2d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
    } else {
      colnames(s1d)[4] <- "value"
      colnames(s2d)[4] <- "value"
      s1d <-
        s1d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
      s2d <-
        s2d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
    }
    
    for (i in 1:nrow(sites1_data[["T23_dat"]])) {
      T23_col[as.character(sites1_data[["T23_dat"]][i, 1]), as.character(sites1_data[["T23_dat"]][i, 2])] <-
        ifelse(method == "prob", ifelse((s2d[i, 3] - s1d[i, 3]) == 0, NA, (s2d[i, 3] - s1d[i, 3])),
               ifelse((s2d[i, 4] - s1d[i, 4]) == 0, NA, (s2d[i, 4] - s1d[i, 4])))
    }
    
    
    T23_vals <- T23_col[which(!is.na(T23_col))]
    
    if (method == "prob") {
      pos_T23_vals <- T23_vals[which(T23_vals > 0)]
    } else {
      pos_T23_vals <- T23_vals[which(T23_vals < 0)]
    }
    
    pos_T23_vals <- pos_T23_vals[order(pos_T23_vals)]
    pos_T23_cols <- my_colour_2p(length(pos_T23_vals))
    
    if (method == "prob") {
      neg_T23_vals <- T23_vals[which(T23_vals < 0)]
    } else {
      neg_T23_vals <- T23_vals[which(T23_vals > 0)]
    }
    neg_T23_vals <- neg_T23_vals[order(neg_T23_vals)]
    neg_T23_cols <- my_colour_2n(length(neg_T23_vals))
    
    colours_T23 <-
      data.frame(
        value = c(neg_T23_vals, pos_T23_vals),
        color.name = c(neg_T23_cols, pos_T23_cols)
      )
    
    ###rows are x
    ### columns are y
    
    
    T23_col_mat[which(!is.na(T23_col))[order(T23_vals)]] <-
      colours_T23$color.name
    
    
    T23_col_mat <- as.character(T23_col_mat)
    
    
    #######################################
    ###################################
    
    xmax <- max(sites1_data[["pl_dat"]]$T2)
    xmin <- min(sites1_data[["pl_dat"]]$T2)
    ymax <- max(sites1_data[["pl_dat"]]$T1)
    ymin <- min(sites1_data[["pl_dat"]]$T1)
    zmax <- max(sites1_data[["pl_dat"]]$T3)
    zmin <- min(sites1_data[["pl_dat"]]$T3)
    
    scale <- mean(c(dist(c(xmin, xmax))[1],
                    dist(c(ymin, ymax))[1],
                    dist(c(zmin, zmax))[1])) / 100
    
    clear3d()
    rgl.viewpoint(theta = 50,
                  phi = 25,
                  zoom = 13 / 16)
    plot3d(
      x,
      y,
      z,
      box = FALSE,
      xlab = "",
      ylab = "",
      zlab = "",
      type = "s",
      radius = scale,
      alpha = 0.8,
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      zlim = c(zmin, zmax),
      col = filled_cells_col$color.name,
      lwd = 0.1
    )
    surface3d(
      x = unique(sites1_data[["T21_dat"]]$T1),
      y = unique(sites1_data[["T21_dat"]]$T2),
      z = matrix(rep(min(sites1_data[["pl_dat"]]$T3), 50 * 50),
                 nrow = 50, ncol = 50),
      alpha = 0.5,
      lit = FALSE,
      front = "lines",
      back = "lines",
      smooth = FALSE
    )
    surface3d(
      x = unique(sites1_data[["T21_dat"]]$T1),
      y = unique(sites1_data[["T21_dat"]]$T2),
      z = matrix(rep(min(sites1_data[["pl_dat"]]$T3), 50 * 50),
                 nrow = 50, ncol = 50),
      color = T21_col_mat,
      smooth = FALSE,
      lit = FALSE
    )
    surface3d(
      x = rep(min(sites1_data[["T21_dat"]]$T1), 50),
      y = unique(sites1_data[["T21_dat"]]$T2),
      z = matrix(rep(unique(sites1_data[["pl_dat"]]$T3), 50),
                 nrow = 50, ncol = 50),
      alpha = 0.5,
      lit = FALSE,
      front = "lines",
      back = "lines",
      smooth = FALSE
    )
    surface3d(
      x = rep(min(sites1_data[["T21_dat"]]$T1), 50),
      y = unique(sites1_data[["T21_dat"]]$T2),
      z = matrix(rep(unique(sites1_data[["pl_dat"]]$T3), 50),
                 nrow = 50, ncol = 50),
      color = T31_col_mat,
      lit = FALSE,
      smooth = FALSE
    )
    surface3d(
      x = unique(sites1_data[["T21_dat"]]$T1),
      y = rep(min(sites1_data[["T21_dat"]]$T2), 50),
      z = t(matrix(
        rep(unique(sites1_data[["pl_dat"]]$T3), 50),
        nrow = 50, ncol = 50
      )),
      alpha = 0.5,
      lit = FALSE,
      front = "lines",
      back = "lines",
      smooth = FALSE
    )
    surface3d(
      x = unique(sites1_data[["T21_dat"]]$T1),
      y = rep(min(sites1_data[["T21_dat"]]$T2), 50),
      z = t(matrix(
        rep(unique(sites1_data[["pl_dat"]]$T3), 50),
        nrow = 50, ncol = 50
      )),
      color = T23_col_mat,
      lit = FALSE,
      smooth = FALSE
    )
    title3d(
      main = title,
      xlab = T2lab,
      ylab = T1lab,
      zlab = T3lab
    )
    
    if (save) {
      rgl.snapshot(file, fmt = 'png')
    }
    
  }

#TPD_Diff_plot(sites1 = primary,sites2 = secondary,data = PREDICTS_tpds, T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body")


TPD_ranDiff_plot <-
  function(data,
           randata,
           sites,
           T1lab = "",
           T2lab = "",
           T3lab = "",
           method,
           threshold,
           title = "",
           save = FALSE,
           file,
           grid = FALSE,
           observed = FALSE,
           free_limits = FALSE) {
    ### extract the TPD data for teh two sets of sites
    
    if (!all(sites %in% names(data)) | !all(sites %in% names(randata))) {
      sites <- sites[which(sites %in% names(data))]
      sites <- sites[which(sites %in% names(randata))]
    }
    
    sites1_data <- TPD_plot_data(data, sites)
    sites2_data <- TPD_plot_data(randata, sites)
    
    ## for the 3D plot need to just get the cells which are functionally occupied
    
    filled_cells_1 <-
      sites1_data[["pl_dat"]] %>% dplyr::group_by(T2, T1, T3) %>%
      dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% filter(!is.na(prob)) %>% data.frame()
    
    filled_cells_1 <- percentile_cells(filled_cells_1)
    
    filled_cells_2 <-
      sites2_data[["pl_dat"]] %>% dplyr::group_by(T2, T1, T3) %>%
      dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% filter(!is.na(prob)) %>% data.frame()
    
    filled_cells_2 <-
      percentile_cells(filled_cells_2) %>% dplyr::rename(prob_2 = prob, percentile_2 = percentile) %>%
      dplyr::filter(percentile_2 <= threshold)
    
    
    ####################################################
    ## Construct hull
    
    
    obs_cvh <- convhulln(filled_cells_1[, c("T1", "T2", "T3")])
    
    filled_cells_2 <-
      filled_cells_2 %>% dplyr::mutate(in_hull = inhulln(obs_cvh, as.matrix(filled_cells_2[, c("T1", "T2", "T3")]))) %>%
      dplyr::filter(in_hull)
    
    
    
    if (method == "prob") {
      filled_cells_1 <- filled_cells_1 %>% dplyr::rename(value = prob)
      filled_cells_2 <-
        filled_cells_2 %>% dplyr::rename(value_2 = prob_2)
      my.colors_high <-
        colorRampPalette(c("blue", "steelblue4", "steelblue1", "lightskyblue1"))
      my.colors_low <-
        colorRampPalette(c("tomato", "orangered", "orangered4", "red"))
      my_colour_2p <-
        colorRampPalette(c("blue", "steelblue1", "lightskyblue1"))
      my_colour_2n <- colorRampPalette(c("yellow", "orange", "red"))
    } else {
      filled_cells_1 <-
        filled_cells_1 %>% dplyr::rename(value = percentile)
      filled_cells_2 <-
        filled_cells_2 %>% dplyr::rename(value_2 = percentile_2)
      my.colors_high <-
        colorRampPalette(c("lightskyblue1", "steelblue1", "steelblue4", "blue"))
      my.colors_low <-
        colorRampPalette(c("red", "orangered4", "orangered", "tomato"))
      my_colour_2p <-
        colorRampPalette(c("lightskyblue1", "steelblue1", "blue"))
      my_colour_2n <- colorRampPalette(c("red", "orange", "yellow"))
    }
    
    
    
    ####################################
    # WORK OUT THE DIFFERENCE ##########
    ####################################
    
    ## join the two data frames together with the secondary sites prob and percentile renamed
    
    diff_cells <-
      filled_cells_1 %>% merge(filled_cells_2, by = c("T1", "T2", "T3")) %>% data.frame()
    
    ##### work out with difference in prob/percentile of cells and whether the occupancy of the cell is new (functional gain) or lost (functional loss)
    
    cells_frame <-
      filled_cells_2 %>% dplyr::left_join(filled_cells_1, by = c("T1", "T2", "T3")) %>% dplyr::filter(is.na(value)) %>%
      dplyr::relocate(value, .before = value_2) %>% rbind(diff_cells) %>% dplyr::mutate(
        lost_cell = ifelse(is.na(value_2), TRUE, FALSE),
        gain_cell = ifelse(is.na(value), TRUE, FALSE),
        value = ifelse(is.na(value), 0, value),
        value_2 = ifelse(is.na(value_2), 0, value_2),
        diff = value_2 - value
      )
    
    
    if (!observed) {
      cells_frame <- cells_frame %>% dplyr::filter(value == 0)
    }
    
    #### for percentile the difference is a bit complicated if the difference between sites 1 and two is positive that means that functional occupancy has gone down and vice versa.
    # with increaseing high values means increasing less occupancy and increasing low values means increasing functional occupancy gain
    
    
    ##### HIGH CELLS --- looking at percentil if the difference is less then that means the occupancy goes up ( cells in which occupancy goes up)
    
    
    ## identify cells which have a negative diff
    if (method == "prob") {
      pos_cells <-
        cells_frame %>% dplyr::filter(diff > 0,!gain_cell,!lost_cell) %>% pull(diff)
    } else {
      pos_cells <-
        cells_frame %>% dplyr::filter(diff < 0,!gain_cell,!lost_cell) %>% pull(diff)
    }
    ### blue is going to represent gain here ### Ramp goes from high to low so for negative(positive cells) the scale should be
    
    color.df_high <-
      data.frame(diff = pos_cells[order(pos_cells)], color.name = my.colors_high(length(pos_cells))) %>% distinct(diff, .keep_all = TRUE)#generates 2001 colors from the color ramp
    
    
    ##### LOW CELLS
    
    ## identify cells which have a positive diff(lose occupancy)
    if (method == "prob") {
      neg_cells <-
        cells_frame %>% dplyr::filter(diff < 0,!gain_cell,!lost_cell) %>% pull(diff)
    } else {
      neg_cells <-
        cells_frame %>% dplyr::filter(diff > 0,!gain_cell,!lost_cell) %>% pull(diff)
    }
    ## red is going to represent loss ## ramp is high to low
    
    color.df_low <-
      data.frame(diff = neg_cells[order(neg_cells)], color.name = my.colors_low(length(neg_cells))) %>% distinct(diff, .keep_all = TRUE)#generates 2001 colors from the color ramp
    
    color.df <- rbind(color.df_high, color.df_low)
    
    filled_cells_col <-
      cells_frame %>% dplyr::left_join(color.df, by = "diff") %>%
      dplyr::mutate(color.name = ifelse(is.na(color.name), "#FFFFFF", color.name))
    
    ### If the cells are gained between the two sets of sites then the cell will colour green
    ### if the cells are lost between the two sets of sites then the cell will colour black
    
    filled_cells_col <-
      filled_cells_col %>% dplyr::mutate(color.name = ifelse(gain_cell, "#990000", paste(color.name)))
    filled_cells_col <-
      filled_cells_col %>% dplyr::mutate(color.name = ifelse(lost_cell, "#000000", paste(color.name)))
    
    
    
    ############################################################
    ############################################################
    
    x <- filled_cells_col$T2
    y <- filled_cells_col$T1
    z <- filled_cells_col$T3
    
    
    
    #######################################
    ###########################################
    if (grid) {
      T21_col_mat <- matrix(rep("#FFFFFF", 50 * 50),
                            nrow = 50, ncol = 50)
      
      T21_col <- matrix(rep(min(sites1_data[["pl_dat"]]$T3), 50 * 50),
                        nrow = 50, ncol = 50)
      rownames(T21_col) <- unique(sites1_data[["T21_dat"]]$T1)
      colnames(T21_col) <- unique(sites1_data[["T21_dat"]]$T2)
      
      s1d <- sites1_data[["T21_dat"]]
      s2d <- sites2_data[["T21_dat"]]
      
      
      if (method == "prob") {
        colnames(s1d)[3] <- "value"
        colnames(s2d)[3] <- "value"
        s1d <-
          s1d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
        s2d <-
          s2d %>% dplyr::mutate(value = ifelse(is.na(value) |
                                                 percentile >= threshold, 0, value))
      } else {
        colnames(s1d)[4] <- "value"
        colnames(s2d)[4] <- "value"
        s1d <-
          s1d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
        s2d <-
          s2d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
      }
      
      
      
      
      
      for (i in 1:nrow(sites1_data[["T21_dat"]])) {
        T21_col[as.character(sites1_data[["T21_dat"]][i, 1]), as.character(sites1_data[["T21_dat"]][i, 2])] <-
          ifelse(method == "prob", ifelse((s2d[i, 3] - s1d[i, 3]) == 0, NA, (s2d[i, 3] - s1d[i, 3])),
                 ifelse((s2d[i, 4] - s1d[i, 4]) == 0, NA, (s2d[i, 4] - s1d[i, 4])))
      }
      
      T21_vals <- T21_col[which(!is.na(T21_col))]
      
      
      if (method == "prob") {
        pos_T21_vals <- T21_vals[which(T21_vals > 0)]
      } else {
        pos_T21_vals <- T21_vals[which(T21_vals < 0)]
      }
      
      pos_T21_vals <- pos_T21_vals[order(pos_T21_vals)]
      pos_T21_cols <- my_colour_2p(length(pos_T21_vals))
      
      if (method == "prob") {
        neg_T21_vals <- T21_vals[which(T21_vals < 0)]
      } else {
        neg_T21_vals <- T21_vals[which(T21_vals > 0)]
      }
      
      neg_T21_vals <- neg_T21_vals[order(neg_T21_vals)]
      neg_T21_cols <- my_colour_2n(length(neg_T21_vals))
      
      colours_T21 <-
        data.frame(
          value = c(neg_T21_vals, pos_T21_vals),
          color.name = c(neg_T21_cols, pos_T21_cols)
        )
      
      ###rows are x
      ### columns are y
      
      
      T21_col_mat[which(!is.na(T21_col))[order(T21_vals)]] <-
        colours_T21$color.name
      
      
      T21_col_mat <- as.character(T21_col_mat)
      
      
      #######################################
      ###########################################
      
      
      T31_col_mat <- matrix(rep("#FFFFFF", 50 * 50),
                            nrow = 50, ncol = 50)
      
      T31_col <- matrix(rep(min(sites1_data[["pl_dat"]]$T3), 50 * 50),
                        nrow = 50, ncol = 50)
      rownames(T31_col) <- unique(sites1_data[["T31_dat"]]$T1)
      colnames(T31_col) <- unique(sites1_data[["T31_dat"]]$T2)
      
      s1d <- sites1_data[["T31_dat"]]
      s2d <- sites2_data[["T31_dat"]]
      
      
      if (method == "prob") {
        colnames(s1d)[3] <- "value"
        colnames(s2d)[3] <- "value"
        s1d <-
          s1d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
        s2d <-
          s2d %>% dplyr::mutate(value = ifelse(is.na(value) |
                                                 percentile >= threshold, 0, value))
      } else {
        colnames(s1d)[4] <- "value"
        colnames(s2d)[4] <- "value"
        s1d <-
          s1d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
        s2d <-
          s2d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
      }
      
      
      
      
      
      for (i in 1:nrow(sites1_data[["T31_dat"]])) {
        T31_col[as.character(sites1_data[["T31_dat"]][i, 1]), as.character(sites1_data[["T31_dat"]][i, 2])] <-
          ifelse(method == "prob", ifelse((s2d[i, 3] - s1d[i, 3]) == 0, NA, (s2d[i, 3] - s1d[i, 3])),
                 ifelse((s2d[i, 4] - s1d[i, 4]) == 0, NA, (s2d[i, 4] - s1d[i, 4])))
      }
      
      T31_vals <- T31_col[which(!is.na(T31_col))]
      
      
      if (method == "prob") {
        pos_T31_vals <- T31_vals[which(T31_vals > 0)]
      } else {
        pos_T31_vals <- T31_vals[which(T31_vals < 0)]
      }
      
      pos_T31_vals <- pos_T31_vals[order(pos_T31_vals)]
      pos_T31_cols <- my_colour_2p(length(pos_T31_vals))
      
      if (method == "prob") {
        neg_T31_vals <- T31_vals[which(T31_vals < 0)]
      } else {
        neg_T31_vals <- T31_vals[which(T31_vals > 0)]
      }
      
      neg_T31_vals <- neg_T31_vals[order(neg_T31_vals)]
      neg_T31_cols <- my_colour_2n(length(neg_T31_vals))
      
      colours_T31 <-
        data.frame(
          value = c(neg_T31_vals, pos_T31_vals),
          color.name = c(neg_T31_cols, pos_T31_cols)
        )
      
      ###rows are x
      ### columns are y
      
      
      T31_col_mat[which(!is.na(T31_col))[order(T31_vals)]] <-
        colours_T31$color.name
      
      
      T31_col_mat <- as.character(T31_col_mat)
      
      #######################################
      ###########################################
      
      
      T23_col_mat <- matrix(rep("#FFFFFF", 50 * 50),
                            nrow = 50, ncol = 50)
      
      T23_col <- matrix(rep(min(sites1_data[["pl_dat"]]$T3), 50 * 50),
                        nrow = 50, ncol = 50)
      rownames(T23_col) <- unique(sites1_data[["T23_dat"]]$T1)
      colnames(T23_col) <- unique(sites1_data[["T23_dat"]]$T2)
      
      s1d <- sites1_data[["T23_dat"]]
      s2d <- sites2_data[["T23_dat"]]
      
      
      if (method == "prob") {
        colnames(s1d)[3] <- "value"
        colnames(s2d)[3] <- "value"
        s1d <-
          s1d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
        s2d <-
          s2d %>% dplyr::mutate(value = ifelse(is.na(value) |
                                                 percentile >= threshold, 0, value))
      } else {
        colnames(s1d)[4] <- "value"
        colnames(s2d)[4] <- "value"
        s1d <-
          s1d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
        s2d <-
          s2d %>% dplyr::mutate(value = ifelse(is.na(value), 0, value))
      }
      
      
      
      
      
      for (i in 1:nrow(sites1_data[["T23_dat"]])) {
        T23_col[as.character(sites1_data[["T23_dat"]][i, 1]), as.character(sites1_data[["T23_dat"]][i, 2])] <-
          ifelse(method == "prob", ifelse((s2d[i, 3] - s1d[i, 3]) == 0, NA, (s2d[i, 3] - s1d[i, 3])),
                 ifelse((s2d[i, 4] - s1d[i, 4]) == 0, NA, (s2d[i, 4] - s1d[i, 4])))
      }
      
      T23_vals <- T23_col[which(!is.na(T23_col))]
      
      
      if (method == "prob") {
        pos_T23_vals <- T23_vals[which(T23_vals > 0)]
      } else {
        pos_T23_vals <- T23_vals[which(T23_vals < 0)]
      }
      
      pos_T23_vals <- pos_T23_vals[order(pos_T23_vals)]
      pos_T23_cols <- my_colour_2p(length(pos_T23_vals))
      
      if (method == "prob") {
        neg_T23_vals <- T23_vals[which(T23_vals < 0)]
      } else {
        neg_T23_vals <- T23_vals[which(T23_vals > 0)]
      }
      
      neg_T23_vals <- neg_T23_vals[order(neg_T23_vals)]
      neg_T23_cols <- my_colour_2n(length(neg_T23_vals))
      
      colours_T23 <-
        data.frame(
          value = c(neg_T23_vals, pos_T23_vals),
          color.name = c(neg_T23_cols, pos_T23_cols)
        )
      
      ###rows are x
      ### columns are y
      
      
      T23_col_mat[which(!is.na(T23_col))[order(T23_vals)]] <-
        colours_T23$color.name
      
      
      T23_col_mat <- as.character(T23_col_mat)
      
    }
    #######################################
    ###################################
    
    if (free_limits) {
      xmax = max(filled_cells$T2)
      xmin = min(filled_cells_1$T2)
      ymax = max(filled_cells_1$T1)
      ymin = min(filled_cells_1$T1)
      zmax = max(filled_cells_1$T3)
      zmin = min(filled_cells_1$T3)
      
      scale <- mean(c(dist(c(xmin, xmax))[1],
                      dist(c(ymin, ymax))[1],
                      dist(c(zmin, zmax))[1])) / 75
      
      
    } else {
      xmax <- max(sites1_data[["pl_dat"]]$T2)
      xmin <- min(sites1_data[["pl_dat"]]$T2)
      ymax <- max(sites1_data[["pl_dat"]]$T1)
      ymin <- min(sites1_data[["pl_dat"]]$T1)
      zmax <- max(sites1_data[["pl_dat"]]$T3)
      zmin <- min(sites1_data[["pl_dat"]]$T3)
      
      scale <- mean(c(dist(c(xmin, xmax))[1],
                      dist(c(ymin, ymax))[1],
                      dist(c(zmin, zmax))[1])) / 100
      
    }
    
    
    clear3d()
    rgl.viewpoint(theta = 50,
                  phi = 25,
                  zoom = 13 / 16)
    plot3d(
      x,
      y,
      z,
      box = FALSE,
      xlab = "",
      ylab = "",
      zlab = "",
      type = "s",
      radius = scale,
      alpha = 0.8,
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      zlim = c(zmin, zmax),
      col = filled_cells_col$color.name,
      lwd = 0.1
    )
    
    if (grid) {
      surface3d(
        x = unique(sites1_data[["T21_dat"]]$T1),
        y = unique(sites1_data[["T21_dat"]]$T2),
        z = matrix(rep(min(
          sites1_data[["pl_dat"]]$T3
        ), 50 * 50),
        nrow = 50, ncol = 50),
        alpha = 0.5,
        lit = FALSE,
        front = "lines",
        back = "lines",
        smooth = FALSE
      )
      surface3d(
        x = unique(sites1_data[["T21_dat"]]$T1),
        y = unique(sites1_data[["T21_dat"]]$T2),
        z = matrix(rep(min(
          sites1_data[["pl_dat"]]$T3
        ), 50 * 50),
        nrow = 50, ncol = 50),
        color = T21_col_mat,
        smooth = FALSE,
        lit = FALSE
      )
      surface3d(
        x = rep(min(sites1_data[["T21_dat"]]$T1), 50),
        y = unique(sites1_data[["T21_dat"]]$T2),
        z = matrix(rep(unique(
          sites1_data[["pl_dat"]]$T3
        ), 50),
        nrow = 50, ncol = 50),
        alpha = 0.5,
        lit = FALSE,
        front = "lines",
        back = "lines",
        smooth = FALSE
      )
      surface3d(
        x = rep(min(sites1_data[["T21_dat"]]$T1), 50),
        y = unique(sites1_data[["T21_dat"]]$T2),
        z = matrix(rep(unique(
          sites1_data[["pl_dat"]]$T3
        ), 50),
        nrow = 50, ncol = 50),
        color = T31_col_mat,
        lit = FALSE,
        smooth = FALSE
      )
      surface3d(
        x = unique(sites1_data[["T21_dat"]]$T1),
        y = rep(min(sites1_data[["T21_dat"]]$T2), 50),
        z = t(matrix(
          rep(unique(sites1_data[["pl_dat"]]$T3), 50),
          nrow = 50, ncol = 50
        )),
        alpha = 0.5,
        lit = FALSE,
        front = "lines",
        back = "lines",
        smooth = FALSE
      )
      surface3d(
        x = unique(sites1_data[["T21_dat"]]$T1),
        y = rep(min(sites1_data[["T21_dat"]]$T2), 50),
        z = t(matrix(
          rep(unique(sites1_data[["pl_dat"]]$T3), 50),
          nrow = 50, ncol = 50
        )),
        color = T23_col_mat,
        lit = FALSE,
        smooth = FALSE
      )
    }
    title3d(
      main = title,
      xlab = T2lab,
      ylab = T1lab,
      zlab = T3lab
    )
    
    if (save) {
      rgl.snapshot(file, fmt = 'png')
    }
    
    
  }


##############################################################
##############################################################
##############################################################



Calc_FRich <-
  function(data,
           sites = NULL,
           sites1 = NULL,
           sites2 = NULL) {
    
    
    if (!is.null(data[["data"]])) {
      cell_volume <- data[["data"]][["cell_volume"]]
    }
    
    if (!is.null(sites1) & !is.null(sites2)) {
      sites1_data <- TPD_plot_data(data, sites1)
      sites2_data <- TPD_plot_data(data, sites2)
      
      site_data <- data[c(sites1, sites2)]
      site_data$sites1$TPDc$RelativeAbundance <-
        sites1_data[["pl_dat"]][["prob"]]
      site_data$sites2$TPDc$RelativeAbundance <-
        sites2_data[["pl_dat"]][["prob"]]
    } else {
      site_data <- data[c(sites)]
      if (!is.null(data[["data"]])) {
        sites_data <- TPD_plot_data(data, sites)
        site_data$total_sites$TPDc$RelativeAbundance <-
          sites_data[["pl_dat"]][["prob"]]
      }
      
    }
    
    
    FRich <- data.frame(SSBS = names(site_data), FRich = NA)
    j <- 1
    for (i in 1:length(site_data)) {
      if (is.null(data[["data"]])) {
        cell_volume <- data[[i]][["data"]][["cell_volume"]]
      }
      
      
      if (!is.null(site_data[[i]][["TPDc"]])) {
        TPD_aux <- site_data[[i]][["TPDc"]][["RelativeAbundance"]]
      } else {
        TPD_aux <- site_data[[i]][["obsTPDc"]][["RelativeAbundance"]]
      }
      TPD_aux[TPD_aux > 0] <- cell_volume
      FRich[j, "FRich"] <- sum(TPD_aux)
      j <- j + 1
    }
    
    return(FRich)
  }



################################################
####### Functional Evenness calc
Calc_FEve <-
  function(data,
           sites = NULL,
           sites1 = NULL,
           sites2 = NULL) {
    if (!is.null(sites1) & !is.null(sites2)) {
      sites1_data <- TPD_plot_data(data, sites1)
      sites2_data <- TPD_plot_data(data, sites2)
      
      site_data <- data[c(sites1, sites2)]
      site_data$sites1$TPDc$RelativeAbundance <-
        sites1_data[["pl_dat"]][["prob"]]
      site_data$sites2$TPDc$RelativeAbundance <-
        sites2_data[["pl_dat"]][["prob"]]
    } else {
      site_data <- data[c(sites)]
      if (!is.null(data[["data"]])) {
        sites_data <- TPD_plot_data(data, sites)
        site_data$total_sites$TPDc$RelativeAbundance <-
          sites_data[["pl_dat"]][["prob"]]
      }
    }
    
    j <- 1
    FEve <- data.frame(SSBS = names(site_data), FEve = NA)
    for (i in 1:length(site_data)) {
      if (!is.null(site_data[[i]][["TPDc"]])) {
        TPD <- site_data[[i]][["TPDc"]][["RelativeAbundance"]]
      } else {
        TPD <- site_data[[i]][["obsTPDc"]][["RelativeAbundance"]]
      }
      TPD_aux <- TPD[TPD > 0]
      TPD_eve <- rep((1 / length(TPD_aux)), times = length(TPD_aux))
      FEve[j, "FEve"] <- sum(pmin(TPD_aux, TPD_eve))
      j <- j + 1
    }
    return(FEve)
  }


##############################################################
#############################################################
# Functional Divergence



Calc_FDiv <-
  function(data,
           sites = NULL,
           sites1 = NULL,
           sites2 = NULL) {
    if (!is.null(data[["data"]])) {
      evaluation_grid <- data$data$evaluation_grid
      cell_volume <- data$data$cell_volume
    }
    
    
    if (!is.null(sites1) & !is.null(sites2)) {
      sites1_data <- TPD_plot_data(data, sites1)
      sites2_data <- TPD_plot_data(data, sites2)
      
      site_data <- data[c(sites1, sites2)]
      site_data$sites1$TPDc$RelativeAbundance <-
        sites1_data[["pl_dat"]][["prob"]]
      site_data$sites2$TPDc$RelativeAbundance <-
        sites2_data[["pl_dat"]][["prob"]]
    } else {
      site_data <- data[c(sites)]
      
      if (!is.null(data[["data"]])) {
        sites_data <- TPD_plot_data(data, sites)
        site_data$total_sites$TPDc$RelativeAbundance <-
          sites_data[["pl_dat"]][["prob"]]
      }
    }
    
    FDiv <- data.frame(SSBS = names(site_data), FDiv = NA)
    k <- 1
    for (i in 1:length(site_data)) {
      if (is.null(data[["data"]])) {
        cell_volume <- site_data[[i]][["data"]][["cell_volume"]]
        evaluation_grid <-
          site_data[[i]][["data"]][["evaluation_grid"]]
      }
      
      if (!is.null(site_data[[i]][["TPDc"]])) {
        TPD <- site_data[[i]][["TPDc"]][["RelativeAbundance"]]
      } else {
        TPD <- site_data[[i]][["obsTPDc"]][["RelativeAbundance"]]
      }
      
      functional_volume <- evaluation_grid[TPD > 0,
                                           , drop = F]
      TPD_aux <- TPD[TPD > 0]
      
      
      for (j in 1:ncol(functional_volume)) {
        # relative distances along each axes from the minimum value
        functional_volume[, j] <-
          (functional_volume[, j] - min(functional_volume[, j])) / (max(functional_volume[, j]) - min(functional_volume[, j]))
      }
      
      
      ### Centre of gravity of all the points
      
      COG <- colMeans(functional_volume, na.rm = TRUE)
      
      
      ## calculate teh distances between the coordinates and the centre of gravity of the hypervolume
      dist_COG <- function(x, COG) {
        result_aux <- stats::dist(rbind(x, COG))
        return(result_aux)
      }
      ### get the distances
      COGDist <- apply(functional_volume, 1, dist_COG,
                       COG)
      
      ## get the mean distance between points and centre of gravity
      meanCOGDist <- mean(COGDist)
      ### deviances away from the mean centre of gravity
      distDeviances <- COGDist - meanCOGDist
      ##abundance weighted deviances
      AWdistDeviances <- sum(TPD_aux * distDeviances)
      ## absolute deviances
      absdistDeviances <- abs(COGDist - meanCOGDist)
      ## abundance weighted absolute distances from the centre of gravity
      AWabsdistDeviances <- sum(TPD_aux * absdistDeviances)
      ## Functional Divergence equals
      FDiv[k, "FDiv"] <-
        (AWdistDeviances + meanCOGDist) / (AWabsdistDeviances +
                                             meanCOGDist)
      k <- k + 1
    }
    return(FDiv)
  }

##################################
###### Functional Dispersion #####
##################################


Calc_FDis <-
  function(data,
           sites = NULL,
           sites1 = NULL,
           sites2 = NULL) {
    if (!is.null(data[["data"]])) {
      evaluation_grid <- data$data$evaluation_grid
      cell_volume <- data$data$cell_volume
    }
    
    if (!is.null(sites1) & !is.null(sites2)) {
      sites1_data <- TPD_plot_data(data, sites1)
      sites2_data <- TPD_plot_data(data, sites2)
      
      site_data <- data[c(sites1, sites2)]
      site_data$sites1$TPDc$RelativeAbundance <-
        sites1_data[["pl_dat"]][["prob"]]
      site_data$sites2$TPDc$RelativeAbundance <-
        sites2_data[["pl_dat"]][["prob"]]
    } else {
      site_data <- data[c(sites)]
      
      if (!is.null(data[["data"]])) {
        sites_data <- TPD_plot_data(data, sites)
        site_data$total_sites$TPDc$RelativeAbundance <-
          sites_data[["pl_dat"]][["prob"]]
      }
    }
    
    FDis <- data.frame(SSBS = names(site_data), FDis = NA)
    k <- 1
    for (i in 1:length(site_data)) {
      if (is.null(data[["data"]])) {
        cell_volume <- site_data[[i]][["data"]][["cell_volume"]]
        evaluation_grid <-
          site_data[[i]][["data"]][["evaluation_grid"]]
      }
      
      if (!is.null(site_data[[i]][["TPDc"]])) {
        TPD <- site_data[[i]][["TPDc"]][["RelativeAbundance"]]
      } else {
        TPD <- site_data[[i]][["obsTPDc"]][["RelativeAbundance"]]
      }
      
      functional_volume <- evaluation_grid[TPD > 0,
                                           , drop = F]
      TPD_aux <- TPD[TPD > 0]
      
      
      for (j in 1:ncol(functional_volume)) {
        # relative distances along each axes from the minimum value
        functional_volume[, j] <-
          (functional_volume[, j] - min(functional_volume[, j])) / (max(functional_volume[, j]) - min(functional_volume[, j]))
      }
      
      
      ### weighted centre point of all the points -- this is just the weighted mean but it is what we have here.
      
      weighted_COG_fun <- function(x, w) {
        results <- c()
        for (weight in 1:length(x)) {
          results <- c(results, x[weight] * w[weight])
        }
        results <- sum(results)
        return(results)
      }
      
      
      weighted_COG <-
        apply(functional_volume, 2, weighted_COG_fun, TPD_aux) / sum(TPD_aux)
      
      
      ## calculate teh distances between the coordinates and the centre of gravity of the hypervolume
      dist_COG <- function(x, COG) {
        result_aux <- stats::dist(rbind(x, COG))
        return(result_aux)
      }
      
      
      ### get the distances
      COGDist <- as.matrix(c(apply(
        functional_volume, 1, dist_COG,
        weighted_COG
      )))
      
      weighted_dist <-  t(as.matrix(TPD_aux)) %*% COGDist
      
      
      FDis[k, "FDis"] <- as.numeric(weighted_dist) / sum(TPD_aux)
      k <- k + 1
    }
    return(FDis)
  }


#############
#############
############# Measure of "Roundness"

Calc_roundness <-
  function(data,
           sites = NULL,
           sites1 = NULL,
           sites2 = NULL) {
    if (!is.null(data[["data"]])) {
      evaluation_grid <- data$data$evaluation_grid
    }
    
    if (!is.null(sites1) & !is.null(sites2)) {
      sites1_data <- TPD_plot_data(data, sites1)
      sites2_data <- TPD_plot_data(data, sites2)
      
      site_data <- data[c(sites1, sites2)]
      site_data$sites1$TPDc$RelativeAbundance <-
        sites1_data[["pl_dat"]][["prob"]]
      site_data$sites2$TPDc$RelativeAbundance <-
        sites2_data[["pl_dat"]][["prob"]]
    } else {
      site_data <- data[c(sites)]
      
      if (!is.null(data[["data"]])) {
        sites_data <- TPD_plot_data(data, sites)
        site_data$total_sites$TPDc$RelativeAbundance <-
          sites_data[["pl_dat"]][["prob"]]
      }
    }
    
    FRound <- data.frame(SSBS = names(site_data), FRound = NA)
    k <- 1
    for (i in 1:length(site_data)) {
      if (is.null(data[["data"]])) {
        evaluation_grid <- site_data[[i]][["data"]][["evaluation_grid"]]
      }
      
      if (!is.null(site_data[[i]][["TPDc"]])) {
        TPD <- site_data[[i]][["TPDc"]][["RelativeAbundance"]]
      } else {
        TPD <- site_data[[i]][["obsTPDc"]][["RelativeAbundance"]]
      }
      functional_volume <- evaluation_grid[TPD > 0,
                                           , drop = F]
      for (j in 1:ncol(functional_volume)) {
        # relative distances along each axes from the minimum value
        functional_volume[, j] <-
          (functional_volume[, j] - min(functional_volume[, j])) / (max(functional_volume[, j]) - min(functional_volume[, j]))
      }
      TPD_aux <- TPD[TPD > 0]
      
      ### centre point of all the points
      COG <- colMeans(functional_volume, na.rm = T)
      
      ## calculate teh distances between the coordinates and the centre of gravity of the hypervolume
      dist_COG <- function(x, COG) {
        result_aux <- stats::dist(rbind(x, COG))
        return(result_aux)
      }
      
      
      if (nrow(functional_volume) < 4 |
          any(is.nan(as.matrix(functional_volume)))) {
        FRound[k, "FRound"] <- NA
        next()
      }
      
      obs_cvh <- t(convhulln(functional_volume))
      
      hull_points <- cbind(functional_volume[obs_cvh, 1],
                           functional_volume[obs_cvh, 2],
                           functional_volume[obs_cvh, 3])
      
      
      
      
      COGDist <- apply(hull_points, 1, dist_COG, COG)
      
      
      ### get the distances
      COGDist <- apply(functional_volume, 1, dist_COG,
                       COG)
      
      
      
      roundness <- (sd(COGDist) / mean(COGDist)) * 100
      
      FRound[k, "FRound"] <- roundness
      
      k <- k + 1
    }
    return(FRound)
  }

### Combine all the functional trait space metrics together

TPD_FD_metrics <-
  function(data,
           sites = NULL,
           sites1 = NULL,
           sites2 = NULL) {
    if (!is.null(sites1) & !is.null(sites2)) {
      FRich <- Calc_FRich(data, sites1, sites2)
      FEve <- Calc_FEve(data, sites1, sites2)
      FDiv <- Calc_FDiv(data, sites1, sites2)
      FRound <- Calc_roundness(data, sites1, sites2)
      FDis <- Calc_FDis(data, sites1, sites2)
    } else {
      FRich <- Calc_FRich(data, sites)
      FEve <- Calc_FEve(data, sites)
      FDiv <- Calc_FDiv(data, sites)
      FRound <- Calc_roundness(data, sites)
      FDis <- Calc_FDis(data, sites)
    }
    
    
    FD_metrics <- FRich %>% dplyr::left_join(FEve, by = "SSBS") %>%
      dplyr::left_join(FDiv, by = "SSBS") %>%
      dplyr::left_join(FRound, by = "SSBS") %>%
      dplyr::left_join(FDis, by = "SSBS") %>%
      set_rownames(FRich$SSBS) %>% dplyr::select(-SSBS)
    
    return(FD_metrics)
  }



############################################
############################################
obj_2_string <- function(x) {
  str <- deparse(substitute(x))
  return(str)
}


Calc_dissim <- function(data, sites1, sites2) {
  results_samp <- list()
  
  sites1_data <- TPD_plot_data(data, sites1)
  sites2_data <- TPD_plot_data(data, sites2)
  
  
  TPD_i <- sites1_data[["pl_dat"]][["prob"]]
  TPD_j <- sites2_data[["pl_dat"]][["prob"]]
  
  
  ## sum of all the minimum
  O_aux <- sum(pmin(TPD_i, TPD_j))
  
  
  shared_aux <- which(TPD_i > 0 & TPD_j > 0)
  A_aux <- sum(pmax(TPD_i[shared_aux], TPD_j[shared_aux])) -
    O_aux
  only_in_i_aux <- which(TPD_i > 0 & TPD_j ==
                           0)
  B_aux <- sum(TPD_i[only_in_i_aux])
  only_in_j_aux <- which(TPD_i == 0 & TPD_j >
                           0)
  C_aux <- sum(TPD_j[only_in_j_aux])
  results_samp$dissim$dissimilarity <- 1 - O_aux
  results_samp$dissim$P_non_shared     <-
    (2 * min(B_aux, C_aux)) / (A_aux + 2 * min(B_aux, C_aux))
  results_samp$dissim$P_shared <-
    1 - results_samp$dissim$P_non_shared
  
  
  return(results_samp)
  
}

Calc_dissim_random <- function(data, randata, sites, threshold) {
  if (!all(sites %in% names(data)) | !all(sites %in% names(randata))) {
    sites <- sites[which(sites %in% names(data))]
    sites <- sites[which(sites %in% names(randata))]
  }
  
  
  results_samp <- list()
  results_samp$dissim$dissimilarity <- NA
  results_samp$dissim$P_shared <- NA
  results_samp$dissim$P_non_shared <- NA
  
  sites1_data <- TPD_plot_data(data, sites)
  ransites_data <- TPD_plot_data(randata, sites)
  
  
  TPD_i <- sites1_data[["pl_dat"]][["prob"]]
  
  
  TPD_j <- ransites_data[["pl_dat"]] %>%
    percentile_cells() %>% dplyr::mutate(prob = ifelse(percentile >= threshold, 0, prob)) %>% pull(prob)
  
  
  
  
  
  
  O_aux <- sum(pmin(TPD_i, TPD_j))
  shared_aux <- which(TPD_i > 0 & TPD_j > 0)
  A_aux <- sum(pmax(TPD_i[shared_aux], TPD_j[shared_aux])) -
    O_aux
  only_in_i_aux <- which(TPD_i > 0 & TPD_j ==
                           0)
  B_aux <- sum(TPD_i[only_in_i_aux])
  only_in_j_aux <- which(TPD_i == 0 & TPD_j >
                           0)
  C_aux <- sum(TPD_j[only_in_j_aux])
  results_samp$dissim$dissimilarity <-
    results_samp$dissim$dissimilarity <- 1 - O_aux
  results_samp$dissim$P_non_shared <-
    results_samp$dissim$P_non_shared     <-
    (2 * min(B_aux, C_aux)) / (A_aux +
                                 2 * min(B_aux, C_aux))
  results_samp$dissim$P_shared <-
    results_samp$dissim$P_shared  <-
    1 - results_samp$dissim$P_non_shared
  
  
  return(results_samp)
  
}

######################################
######################################
#####################################



species_fit <- function(cells) {
  species_cells_frame <-
    data.frame(cells[, c(1:3)], occupying_species = NA)
  
  
  
  for (i in 1:nrow(species_cells_frame)) {
    print(i)
    
    x <- cells[i, 1]
    y <- cells[i, 2]
    z <- cells[i, 3]
    
    
    occ_frame <-
      species_TPD %>% dplyr::filter(locomotion == x, foraging == y, body == z)
    
    
    
    
    species_cells_frame[i, 4] <-
      paste(gsub(
        names(which(colSums(occ_frame[, -c(1:3)]) > 0)),
        pattern = "\\.",
        replacement = " "
      ), collapse = "/")
    
    
  }
  return(species_cells_frame)
}

#########################################
#########################################
#########################################
#########################################


TPD_species_occupancy <- function(data, sites1, sites2, randata) {
  TPD_lost_gain <- list()
  cell_volume <- data[["data"]][["cell_volume"]]
  
  
  
  if (!all(sites1 %in% names(data)) |
      !all(sites1 %in% names(randata))) {
    sites1 <- sites1[which(sites1 %in% names(data))]
    sites1 <- sites1[which(sites1 %in% names(randata))]
  }
  
  if (!all(sites2 %in% names(data)) |
      !all(sites2 %in% names(randata))) {
    sites2 <- sites2[which(sites2 %in% names(data))]
    sites2 <- sites2[which(sites2 %in% names(randata))]
  }
  
  
  
  sites1_data <- TPD_plot_data(data, sites1)
  sites2_data <- TPD_plot_data(data, sites2)
  
  
  ####### species pools - should only consider areas of trait space that are occupied by each LU null commmunities
  
  null_sites1 <- TPD_plot_data(randata, sites1)
  null_sites2 <- TPD_plot_data(randata, sites2)
  
  
  both_occ <-
    which(null_sites1[["pl_dat"]][["prob"]] > 0 &
            null_sites2[["pl_dat"]][["prob"]] > 0)
  
  
  filled_cells_1 <-
    sites1_data[["pl_dat"]] %>% dplyr::slice(both_occ) %>% dplyr::group_by(T2, T1, T3) %>%
    dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% filter(!is.na(prob)) %>% data.frame() %>% percentile_cells()
  
  filled_cells_2 <-
    sites2_data[["pl_dat"]] %>% dplyr::slice(both_occ) %>% dplyr::group_by(T2, T1, T3) %>%
    dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% filter(!is.na(prob)) %>% data.frame() %>% percentile_cells() %>%
    dplyr::rename(prob_2 = prob, percentile_2 = percentile)
  
  ## join the two data frames together with the secondary sites prob and percentile renamed
  
  diff_cells <-
    filled_cells_1 %>% dplyr::left_join(filled_cells_2, by = c("T1", "T2", "T3")) %>% data.frame()
  
  ##### work out with difference in prob/percentile of cells and whether the occupancy of the cell is new (functional gain) or lost (functional loss)
  
  cells_frame <-
    filled_cells_2 %>% dplyr::left_join(filled_cells_1, by = c("T1", "T2", "T3")) %>% dplyr::filter(is.na(prob)) %>%
    dplyr::relocate(prob, .before = prob_2) %>% rbind(diff_cells) %>% dplyr::mutate(
      lost_cell = ifelse(is.na(prob_2), TRUE, FALSE),
      gain_cell = ifelse(is.na(prob), TRUE, FALSE),
      prob = ifelse(is.na(prob), 0, prob),
      prob_2 = ifelse(is.na(prob_2), 0, prob_2),
      diff = prob_2 - prob
    )
  
  
  gain_species <-
    cells_frame %>% dplyr::filter(gain_cell) %>% dplyr::select(T1, T2, T3)
  gain_sp <- species_fit(gain_species)
  
  lost_species <-
    cells_frame %>% dplyr::filter(lost_cell) %>% dplyr::select(T1, T2, T3)
  lost_sp <- species_fit(lost_species)
  
  
  
  TPD_lost_gain$gain$gain_species <- gain_sp
  TPD_lost_gain$gain$gain_volume <-
    nrow(gain_species) * cell_volume
  
  TPD_lost_gain$lost$lost_species <- lost_sp
  TPD_lost_gain$lost$lost_volume <- nrow(lost_species) * cell_volume
  
  
  return(TPD_lost_gain)
}


################################################################################
#################################################################################
## Function to detect holes ###################################################




TPD_holes <-
  function(data = NULL,
           randata = NULL,
           sites,
           threshold) {
    ### extract the TPD data for teh two sets of sites
    
    TPD_holes_list <- list()
    
    
    if (is.null(data)) {
      cellvolume <- randata[[sites]][["data"]][["cell_volume"]]
    } else {
      cellvolume <- data[["data"]][["cell_volume"]]
    }
    distancecharacteristic <- 1 / (cellvolume ^ (1 / 3))
    
    
    if (!is.null(data)) {
      sites <- sites[which(sites %in% names(data))]
    }
    
    sites <- sites[which(sites %in% names(randata))]
    
    
    if (is_empty(sites)) {
      return(TPD_holes_list)
    }
    
    
    if (!is.null(data)) {
      sites1_data <- TPD_plot_data(data, sites)
      sites2_data <- TPD_plot_data(randata, sites)
      
      
      
      
      ## for the 3D plot need to just get the cells which are functionally occupied
      
      filled_cells_1 <-
        sites1_data[["pl_dat"]] %>% dplyr::group_by(T2, T1, T3) %>%
        dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% filter(!is.na(prob)) %>% data.frame()
      
      filled_cells_1 <- percentile_cells(filled_cells_1)
      
      filled_cells_2 <-
        sites2_data[["pl_dat"]] %>% dplyr::group_by(T2, T1, T3) %>%
        dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% filter(!is.na(prob)) %>% data.frame()
      
      filled_cells_2 <-
        percentile_cells(filled_cells_2) %>% dplyr::rename(prob_2 = prob, percentile_2 = percentile) %>%
        dplyr::filter(percentile_2 <= threshold)
      
      
      
    } else {
      sites1_data <-
        data.frame(
          T1 = randata[[sites]][["data"]][["evaluation_grid"]][[1]],
          T2 = randata[[sites]][["data"]][["evaluation_grid"]][[2]],
          T3 = randata[[sites]][["data"]][["evaluation_grid"]][[3]],
          prob = randata[[sites]][["obsTPDc"]][["RelativeAbundance"]]
        )
      
      sites2_data <-
        data.frame(
          T1 = randata[[sites]][["data"]][["evaluation_grid"]][[1]],
          T2 = randata[[sites]][["data"]][["evaluation_grid"]][[2]],
          T3 = randata[[sites]][["data"]][["evaluation_grid"]][[3]],
          prob = randata[[sites]][["ranTPDc"]][["RelativeAbundance"]]
        )
      
      
      
      filled_cells_1 <-
        sites1_data %>% dplyr::group_by(T2, T1, T3) %>%
        dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% filter(!is.na(prob)) %>% data.frame()
      
      filled_cells_1 <- percentile_cells(filled_cells_1)
      
      filled_cells_2 <-
        sites2_data %>% dplyr::group_by(T2, T1, T3) %>%
        dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% filter(!is.na(prob)) %>% data.frame()
      
      filled_cells_2 <-
        percentile_cells(filled_cells_2) %>% dplyr::rename(prob_2 = prob, percentile_2 = percentile) %>%
        dplyr::filter(percentile_2 <= threshold)
      
    }
    
    ####################################
    # WORK OUT THE DIFFERENCE ##########
    ####################################
    
    ## join the two data frames together with the secondary sites prob and percentile renamed
    
    diff_cells <-
      filled_cells_1 %>% dplyr::left_join(filled_cells_2, by = c("T1", "T2", "T3")) %>% data.frame()
    
    
    ##### work out with difference in prob/percentile of cells and whether the occupancy of the cell is new (functional gain) or lost (functional loss)
    
    cells_frame <-
      filled_cells_2 %>% dplyr::left_join(filled_cells_1, by = c("T1", "T2", "T3")) %>% dplyr::filter(is.na(prob)) %>%
      dplyr::relocate(prob, .before = prob_2) %>% rbind(diff_cells) %>% dplyr::mutate(
        lost_cell = ifelse(is.na(prob_2), TRUE, FALSE),
        gain_cell = ifelse(is.na(prob), TRUE, FALSE),
        prob = ifelse(is.na(prob), 0, prob),
        prob_2 = ifelse(is.na(prob_2), 0, prob_2),
        diff = prob_2 - prob
      )
    ###########################
    # minimum convex hull surrounding observed site
    
    if (nrow(filled_cells_1) < 4) {
      return(TPD_holes_list)
    }
    
    cvh <-
      convhulln(filled_cells_1[, c("T1", "T2", "T3")], options = "FA")
    cvh_vol <- cvh$vol
    obs_cvh <- convhulln(filled_cells_1[, c("T1", "T2", "T3")])
    
    
    ########################
    #### Total cells within
    #### convex hull
    ########################
    
    total_int_cells <-
      cells_frame[inhulln(obs_cvh, as.matrix(cells_frame[, c("T1", "T2", "T3")])), ]
    
    ###################################
    
    absent_cells <-
      cells_frame %>% dplyr::filter(gain_cell) %>% dplyr::select(T1, T2, T3)
    absent_cells$internal <- inhulln(obs_cvh, as.matrix(absent_cells))
    absent_cells$external <- !absent_cells$internal
    
    abs_int <- absent_cells %>% dplyr::filter(internal)
    abs_ext <- absent_cells %>% dplyr::filter(external)
    
    if (nrow(abs_int) == 0) {
      return(TPD_holes_list)
    }
    
    if (nrow(abs_ext) == 0) {
      return(TPD_holes_list)
    }
    
    
    if (nrow(abs_int) < 2) {
      holes_int <- data.frame(abs_int[, c(1:3)], holes = 1)
    } else {
      clusters_int <- dist(abs_int[, c(1:3)], method = "euclidean") ^ 2
      tree_int <- fastcluster::hclust(clusters_int, method = "ward.D2")
      holes_int <-
        data.frame(abs_int[, c(1:3)], holes = cutree(tree_int, h = distancecharacteristic))
      
    }
    
    if (nrow(abs_ext) < 2) {
      holes_ext <- data.frame(abs_ext[, c(1:3)], holes = 1)
    } else {
      clusters_ext <- dist(abs_ext[, c(1:3)], method = "euclidean") ^ 2
      tree_ext <- fastcluster::hclust(clusters_ext, method = "ward.D2")
      holes_ext <-
        data.frame(abs_ext[, c(1:3)], holes = cutree(tree_ext, h = distancecharacteristic))
      
    }
    
    
    
    ########
    #######
    #### so some metrics to get from these 1) number of holes 2) mean size of holes 3) proportion of internal volume to holes 4)
    
    holes_int_metrics <-
      holes_int %>% dplyr::group_by(holes) %>% dplyr::summarise(size = n()) %>% ungroup() %>%
      dplyr::mutate(
        hole__richness = cellvolume * size,
        proportion = (cellvolume * size) / (nrow(total_int_cells) * cellvolume)
      )
    
    holes_ext_metrics <-
      holes_ext %>% dplyr::group_by(holes) %>% dplyr::summarise(size = n()) %>% ungroup() %>%
      dplyr::mutate(
        hole__richness = cellvolume * size,
        proportion = (cellvolume * size) / (nrow(cells_frame) * cellvolume)
      )
    
    TPD_holes_list$internal$internal_hole_cells <- abs_int[, c(1:3)]
    TPD_holes_list$internal$number_of_holes <- nrow(holes_int_metrics)
    TPD_holes_list$internal$total_hole_volume <-
      sum(holes_int_metrics$hole__richness)
    TPD_holes_list$internal$proportion_holes_volume <-
      sum(holes_int_metrics$proportion)
    
    TPD_holes_list$external$external_hole_cells <- abs_ext[, c(1:3)]
    TPD_holes_list$external$number_of_holes <- nrow(holes_ext_metrics)
    TPD_holes_list$external$total_hole_volume <-
      sum(holes_ext_metrics$hole__richness)
    TPD_holes_list$external$proportion_holes_volume <-
      sum(holes_ext_metrics$proportion)
    
    
    return(TPD_holes_list)
    
  }


TPD_site_check <- function(data, LU, realm) {
  TPD_sample_check <- list()
  
  lu_sites <-
    TPD_LU %>% dplyr::filter(Realm == realm, Predominant_habitat == LU) %>% dplyr::distinct(SSBS) %>% pull()
  
  if (length(lu_sites) == 0) {
    return(TPD_sample_check)
  } else {
    cumulative_df <-
      matrix(rep(NA, (4 * length(lu_sites))), ncol = 4) %>% data.frame()
    colnames(cumulative_df) <- c("Sample", "FRich", "FEve", "FDiv")
    cumulative_df$Sample <- 1:length(lu_sites)
    
    cumul_site <- c()
    i <- 1
    for (sit in sample(lu_sites, replace = FALSE)) {
      cumul_site <- c(cumul_site, sit)
      
      fd_met <- TPD_FD_metrics(data, sites = cumul_site)
      fd_met <- data.frame(fd_met["total_sites", ])
      cumulative_df[i, c(2:4)] <- as.numeric(fd_met[1, 1:3])
      i <- i + 1
    }
    
    
    
    line_plot_FRich <-
      ggplot(data = cumulative_df, aes(x = Sample, y = FRich)) +
      geom_line() +
      geom_point()
    
    line_plot_FDiv <-
      ggplot(data = cumulative_df, aes(x = Sample, y = FDiv)) +
      geom_line() +
      geom_point()
    
    line_plot_FEve <-
      ggplot(data = cumulative_df, aes(x = Sample, y = FEve)) +
      geom_line() +
      geom_point()
    
    multi_plot <-
      ggarrange(line_plot_FRich, line_plot_FEve, line_plot_FDiv)
    
    
    
    TPD_sample_check$cumulative_df <- cumulative_df
    TPD_sample_check$ggplots <- multi_plot
    
    return(TPD_sample_check)
    
  }
}

###########################################################
###########################################################
###########################################################

########################################################
########################################################



pool_sp <- function(df, sites, pool) {
  pool_species <- c()
  for (sit in sites) {
    pool_species <-
      unique(c(pool_species, species_pools[[sit]][["0.9"]]))
  }
  
  for (i in 1:nrow(df)) {
    sp <- unlist(strsplit(df[i, "occupying_species"], split = "/"))
    df[i, "occupying_species"] <-
      paste(sp[which(sp %in% pool_species)], collapse = "/")
  }
  
  return(df)
  
}


###################################################
###################################################

split_sp_func <- function(df) {
  species <- c()
  for (i in 1:nrow(df)) {
    species <-
      unique(c(species, unlist(strsplit(df[i, "occupying_species"], split = "/"))))
  }
  return(species)
}

##############################################################
##############################################################


TPD_forage_mapping_data <-
  function(data,
           randata = NULL,
           for_data,
           sites,
           guilds = c("Fr",
                      "Gr",
                      "Ne",
                      "In",
                      "Vt",
                      "Aq.p",
                      "Sc",
                      "Hb.A",
                      "Hb.T",
                      "Om",
                      "Unclassified")) {
    TPD_for_map_data <- list()
    ## Create a colour list for each of the dietary guilds
    
    trophic_niche_colours <- list()
    trophic_niche_colours$Fr <- "red2"
    trophic_niche_colours$Gr <- "orange2"
    trophic_niche_colours$Ne <- "olivedrab"
    trophic_niche_colours$In <- "mediumblue"
    trophic_niche_colours$Vt <- "purple1"
    trophic_niche_colours$Aq.p <- "steelblue"
    trophic_niche_colours$Sc <- "goldenrod4"
    trophic_niche_colours$Hb.A <- "turquoise4"
    trophic_niche_colours$Hb.T <- "springgreen4"
    trophic_niche_colours$Om <- "white"
    trophic_niche_colours$Unclassified <- "grey25"
    
    ### some sites didn't make the cut to the data due to deficiencies so drop them
    
    if (!is.null(randata) &
        (!all(sites %in% names(data)) | !all(sites %in% names(randata)))) {
      sites <- sites[which(sites %in% names(data))]
      sites <- sites[which(sites %in% names(randata))]
    }
    
    if (!all(sites %in% names(data))) {
      sites <- sites[which(sites %in% names(data))]
      
    }
    
    
    ## collate the TPD plot data from the TPD
    
    data_3d <- TPD_plot_data(data, sites)
    
    
    # first identify which cells in 3D space are occupied
    
    filled_cells <-
      data_3d[["pl_dat"]] %>% dplyr::group_by(T2, T1, T3) %>%
      dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% filter(!is.na(prob)) %>% data.frame()
    
    
    if (!is.null(randata)) {
      ran_data_3d <- TPD_plot_data(randata, sites)
      
      ran_filled_cells <-
        ran_data_3d[["pl_dat"]] %>% dplyr::group_by(T2, T1, T3) %>%
        dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% filter(!is.na(prob)) %>% data.frame()
      
      species_tpd <- species_fit(ran_filled_cells)
    }
    
    
    ## work out what species could occupy areas of trait space
    
    
    
    #### of all the species that could occupy those areas of trait space which are in those sites species pools.
    
    species_tpd <-
      pool_sp(df = species_tpd, sites = sites, pool = "0.9")
    
    ## add columns to the df to record how many species of each dietary guil are present in these areas .
    
    species_tpd$Fr <- NA
    species_tpd$Gr <- NA
    species_tpd$Ne <- NA
    species_tpd$In <- NA
    species_tpd$Vt <- NA
    species_tpd$Aq.p <- NA
    species_tpd$Sc <- NA
    species_tpd$Hb.A <- NA
    species_tpd$Hb.T <- NA
    species_tpd$Om <- NA
    species_tpd$Trophic_niche <- NA
    species_tpd$colour <- NA
    
    
    ## for each row of the df I am going to determine the frequency of species of different dietary guilds and whether that particular cell
    ## can be classified to be a certain deitary guild and in line with classifying species to a dietary guild if >70% of species are of
    ## a particular guil then the cell in trait space is classified as such.
    
    for (i in 1:nrow(species_tpd)) {
      # data frame of the species and join on foraging data
      
      spp <-
        data.frame(Birdlife_Name = split_sp_func(species_tpd[i, ])) %>% dplyr::left_join(for_data[, c("Birdlife_Name", "Trophic_Niche", "Foraging_Niche")], by = "Birdlife_Name") %>%
        dplyr::group_by(Trophic_Niche) %>% dplyr::summarise(count = n())
      
      # if there is only a single row then you can classify that cell as that guild and record the count
      p_guilds <- spp$Trophic_Niche
      
      for (g in p_guilds) {
        if (nrow(spp) == 1) {
          species_tpd[i, "Trophic_niche"] <- g
          
        }
        ## guild count
        
        species_tpd[i, g] <-
          spp %>% filter(Trophic_Niche == g) %>% pull(count)
      }
      
      
      if (is.na(species_tpd[i, "Trophic_niche"])) {
        species_tpd[i, "Trophic_niche"] <- "Unclassified"
        
        
        for (g in p_guilds) {
          if (species_tpd[i, g] / sum(species_tpd[i, p_guilds]) >= 0.7) {
            species_tpd[i, "Trophic_niche"] <- g
          }
        }
        
      }
      
      
      ##### colours will be assigned if over 70% of species is of that trophic niche
      
      
      species_tpd[i, "colour"] <-
        trophic_niche_colours[[species_tpd[i, "Trophic_niche"]]]
      
    }
    
    if (!is.null(randata)) {
      TPD_for_map_data$random <- species_tpd
      
      obs_species_tpd <-
        filled_cells %>% merge(species_tpd, by = c("T1", "T2", "T3"))
      
      TPD_for_map_data$observed <- obs_species_tpd
      
    } else {
      TPD_for_map_data$observed <- species_tpd
    }
    
    return(TPD_for_map_data)
    
  }




TPD_forage_mapping_plot <-
  function(data,
           tpd_for_map,
           T1lab,
           T2lab,
           T3lab,
           guilds = "all",
           title,
           save = FALSE,
           s_file,
           animation = FALSE,
           a_file) {
    data_3d <- data[["data"]][["evaluation_grid"]]
    
    
    if (guilds != "all") {
      data <- tpd_for_map %>% dplyr::filter(Trophic_niche %in% guilds)
    } else {
      data <- tpd_for_map
    }
    
    
    x <- data$T2
    y <- data$T1
    z <- data$T3
    
    
    
    #####################################
    
    
    #######################################
    
    
    xmax <- max(data_3d[, 2])
    xmin <- min(data_3d[, 2])
    ymax <- max(data_3d[, 1])
    ymin <- min(data_3d[, 1])
    zmax <- max(data_3d[, 3])
    zmin <- min(data_3d[, 3])
    
    scale <- mean(c(dist(c(xmin, xmax))[1],
                    dist(c(ymin, ymax))[1],
                    dist(c(zmin, zmax))[1])) / 100
    
    
    clear3d()
    rgl.viewpoint(theta = 50,
                  phi = 25,
                  zoom = 13 / 16)
    plot3d(
      x,
      y,
      z,
      box = FALSE,
      xlab = "",
      ylab = "",
      zlab = "",
      type = "s",
      radius = scale,
      alpha = 0.8,
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      zlim = c(zmin, zmax),
      col = data$colour,
      lwd = 0.1
    )
    title3d(
      main = title,
      xlab = T2lab,
      ylab = T1lab,
      zlab = T3lab
    )
    
    
    if (save) {
      rgl.snapshot(s_file, fmt = 'png')
    }
    
    if (animation) {
      movie3d(
        movie = a_file,
        spin3d(
          axis = c(0, 1, 0),
          rpm = 10,
          dev = cur3d()
        ),
        startTime = 0,
        duration = 12,
        dir = ".",
        type = "gif",
        clean = T,
        fps = 2
      )
    }
    
    
  }

#####################################
############################################
#### Next function is partitioning trait space to areas realted to a specific dietary guil and calculating beta diversity to see how
#### difference areas of trait space are impacted by land use change compared to null expectations


proportional_occupancy <- function(data,
                                   randata,
                                   fordata,
                                   ranfordata,
                                   sites,
                                   guilds = c("Fr",
                                              "Gr",
                                              "Ne",
                                              "In",
                                              "Vt",
                                              "Aq.p",
                                              "Sc",
                                              "Hb.A",
                                              "Hb.T",
                                              "Om",
                                              "Unclassified"),
                                   standardise = FALSE) {
  if (!all(sites %in% names(data)) | !all(sites %in% names(randata))) {
    sites <- sites[which(sites %in% names(data))]
    sites <- sites[which(sites %in% names(randata))]
  }
  
  
  prop_occ_df <- c()
  
  sites1_data <- TPD_plot_data(data, sites)
  ransites_data <- TPD_plot_data(randata, sites)
  
  
  
  
  for (g in guilds) {
    ran_guild_data <- ranfordata %>% dplyr::filter(Trophic_niche == g)
    
    guild_cells <- c()
    for (i in 1:nrow(ran_guild_data)) {
      cell <-
        as.numeric(
          which(
            sites1_data[["pl_dat"]]$T1 == ran_guild_data[i, "T1"] &
              sites1_data[["pl_dat"]]$T2 == ran_guild_data[i, "T2"] &
              sites1_data[["pl_dat"]]$T3 == ran_guild_data[i, "T3"]
          )
        )
      
      guild_cells <- c(guild_cells, cell)
    }
    
    
    TPD_i <- sites1_data[["pl_dat"]][guild_cells, "prob"]
    
    TPD_j <- ransites_data[["pl_dat"]] %>%
      percentile_cells()  %>% slice(guild_cells) %>% pull(prob)
    
    if (standardise) {
      TPD_i <- TPD_i / (sum(TPD_i))
      TPD_j <- TPD_j / (sum(TPD_j))
    }
    
    
    O_aux <- sum(pmin(TPD_i, TPD_j))
    
    results_samp <- data.frame(
      similarity = O_aux,
      proportional_occupancy = sum(TPD_i > 0) / sum(TPD_j > 0),
      Trophic_niche = g
    )
    
    
    prop_occ_df <- rbind(prop_occ_df, results_samp)
  }
  return(prop_occ_df)
}

#########################################
########################################
###### re ordering combinations


reorder_combinations <- function(combo_df) {
  priority_order <- c(
    "Primary vegetation/Secondary vegetation",
    "Primary vegetation/Plantation forest",
    "Primary vegetation/Pasture",
    "Primary vegetation/Cropland",
    "Primary vegetation/Minimal agriculture",
    "Primary vegetation/Intensive agriculture",
    "Primary vegetation/Urban",
    "Secondary vegetation/Plantation forest",
    "Secondary vegetation/Pasture",
    "Secondary vegetation/Cropland",
    "Secondary vegetation/Minimal agriculture",
    "Secondary vegetation/Intensive agriculture",
    "Secondary vegetation/Urban",
    "Plantation forest/Pasture",
    "Plantation forest/Cropland",
    "Plantation forest/Minimal agriculture",
    "Plantation forest/Intensive agriculture",
    "Plantation forest/Urban",
    "Pasture/Cropland",
    "Pasture/Urban",
    "Minimal agriculture/Intensive agriculture",
    "Minimal agriculture/Urban",
    "Intensive agriculture/Urban",
    "Cropland/Urban"
  )
  
  for (i in 1:nrow(combo_df)) {
    LU1 <- combo_df[i, 1]
    LU2 <- combo_df[i, 2]
    
    if (LU1 == LU2) {
      next()
    }
    
    priority <-
      grep(priority_order, pattern = LU1, value = TRUE)[grep(priority_order, pattern = LU1) %in% grep(priority_order, pattern = LU2)]
    
    combo_df[i, 1] <-  unlist(str_split(priority, pattern = "/"))[1]
    combo_df[i, 2] <-  unlist(str_split(priority, pattern = "/"))[2]
    
    
  }
  return(combo_df)
  
}



####################################
####################################
# ####################################


centre_of_mass_for <-
  function(data, sites = NULL, tpd_for_map, threshold) {
    sites_data <- TPD_plot_data(data, sites)
    
    
    prop_COM <- c()
    COM_coords <- c()
    
    
    functional_volume <-
      sites_data[["pl_dat"]] %>% dplyr::group_by(T2, T1, T3) %>%
      dplyr::mutate(prob = ifelse(prob == 0, NA, prob))  %>% data.frame() %>% percentile_cells() %>% dplyr::filter(!is.na(percentile))
    
    full_guild_prop <- functional_volume %>%
      merge(tpd_for_map[, colnames(tpd_for_map) != "prob"], by = c("T1", "T2", "T3")) %>% dplyr::group_by(Trophic_niche) %>%
      dplyr::mutate(
        full_probability = sum(prob),
        n_cells = n(),
        total_cells = nrow(functional_volume),
        full_prop_COM = n_cells / total_cells
      ) %>%
      dplyr::distinct(Trophic_niche, full_prop_COM, full_probability)
    
    
    for (i in 1:length(threshold)) {
      percentile_vol <-
        functional_volume %>% dplyr::filter(percentile <= threshold[i]) %>%
        dplyr::select(T1, T2, T3, prob) %>% dplyr::mutate(percentile = threshold[i])
      
      if (nrow(percentile_vol) == 0) {
        next()
      }
      
      
      percentile_vol$prob <-
        percentile_vol$prob / sum(percentile_vol$prob)
      
      weights <- percentile_vol %>% pull(prob)
      
      
      
      
      
      prop_species <-
        tpd_for_map[, colnames(tpd_for_map) != "prob"] %>% merge(percentile_vol, by = c("T1", "T2", "T3")) %>%
        dplyr::group_by(Trophic_niche) %>%
        dplyr::summarise(n_cells = n()) %>% ungroup() %>% dplyr::mutate(prop_COM = n_cells /
                                                                          sum(n_cells),
                                                                        percentile = threshold[i]) %>%
        merge(full_guild_prop) %>% dplyr::mutate(relative_prop_COM = prop_COM /
                                                   full_prop_COM) %>%
        dplyr::select(-n_cells) %>%
        data.frame()
      
      prob_species <-
        tpd_for_map[, colnames(tpd_for_map) != "prob"] %>% merge(percentile_vol, by = c("T1", "T2", "T3")) %>%
        group_by(Trophic_niche) %>%
        dplyr::summarise(prob_COM = sum(prob), percentile = threshold[i]) %>%
        data.frame()
      
      
      percentile_species <-
        cbind(prop_species, prob_COM = prob_species[, c("prob_COM")]) %>%
        dplyr::mutate(relative_prob_COM = prob_COM / full_probability)
      
      
      guilds <-
        c("Om",
          "Gr",
          "Fr",
          "In",
          "Ne",
          "Hb.T",
          "Hb.A",
          "Aq.p",
          "Vt",
          "Sc")
      
      if (any(!(guilds %in% percentile_species$Trophic_niche))) {
        for (g in guilds[!(guilds %in% percentile_species$Trophic_niche)]) {
          g_row <- percentile_species[1, ]
          g_row$Trophic_niche <- g
          g_row[, !(colnames(g_row) %in% c("Trophic_niche", "percentile"))] <-
            0
          
          percentile_species <- rbind(percentile_species, g_row)
        }
        
      }
      
      
      
      COM_coords <- rbind(COM_coords, percentile_vol)
      
      
      COM_point <-
        c(apply(percentile_vol[, c("T1", "T2", "T3")], 2, weighted.mean, w = weights)) /
        sum(weights)
      
      
      
      COM_sd <- c()
      for (col in 1:ncol(percentile_vol)) {
        COM_sd <-
          c(COM_sd,
            Hmisc::wtd.var(
              x = percentile_vol[, col],
              weights = weights,
              normwt = TRUE
            ))
      }
      
      
      percentile_species$T1 <- COM_point[1]
      percentile_species$sdT1 <- COM_sd[1]
      
      percentile_species$T2 <- COM_point[2]
      percentile_species$sdT2 <- COM_sd[2]
      
      percentile_species$T3 <- COM_point[3]
      percentile_species$sdT3 <- COM_sd[3]
      
      
      prop_COM <- rbind(prop_COM, percentile_species)
      
      
    }
    
    COM_list <- list(prop_COM, COM_coords)
    names(COM_list) <- c("prop_COM", "COM_coords")
    
    return(COM_list)
  }



#################
#################
#################
# Battacharyya distance between hypervolumes as a measure of dissimilarity



B_dist <- function(data, sites) {
}
