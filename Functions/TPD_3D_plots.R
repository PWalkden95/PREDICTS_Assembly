
markdown_rds_open <- function(path) {
  open_file <- try(readRDS(path), silent = TRUE)
  
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
      pl_dat[, 4] + data[[sit]][["TPDc"]][[1]]
  }
  
  
  pl_dat[, 4] <- pl_dat[, 4] / length(site)
  colnames(pl_dat) <- c("y", "x", "z", "prob")
  pl_dat <- data.frame(pl_dat)
  
  xy_dat <-
    pl_dat %>% dplyr::group_by(x, y) %>% dplyr::summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% set_colnames(c("x", "y", "prob")) %>% data.frame() %>% percentile_cells() %>%
    dplyr::filter(!is.na(x))
  
  yz_dat <-
    pl_dat %>% dplyr::group_by(y, z) %>% dplyr::summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% set_colnames(c("y", "z", "prob")) %>% data.frame() %>% percentile_cells() %>%
    dplyr::filter(!is.na(y))
  
  xz_dat <-
    pl_dat %>% dplyr::group_by(x, z) %>% dplyr::summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% set_colnames(c("x", "z", "prob")) %>% data.frame() %>% percentile_cells() %>%
    dplyr::filter(!is.na(x))
  
  
  trait_plot_data <- list()
  trait_plot_data$pl_dat <- pl_dat
  trait_plot_data$xy_dat <- xy_dat
  trait_plot_data$yz_dat <- yz_dat
  trait_plot_data$xz_dat <- xz_dat
  
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
  occupied_cells <- data %>% dplyr::group_by(y,x,z) %>%
    dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% dplyr::filter(!is.na(prob)) %>% data.frame
  
  return(occupied_cells)
}


#####################

grid_function <- function(plot_data, axes = NULL){
  
  other_axes <- colnames(plot_data[["pl_dat"]][1:3])[which(!(colnames(plot_data[["pl_dat"]])[1:3] %in% axes))]
  
  
  trait_dataframe <- plot_data[[which(!grepl(names(plot_data)[-1], pattern = other_axes)) + 1]]
  
  ## create a white colour matrix that is our default and blank
  
  axis_colour_matrix <- matrix(rep("#FFFFFF", 50 * 50),
                               nrow = 50, ncol = 50)
  
  #### create a matrix that takes the mean occupancy of cells 
  
  axis_value_matrix <- matrix(rep(unique(plot_data[["pl_dat"]][,other_axes]), 50 * 50),
                              nrow = 50, ncol = 50, dimnames = list(unique(trait_dataframe[,axes[1]]),
                                                                    unique(trait_dataframe[,axes[2]]))) 
  
  #### for each row in the dataframe fill in the values of the mean value for the two dimensions
  
  for(i in 1:nrow(trait_dataframe)){
    axis_value_matrix[as.character(trait_dataframe[i,axes[1]]),as.character(trait_dataframe[i,axes[2]])] <- 
      ifelse(method == "prob",
             trait_dataframe[i,"prob"],
             trait_dataframe[i,"percentile"])
  }
  
  ## now to get the colours for the grid cells
  
  ## extract the values from the grid
  
  grid_values <- axis_value_matrix[!is.na(axis_value_matrix)]
  
  
  ## rows are x
  ##  columns are y
  
  axes_colours <- data.frame(value = grid_values[order(grid_values)], colour_name = my_colour_2(length(grid_values))) 
  
  axis_colour_matrix[which(!is.na(axis_value_matrix))[order(grid_values)]] <- axes_colours[,"colour_name"]
  
  axis_colour_matrix <- as.character(axis_colour_matrix)
  
  
  return(colour_matrix = axis_colour_matrix)
  
  
}

######################################

limits_function <- function(data){
  
    limits <- c(
      xmax = max(data$x),
      xmin = min(data$x),
      ymax = max(data$y),
      ymin = min(data$y),
      zmax = max(data$z),
      zmin = min(data$z))

  
  return(limits)
}

###########################
scale_function <- function(limits, factor){
  scale <- mean(c(dist(c(
    limits["xmin"], limits["xmax"]
  ))[1],
  dist(c(
    limits["ymin"], limits["xmax"]
  ))[1],
  dist(c(
    limits["zmin"], limits["zmax"]
  ))[1])) / factor
  
  return(scale)
}


#########################


TPD_3d_plot <-
  function(data,
           sites = NULL,
           xlab = "",
           ylab = "",
           zlab = "",
           method = "prob",
           title = "",
           save = FALSE,
           s_file,
           limits = limits,
           scale = scale,
           grid = FALSE,
           free_limits = FALSE, 
           animation = FALSE,
           a_file) {
    
    
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
      
      
      
      x <- filled_cells[,"x"]
      y <- filled_cells[,"y"]
      z <- filled_cells[,"z"]
      
      
      
      
      
      #creates a function my.colors which interpolates n colors between blue, white and red
      color.df <-
        data.frame(value = TPD_colours[order(TPD_colours)], color.name =
                     my.colors(length(TPD_colours))) %>% distinct(value, .keep_all = TRUE)#generates 2001 colors from the color ramp
     
      
      filled_cells_col <- filled_cells %>% 
        dplyr::mutate(color.name = color.df$color.name[unlist(lapply(filled_cells$value,find_position, y = color.df$value))])
      

  
      
     
      #######################################
      ###########################################
      if (grid) {
        
        
        xy_colour_matrix <- grid_function(plot_data = data_3d,axes = c("x","y"))
        
        yz_colour_matrix <- grid_function(plot_data = data_3d,axes = c("z","y"))
        
        xz_colour_matrix <- grid_function(plot_data = data_3d, axes = c("x","z")) 
        

      }
    #######################################
    ###################################
    
    
    if (free_limits) {

      limits <- limits_function(data = filled_cells)
      
      
      scale <- scale_function(limits = limits, factor = 75)
      
    } else {
      
      limits <- limits_function(data = data_3d[["pl_dat"]])
      
      scale <- scale_function(limits = limits, factor = 100)
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
        x = unique(data_3d[["xy_dat"]]$x),
        y = unique(data_3d[["xy_dat"]]$y),
        z = matrix(rep(min(data_3d[["pl_dat"]]$z), 50 * 50),
                   nrow = 50, ncol = 50),
        alpha = 0.5,
        lit = FALSE,
        front = "lines",
        back = "lines"
      )
      surface3d(
        x = unique(data_3d[["xy_dat"]]$x),
        y = unique(data_3d[["xy_dat"]]$y),
        z = matrix(rep(min(data_3d[["pl_dat"]]$z), 50 * 50),
                   nrow = 50, ncol = 50),
        color = xy_colour_matrix,
        smooth = FALSE,
        lit = FALSE
      )
      surface3d(
        x = rep(min(data_3d[["xy_dat"]]$x), 50),
        y = unique(data_3d[["xy_dat"]]$y),
        z = matrix(rep(unique(data_3d[["pl_dat"]]$z), 50),
                   nrow = 50, ncol = 50),
        alpha = 0.5,
        lit = FALSE,
        front = "lines",
        back = "lines"
      )
      surface3d(
        x = rep(min(data_3d[["xy_dat"]]$x), 50),
        y = unique(data_3d[["xy_dat"]]$y),
        z = matrix(rep(unique(data_3d[["pl_dat"]]$z), 50),
                   nrow = 50, ncol = 50),
        color = yz_colour_matrix,
        lit = FALSE,
        smooth = FALSE
      )
      surface3d(
        x = unique(data_3d[["xy_dat"]]$x),
        y = rep(min(data_3d[["xy_dat"]]$y), 50),
        z = t(matrix(
          rep(unique(data_3d[["pl_dat"]]$z), 50),
          nrow = 50, ncol = 50
        )),
        alpha = 0.5,
        lit = FALSE,
        front = "lines",
        back = "lines"
      )
      surface3d(
        x = unique(data_3d[["xy_dat"]]$x),
        y = rep(min(data_3d[["xy_dat"]]$y), 50),
        z = t(matrix(
          rep(unique(data_3d[["pl_dat"]]$z), 50),
          nrow = 50, ncol = 50
        )),
        color = xz_colour_matrix,
        lit = FALSE,
        smooth = FALSE
      )
      
    }
    title3d(
      main = title,
      xlab = xlab,
      ylab = ylab,
      zlab = zlab
    )
    
    if (save) {
      rgl.snapshot(s_file, fmt = 'png')
    }
    
    if (animation) {
      movie3d(
        movie = a_file,
        spin3d(
          axis = c(0, 1, 0),
          rpm = 6,
          dev = cur3d()
        ),
        startTime = 0,
        duration = 10,
        dir = ".",
        type = "gif",
        clean = T,
        fps = 10
      )
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



##############################################################
##############################################################
##############################################################

################################################

site_data_functional_metrics <- function(data, sites){
  
  if (!is.null(data[["data"]])) {
    cell_volume <- data[["data"]][["cell_volume"]]
  }
  
  
  site_data <- data[c(sites)]
  
  if (!is.null(data[["data"]])) {
    sites_data <- TPD_plot_data(data, sites)
    site_data$total_sites$TPDc$RelativeAbundance <-
      sites_data[["pl_dat"]][["prob"]]
  }
  
  
  return(site_data)
  
}


############################################

volume_evaluation_extract <- function(data, extract = c("cell_volume","evaluation_grid")){
  
  
  if("cell_volume" %in% extract){
  value <- data[["data"]][["cell_volume"]]
    }
  
  if("evaluation_grid" %in% extract){
    value <- data[["data"]][["evaluation_grid"]]
    
  }
  
  return(value)
  
}


##########################################



Calc_FRich <-
  function(data,
           sites = NULL
           ) {
    
  if(!is.null(data[["data"]])){  
  cell_volume <- volume_evaluation_extract(data = data, extract = "cell_volume")    
  }  
  
  
 site_data <- site_data_functional_metrics(data = data, sites = sites)
    
    
    
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
           sites = NULL
           ) {
    
    
      site_data <- site_data_functional_metrics(data = data, sites = sites)
    
    
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
           sites = NULL) {
    
    
    if (!is.null(data[["data"]])) {
      evaluation_grid <- volume_evaluation_extract(data = data, extract = "evaluation_grid")
      cell_volume <- volume_evaluation_extract(data = data, extract = "cell_volume")
    }
    
    

    site_data <- site_data_functional_metrics(data = data, sites = sites)
    
    
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
           sites = NULL) {
    
    if (!is.null(data[["data"]])) {
      evaluation_grid <- volume_evaluation_extract(data = data, extract = "evaluation_grid")
      cell_volume <- volume_evaluation_extract(data = data, extract = "cell_volume")
    }
    

      site_data <- site_data_functional_metrics(data = data, sites = sites)
      

    
    
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
           sites = NULL) {
    
    if (!is.null(data[["data"]])) {
      evaluation_grid <- volume_evaluation_extract(data = data, extract = "evaluation_grid")
    }
    

      site_data <- site_data_functional_metrics(data = data,sites = sites)
      
      
    
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
           sites = NULL) {
    
      FRich <- Calc_FRich(data, sites)
      FEve <- Calc_FEve(data, sites)
      FDiv <- Calc_FDiv(data, sites)
      FRound <- Calc_roundness(data, sites)
      FDis <- Calc_FDis(data, sites)
    
    
    
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
  
  
  TPD_i <- sites1_data[["pl_dat"]]  %>%
    percentile_cells() %>% dplyr::mutate(prob = ifelse(percentile >= threshold, 0, prob)) %>% pull(prob)
  
  
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

require(dbscan)
require(fpc)
require(factoextra)




TPD_holes <-
  function(data = NULL,
           randata = NULL,
           sites,
           threshold, 
           minimum_points,
          iteration) {
    ### extract the TPD data for teh two sets of sites
    
    TPD_holes_list <- list()
    
    
    
    
    
    
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
      
      
      if (is.null(data)) {
        cellvolume <- randata[[sites]][["data"]][["cell_volume"]]
      } else {
        cellvolume <- data[["data"]][["cell_volume"]]
      }
      
      
      cell_distance <- data.frame(x = c(abs(unique(sites1_data[["pl_dat"]][["x"]][order(sites1_data[["pl_dat"]][["x"]])])[1] -
            unique(sites1_data[["pl_dat"]][["x"]][order(sites1_data[["pl_dat"]][["x"]])])[2]),0), 
      
      y = c(abs(unique(sites1_data[["pl_dat"]][["z"]][order(sites1_data[["pl_dat"]][["z"]])])[1] -
            unique(sites1_data[["pl_dat"]][["z"]][order(sites1_data[["pl_dat"]][["z"]])])[2]), 0),
      
      z = c(abs(unique(sites1_data[["pl_dat"]][["y"]][order(sites1_data[["pl_dat"]][["y"]])])[1] -
            unique(sites1_data[["pl_dat"]][["y"]][order(sites1_data[["pl_dat"]][["y"]])])[2]), 0))
      
      
      cell_distance <- dist(cell_distance)[1]/2
      
      ## for the 3D plot need to just get the cells which are functionally occupied
      
      filled_cells_1 <-
        percentile_cells(hypervolume_occupied_cells(data = sites1_data[["pl_dat"]])) %>% dplyr::filter(percentile <= threshold)
      
      
      
      filled_cells_2 <-
        hypervolume_occupied_cells(data = sites2_data[["pl_dat"]])
      
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
        percentile_cells(hypervolume_occupied_cells(data = sites1_data[["pl_dat"]]))
      
      
      
      filled_cells_2 <-
        hypervolume_occupied_cells(data = sites2_data[["pl_dat"]])
      
      filled_cells_2 <-
        percentile_cells(filled_cells_2) %>% dplyr::rename(prob_2 = prob, percentile_2 = percentile)
      
      
    }
    
    ####################################
    # WORK OUT THE DIFFERENCE ##########
    ####################################
    
    ## join the two data frames together with the secondary sites prob and percentile renamed
    
    diff_cells <-
      filled_cells_1 %>% dplyr::left_join(filled_cells_2, by = c("x", "y", "z")) %>% data.frame()
    
    
    ##### work out with difference in prob/percentile of cells and whether the occupancy of the cell is new (functional gain) or lost (functional loss)
    
    cells_frame <-
      filled_cells_2 %>% dplyr::left_join(filled_cells_1, by = c("x", "y", "z")) %>% dplyr::filter(is.na(prob)) %>%
      dplyr::relocate(prob, .before = prob_2) %>% rbind(diff_cells) %>% dplyr::mutate(
        lost_cell = ifelse(is.na(prob_2), TRUE, FALSE),
        gain_cell = ifelse(is.na(prob), TRUE, FALSE),
        prob = ifelse(is.na(prob), 0, prob),
        prob_2 = ifelse(is.na(prob_2), 0, prob_2),
        diff = prob_2 - prob
      )
    

    
    positive_percentile <- cells_frame %>% dplyr::filter(diff >  0 ) %>% dplyr::arrange(desc(diff)) %>% 
      dplyr::mutate(diff = diff/sum(diff))
    positive_percentile[,"p_percentile"] <- cumsum(positive_percentile[,"diff"]) 
    positive_percentile <- positive_percentile %>% dplyr::filter(p_percentile <= threshold) %>% dplyr::select(-p_percentile)
    

    negative_percentile <- cells_frame %>% dplyr::filter(diff <  0 ) %>% dplyr::arrange(diff) %>% 
      dplyr::mutate(diff = diff/sum(diff))
    negative_percentile[,"p_percentile"] <- cumsum(negative_percentile[,"diff"]) 
    negative_percentile <- negative_percentile %>% dplyr::filter(p_percentile <= threshold) %>% dplyr::select(-p_percentile)
    
    
    
    cells_frame <- rbind(positive_percentile,negative_percentile)
    
    ###########################
    # minimum convex hull surrounding observed site
    
    if (nrow(filled_cells_1) <= 4) {
      return(TPD_holes_list)
    }
    
    cvh <-
      convhulln(filled_cells_1[, c("x", "y", "z")], options = "FA")
    cvh_vol <- cvh$vol
    obs_cvh <- convhulln(filled_cells_1[, c("x", "y", "z")])
    
    
    ########################
    #### Total cells within
    #### convex hull
    ########################
    
    total_int_cells <-
      cells_frame[inhulln(obs_cvh, as.matrix(cells_frame[, c("x", "y", "z")])), ]
    
    ###################################
    
    absent_cells <-
      cells_frame %>% dplyr::filter(gain_cell) %>% dplyr::select(x, y, z, prob_2)
    absent_cells$internal <- inhulln(obs_cvh, as.matrix(absent_cells[,c(1:3)]))
 
    
    abs_int <- absent_cells %>% dplyr::filter(internal)
   
    
    if (nrow(abs_int) == 0) {
      return(TPD_holes_list)
    }
    
    
    
    
    if (nrow(abs_int) < 2) {
      hole_size <-   data.frame(radius = NA,
                                         number_of_holes = NA,
                                         total_hole_volume = NA,
                                         mean_size = NA, 
                                         max_hole_size = NA, 
                                         min_hole_size = NA,
                                         total_proportion = NA,
                                         mean_proportion = NA,
                                         max_proportion = NA,
                                         min_proportion = NA,
                                         total_absence_proportion = NA,
                                         hypervolume_occupancy = NA)
    } else {

        
      hole_size <- c()
      
      hole_check <- 0
        for(i in seq(0,10,iteration)){
      
      hole_data <- data.frame(cluster = dbscan::dbscan(abs_int[,c(1:3)],eps = i, minPts = minimum_points)[["cluster"]]) %>% 
        dplyr::group_by(cluster) %>% dplyr::summarise(size = n()) %>% dplyr::mutate(cluster_volume = size * cellvolume, 
                                                                                    proportion = cluster_volume/ (nrow(total_int_cells)*cellvolume)) %>%
        dplyr::filter(cluster != 0)
      
      if(nrow(hole_data) == 0){
        next()
      }
      
      h_data <- data.frame(radius = i,
        number_of_holes = nrow(hole_data),
        total_hole_volume = sum(hole_data$size) * cellvolume,
                           mean_size = mean(hole_data$cluster_volume), 
                           max_hole_size = max(hole_data$cluster_volume), 
                           min_hole_size = min(hole_data$cluster_volume),
        total_proportion = sum(hole_data$size) * cellvolume / (nrow(total_int_cells)*cellvolume),
        mean_proportion = mean(hole_data$proportion),
        max_proportion = max(hole_data$proportion),
        min_proportion = min(hole_data$proportion),
        total_absence_proportion = nrow(abs_int)*cellvolume/(nrow(total_int_cells)*cellvolume),
        hypervolume_occupancy = nrow(filled_cells_1) * cellvolume)
      
      hole_size <- rbind(hole_size,h_data)
      
      
      
      
      if(hole_check > (h_data$total_hole_volume * 0.99)| h_data$total_proportion > (h_data$total_absence_proportion * 0.99)){
        break()
      }
      
      hole_check <- h_data$total_hole_volume
   
      }
     
      
      
    }
    
    

    
    ########
    #######
    #### so some metrics to get from these 1) number of holes 2) mean size of holes 3) proportion of internal volume to holes 4)

    
    TPD_holes_list$internal$internal_hole_cells <- abs_int[, c(1:3)]
    TPD_holes_list$internal$metrics_frame <- hole_size
    
    
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
    trophic_niche_colours$Om <- "orchid"
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
      data_3d[["pl_dat"]] %>% dplyr::group_by(x, y, z) %>%
      dplyr::mutate(prob = ifelse(prob == 0, NA, prob)) %>% filter(!is.na(prob)) %>% data.frame()
    
    
    if (!is.null(randata)) {
      ran_data_3d <- TPD_plot_data(randata, sites)
      
      ran_filled_cells <-
        ran_data_3d[["pl_dat"]] %>% dplyr::group_by(x, y, z) %>%
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
        filled_cells %>% merge(species_tpd, by = c("y", "x", "z"))
      
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
    
    
    x <- data$x
    y <- data$y
    z <- data$z
    
    
    
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
          rpm = 6,
          dev = cur3d()
        ),
        startTime = 0,
        duration = 10,
        dir = ".",
        type = "gif",
        clean = T,
        fps = 12
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
                                   standardise = FALSE,
                                   threshold = 1) {
  
  if (!all(sites %in% names(data)) | !all(sites %in% names(randata))) {
    sites <- sites[which(sites %in% names(data))]
    sites <- sites[which(sites %in% names(randata))]
  }
  
  
  prop_occ_df <- c()
  
  sites1_data <- TPD_plot_data(data, sites)[["pl_dat"]] %>% percentile_cells() %>% 
    dplyr::mutate(prob = ifelse(percentile >= threshold, 0, prob))
  ransites_data <- TPD_plot_data(randata, sites)[["pl_dat"]] %>% percentile_cells() %>%
    dplyr::mutate(prob = ifelse(percentile >= threshold, 0, prob))
  
  
  
  for (g in guilds) {
    ran_guild_data <- ranfordata %>% dplyr::filter(Trophic_niche == g)
    
    guild_cells <- c()
        for (i in 1:nrow(ran_guild_data)) {
      cell <-
        as.numeric(
          which(
            sites1_data$y == ran_guild_data[i, "y"] &
              sites1_data$x == ran_guild_data[i, "x"] &
              sites1_data$z == ran_guild_data[i, "z"]
          )
        )
      
      
      guild_cells <- c(guild_cells, cell)
          }
      
    TPD_i <- sites1_data[guild_cells, "prob"]
    
    TPD_j <- ransites_data %>%
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
    "Primary forest/Secondary vegetation",
    "Primary forest/Primary non-forest",
    "Primary forest/Plantation forest",
    "Primary forest/Pasture",
    "Primary forest/Cropland",
    "Primary forest/Minimal agriculture",
    "Primary forest/Intensive agriculture",
    "Primary forest/Urban",
    "Primary non-forest/Secondary vegetation",
    "Primary non-forest/Plantation forest",
    "Primary non-forest/Pasture",
    "Primary non-forest/Cropland",
    "Primary non-forest/Urban",
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
# 
# data = PREDICTS_tpds
# sites = lu_sites
# tpd_for_map = TPD_map
# threshold = seq(0.1, 1, 0.1)


centre_of_mass_for <-
  function(data, sites = NULL, tpd_for_map, percentile) {
    sites_data <- TPD_plot_data(data, sites)
    
    
    prop_COM <- c()
    COM_coords <- c()
    
    
    functional_volume <-
      sites_data[["pl_dat"]] %>% dplyr::group_by(x, y, z) %>%
      dplyr::mutate(prob = ifelse(prob == 0, NA, prob))  %>% data.frame() %>% percentile_cells() %>% dplyr::filter(!is.na(percentile))
    
    full_guild_prop <- functional_volume %>%
      merge(tpd_for_map[, colnames(tpd_for_map) != "prob"], by = c("y", "x", "z")) %>% dplyr::group_by(Trophic_niche) %>%
      dplyr::mutate(
        full_probability = sum(prob),
        n_cells = n(),
        total_cells = nrow(functional_volume),
        full_prop_COM = n_cells / total_cells
      ) %>%
      dplyr::distinct(Trophic_niche, full_prop_COM, full_probability)
    
    
    for (i in 1:length(percentile)) {
      percentile_vol <-
        functional_volume %>% dplyr::filter(percentile <= percentile[i]) %>%
        dplyr::select(y, x, z, prob) %>% dplyr::mutate(percentile = percentile[i])
      
      if (nrow(percentile_vol) == 0) {
        next()
      }
      
      
      percentile_vol$prob <-
        percentile_vol$prob / sum(percentile_vol$prob)
      
      weights <- percentile_vol %>% pull(prob)
      
      
      
      
      
      prop_species <-
        tpd_for_map[, colnames(tpd_for_map) != "prob"] %>% merge(percentile_vol, by = c("y", "x", "z")) %>%
        dplyr::group_by(Trophic_niche) %>%
        dplyr::summarise(n_cells = n()) %>% ungroup() %>% dplyr::mutate(prop_COM = n_cells /
                                                                          sum(n_cells),
                                                                        percentile = percentile[i]) %>%
        merge(full_guild_prop) %>% dplyr::mutate(relative_prop_COM = prop_COM /
                                                   full_prop_COM) %>%
        dplyr::select(-n_cells) %>%
        data.frame()
      
      prob_species <-
        tpd_for_map[, colnames(tpd_for_map) != "prob"] %>% merge(percentile_vol, by = c("y", "x", "z")) %>%
        group_by(Trophic_niche) %>%
        dplyr::summarise(prob_COM = sum(prob), percentile = percentile[i]) %>%
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
        c(apply(percentile_vol[, c("y", "x", "z")], 2, weighted.mean, w = weights)) /
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
