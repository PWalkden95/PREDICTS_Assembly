################## OVERFLOW FUNCTIONS

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
