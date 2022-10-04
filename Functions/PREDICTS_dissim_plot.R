PREDICTS_dissim_plot <- function(data){

  
  
axis_col <- function(data,x){
  axis_dat <- na.omit(x)
  for(i in 1:ncol(data)){
    val <- na.omit(data[,i])
    if(all(val == axis_dat)){
      col <- colnames(data)[i]
      break()
    }
  }
  return(col) 
}


# Correlation panel
panel.beta <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  x_name <- axis_col(data[,,1],x)
  y_name <- axis_col(data[,,1],y)
  r <- round(data[,,1][x_name,y_name], digits = 2)
  colours <- colorRampPalette(c("white", "yellow", "orange", "red"))
  plot_cols <- colours(length(seq(0,1,0.01)))
  txt <- paste0("Beta = ", r)
  cex.cor <- 0.8/strwidth(txt)
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = plot_cols[(r*100)+1])
  text(0.5, 0.5, txt, cex = cex.cor, font = 2)
}
# Customize upper panel
panel.nest<- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  x_name <- axis_col(data[,,1],x)
  y_name <- axis_col(data[,,1],y)
  r_2 <- round(data[,,2][x_name,y_name], digits = 2)
  r_3 <- round(data[,,3][x_name,y_name], digits = 2)
  txt_2 <- paste0("Shared = ", r_2)
  txt_3 <- paste0("Not_shared = ", r_3)
  cex.cor_2 <- 0.8/strwidth(txt_2)
  cex.cor_3 <- 0.8/strwidth(txt_3)
  text(0.5, 0.7, txt_2, cex = cex.cor_2+0.1, font = 2)
  text(0.5, 0.3, txt_3, cex = cex.cor_3+0.1, font = 2)
}

col_names <- c("Primary vegetation", "Primary forest","Primary non-forest", "Secondary vegetation", "Plantation forest", "Cropland", "Pasture","Intensive agriculture","Minimal agriculture", "Urban")
diag_colours <- c("olivedrab4","olivedrab4","#95A900", "olivedrab3","lightgreen","gold1","khaki","gold1","khaki","ivory4")
names(diag_colours) <- col_names

diag.pan<- function(x){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  diag_name <- colnames(data[,,1])[which(is.na(x))]
  colour <- diag_colours[diag_name] 
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = colour)
  words <- unlist(strsplit(diag_name, split = " "))
  if(length(words) == 1){
  text(0.5, 0.5, words[1], cex = 1.45, font = 2)
  } else {
    text(0.5, 0.7, words[1], cex = 1.2, font = 2)
    text(0.5, 0.3, words[2], cex = 1.2, font = 2)
  }
    }




pairs2 <- function (x, labels, panel = points, ..., lower.panel = panel, 
            upper.panel = panel, diag.panel = NULL, text.panel = textPanel, 
            label.pos = 0.5 + has.diag/3, cex.labels = NULL, font.labels = 1, 
            row1attop = TRUE, gap = 1) 
  {
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                                 y, txt, cex = cex, font = font)
    localPlot <- function(..., main, oma, font.main, cex.main, yaxt = 'n') plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
    localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
      x <- as.data.frame(x)
      for (i in seq_along(names(x))) {
        if (is.factor(x[[i]]) || is.logical(x[[i]])) 
          x[[i]] <- as.numeric(x[[i]])
        if (!is.numeric(unclass(x[[i]]))) 
          stop("non-numeric argument to 'pairs'")
      }
    }
    else if (!is.numeric(x)) 
      stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
      lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
      upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
      diag.panel <- match.fun(diag.panel)
    if (row1attop) {
      tmp <- lower.panel
      lower.panel <- upper.panel
      upper.panel <- tmp
      tmp <- has.lower
      has.lower <- has.upper
      has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2) 
      stop("only one column in the argument to 'pairs'")
    has.labs <- TRUE
    if (missing(labels)) {
      labels <- colnames(x)
      if (is.null(labels)) 
        labels <- paste("var", 1L:nc)
    }
    else if (is.null(labels)) 
      has.labs <- FALSE
    oma <- if ("oma" %in% nmdots) 
      dots$oma
    else NULL
    main <- if ("main" %in% nmdots) 
      dots$main
    else NULL
    if (is.null(oma)) {
      oma <- c(4, 4, 4, 4)
      if (!is.null(main)) 
        oma[3L] <- 6
    }
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    for (i in if (row1attop) 
      1L:nc
      else nc:1L) for (j in 1L:nc) {
        localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                  type = "n", ...)
        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
          box()
          # edited here...
          #           if (i == 1 && (!(j%%2) || !has.upper || !has.lower)) 
          #           localAxis(1 + 2 * row1attop, x[, j], x[, i], 
          #                       ...)
          # draw x-axi
          mfg <- par("mfg")
          if (i == j) {
            if (has.diag) 
              localDiagPanel(as.vector(x[, i]), ...)
            if (has.labs) {
              par(usr = c(0, 1, 0, 1))
              if (is.null(cex.labels)) {
                l.wid <- strwidth(labels, "user")
                cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
              }
              text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                         font = font.labels)
            }
          }
          else if (i < j) 
            localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                                                           i]), ...)
          else localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                                                              i]), ...)
          if (any(par("mfg") != mfg)) 
            stop("the 'panel' function made a new plot")
        }
        else par(new = FALSE)
      }
    if (!is.null(main)) {
      font.main <- if ("font.main" %in% nmdots) 
        dots$font.main
      else par("font.main")
      cex.main <- if ("cex.main" %in% nmdots) 
        dots$cex.main
      else par("cex.main")
      mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    }
  }

pairs2(data[,,1], upper.panel = panel.nest, lower.panel = panel.beta, diag.panel = diag.pan, labels = NULL, cex.labels = NULL)

}



