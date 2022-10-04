### densityProfile: function to estimate functional richness profiles (functional space occupied by probabilistic quantiles) from kde object
#### (x is a 'kde' object)
densityProfile <- function(x, probs=seq(0, 0.99, by=0.01)){
  TPD <- x$estimate
  if(class(x$eval.points) == "list"){
    x$eval.points <- expand.grid(x$eval.points)
  }
  cellEdge <- numeric()
  for(i in 1:ncol(x$eval.points)){
    cellEdge[i] <- unique(x$eval.points[,i])[2] - unique(x$eval.points[,i])[1]
  }
  cellSize <- prod(cellEdge)
  TPD <- TPD/sum(TPD)
  alphaSpace_aux <- TPD[order(TPD, decreasing = T)]
  FRicFunct <- function(TPD, alpha = 0.99) { ### FUNCTION TO ESTIMATE FUNCTIONAL RICHNESS FROM KDE OBJECT
    FRic <- numeric(length(alpha))
    for(i in 1:length(alpha)){
      TPDAux <- TPD
      greater_prob <- max(alphaSpace_aux[which(cumsum(alphaSpace_aux) > alpha[i])])
      TPDAux[TPDAux < greater_prob] <- 0
      FRic[i] <- sum(TPDAux>0)* cellSize
    }
    names(FRic) <- alpha
    return(FRic)
  }
  result <- FRicFunct(TPD = TPD, alpha = probs)
  return(result)
}

### densityProfileTPD: function to estimate functional richness profiles (functional space occupied by probabilistic quantiles) from TPDc object
#### (x is a 'TPDc' object)
densityProfileTPD <- function(x, probs=seq(0, 1, by=0.01)){
  TPDList <- x$TPDc$TPDc
  results <- matrix(NA, nrow=length(TPDList), ncol=length(probs),
                    dimnames = list(names(TPDList), probs))
  for(comm in 1:length(TPDList)){
    TPD <- TPDList[[comm]]
    cellSize <- x$data$cell_volume
    alphaSpace_aux <- TPD[order(TPD, decreasing = T)]
    FRicFunct <- function(TPD, alpha = 1) { ### FUNCTION TO ESTIMATE FUNCTIONAL RICHNESS FROM TPD OBJECT
      FRic <- numeric(length(alpha))
      for(i in 1:length(alpha)){
        TPDAux <- TPD
        greater_prob <- max(alphaSpace_aux[which(cumsum(alphaSpace_aux) >= alpha[i])])
        TPDAux[TPDAux < greater_prob] <- 0
        FRic[i] <- sum(TPDAux>0)* cellSize
      }
      names(FRic) <- alpha
      return(FRic)
    }
    results[comm,] <- FRicFunct(TPD = TPD, alpha = probs)  
  }
  return(results)
}
### densityProfileTPD_large: function to estimate (functional space occupied by probabilistic quantiles) from TPDc_large object
#### x is a TPDc object created with the "large" version
densityProfileTPD_large <- function(x, probs=seq(0, 1, by=0.01)){
  TPDList <- x$TPDc$TPDc
  results <- matrix(NA, nrow=length(TPDList), ncol=length(probs),
                    dimnames = list(names(TPDList), probs))
  for(comm in 1:length(TPDList)){
    TPD <- TPDList[[comm]][,"notZeroProb"]
    cellSize <- x$data$cell_volume
    alphaSpace_aux <- TPD[order(TPD, decreasing = T)]
    FRicFunct <- function(TPD, alpha = 1) { ### FUNCTION TO ESTIMATE FUNCTIONAL RICHNESS FROM TPD OBJECT
      FRic <- numeric(length(alpha))
      for(i in 1:length(alpha)){
        TPDAux <- TPD
        greater_prob <- max(alphaSpace_aux[which(cumsum(alphaSpace_aux) >= alpha[i])])
        TPDAux[TPDAux < greater_prob] <- 0
        FRic[i] <- sum(TPDAux>0)* cellSize
      }
      names(FRic) <- alpha
      return(FRic)
    }
    results[comm,] <- FRicFunct(TPD = TPD, alpha = probs)  
  }
  return(results)
}

### overlapF: function to estimate probabilistic overlap between kde objects
overlapF <- function(x,y){
  if(!identical(x$eval.points, y$eval.points)) stop("Not identical eval points!")
  TPDx <- x$estimate
  if(class(x$eval.points) == "list"){
    x$eval.points <- expand.grid(x$eval.points)
  }
  cellEdge <- numeric()
  for(i in 1:ncol(x$eval.points)){
    cellEdge[i] <- unique(x$eval.points[,i])[2] - unique(x$eval.points[,i])[1]
  }
  cellSizex <- prod(cellEdge)
  TPDx <- TPDx * cellSizex
  TPDx <- TPDx/sum(TPDx)
  
  TPDy<- y$estimate
  TPDy <- TPDy * cellSizex
  TPDy <- TPDy/sum(TPDy)
  OL <- sum(pmin(TPDx, TPDy))
  return(OL)
}

### percentileTPD: function to transform probabilities from TPDc object into quantiles
### x is a TPDc object 
percentileTPD <- function(x){
  TPDList <- x$TPDc$TPDc
  results <- list()
  for(comm in 1:length(TPDList)){
    percentile <- rep(NA, nrow(TPDList[[comm]]))
    TPDList[[comm]]<- cbind(TPDList[[comm]], percentile)
    orderTPD <- order(TPDList[[comm]][,"notZeroProb"], decreasing = T)
    TPDList[[comm]] <- TPDList[[comm]][orderTPD,]
    TPDList[[comm]][,"percentile"] <- cumsum(TPDList[[comm]][,"notZeroProb"])
    TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][,"notZeroIndex"]),]
    results[[comm]] <- TPDList[[comm]]
  }
  names(results) <- names(TPDList)
  return(results)
}

### imageTPD: fucntion to plot a TPDc object of two dimensions
# x is a TPDc object 
imageTPD <- function(x, thresholdPlot = 0.99){
  TPDList <- x$TPDc$TPDc
  imageTPD <- list()
  for(comm in 1:length(TPDList)){
    percentile <- rep(NA, length(TPDList[[comm]]))
    TPDList[[comm]]<- cbind(index = 1:length(TPDList[[comm]]),
                            prob = TPDList[[comm]], percentile)
    orderTPD <- order(TPDList[[comm]][,"prob"], decreasing = T)
    TPDList[[comm]] <- TPDList[[comm]][orderTPD,]
    TPDList[[comm]][,"percentile"] <- cumsum(TPDList[[comm]][,"prob"])
    TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][,"index"]),]
    imageTPD[[comm]] <- TPDList[[comm]]
  }
  names(imageTPD) <- names(TPDList)
  spacePercentiles <- matrix(data=0, nrow = nrow(x$data$evaluation_grid), 
                             ncol = length(TPDList), 
                             dimnames = list(1:nrow(x$data$evaluation_grid),
                                             names(TPDList)))
  trait1Edges <- unique(x$data$evaluation_grid[,1])
  trait2Edges <- unique(x$data$evaluation_grid[,2])
  imageMat <- array(NA, c(length(trait1Edges), 
                          length(trait2Edges),
                          length(imageTPD)),
                    dimnames = list(trait1Edges, trait2Edges, names(TPDList)))
  for(comm in 1:length(TPDList)){
    percentileSpace <- x$data$evaluation_grid
    percentileSpace$percentile <- rep(NA, nrow(percentileSpace))
    percentileSpace[, "percentile"] <-imageTPD[[comm]][,"percentile"]
    for(i in 1:length(trait2Edges)){
      colAux <- subset(percentileSpace, percentileSpace[,2] == trait2Edges[i])  
      imageMat[, i, comm] <- colAux$percentile 
    } 
    imageMat[imageMat > thresholdPlot] <- NA
  }
  return(imageMat)
}

### quantileTPD: function to transform probabilities from TPDc object into quantiles
### x is a TPDc object 
quantileTPD <- function(x, thresholdPlot = 0.99){
  TPDList <- x$TPDc$TPDc
  quantileTPD <- list()
  for(comm in 1:length(TPDList)){
    percentile <- rep(NA, length(TPDList[[comm]]))
    TPDList[[comm]]<- cbind(index = 1:length(TPDList[[comm]]),
                            prob = TPDList[[comm]], percentile)
    orderTPD <- order(TPDList[[comm]][,"prob"], decreasing = T)
    TPDList[[comm]] <- TPDList[[comm]][orderTPD,]
    TPDList[[comm]][,"percentile"] <- cumsum(TPDList[[comm]][,"prob"])
    TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][,"index"]),]
    quantileTPD[[comm]] <- TPDList[[comm]]
  }
  names(quantileTPD) <- names(TPDList)
  spacePercentiles <- matrix(data=0, nrow = nrow(x$data$evaluation_grid), 
                             ncol = length(TPDList), 
                             dimnames = list(1:nrow(x$data$evaluation_grid),
                                             names(TPDList)))
  trait1Edges <- unique(x$data$evaluation_grid[,1])
  trait2Edges <- unique(x$data$evaluation_grid[,2])
  imageMat <- array(NA, c(length(trait1Edges), 
                          length(trait2Edges),
                          length(quantileTPD)),
                    dimnames = list(trait1Edges, trait2Edges, names(TPDList)))
  for(comm in 1:length(TPDList)){
    percentileSpace <- x$data$evaluation_grid
    percentileSpace$percentile <- rep(NA, nrow(percentileSpace))
    percentileSpace[, "percentile"] <-quantileTPD[[comm]][,"percentile"]
    for(i in 1:length(trait2Edges)){
      colAux <- subset(percentileSpace, percentileSpace[,2] == trait2Edges[i])  
      imageMat[, i, comm] <- 1 - colAux$percentile 
    } 
    imageMat[imageMat < (1 - thresholdPlot)] <- 0
  }
  return(imageMat)
}

### MAE_TPD: function to estimate differences of quantiles between two quantileTPD objects
MAE_TPD <- function(x, y){# x and y are made with quantileTPD function
  diffXY <- abs(x - y)
  MAE <- mean(diffXY)
  return(MAE)
}  

### quantileTPD_large: function to transform probabilities from TPDc_large object into quantiles
### x is a TPDc_large object 
quantileTPD_large <- function(x){
  TPDList <- x$TPDc$TPDc
  results <- list()
  for(comm in 1:length(TPDList)){
    percentile <- rep(NA, nrow(TPDList[[comm]]))
    TPDList[[comm]]<- cbind(TPDList[[comm]], percentile)
    orderTPD <- order(TPDList[[comm]][,"notZeroProb"], decreasing = T)
    TPDList[[comm]] <- TPDList[[comm]][orderTPD,]
    TPDList[[comm]][,"percentile"] <- 1- cumsum(TPDList[[comm]][,"notZeroProb"])
    TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][,"notZeroIndex"]),]
    results[[comm]] <- TPDList[[comm]]
  }
  names(results) <- names(TPDList)
  return(results)
}

### MAE_TPD_large: function to estimate differences of quantiles between two quantileTPD_large objects
MAE_TPD_large <- function(x, y){# x and y are made with quantileTPD_large function
  mergeXY <- merge(x, y, by="notZeroIndex", all=T)
  mergeXY[is.na(mergeXY)] <- 0
  mergeXY$change <- abs(mergeXY$percentile.x - mergeXY$percentile.y)
  MAE <- mean(mergeXY$change)
  return(MAE)
}

### TPDsMean_large: function to estimate TPDs functions using a single average value per species and a given bandwidth (standard deviation). This function is equivalent to TPD::TPDsMean, but does not record the cells with zero probability, making it more suitable for higher dimensions.
TPDsMean_large<- function(species, means, sds, alpha = 0.95, samples = NULL,
                          trait_ranges = NULL, n_divisions = NULL, tolerance = 0.05) {
  
  # INITIAL CHECKS:
  #   1. Compute the number of dimensions (traits):
  means <- as.matrix(means)
  dimensions <- ncol(means)
  if (dimensions > 4) {
    stop("No more than 4 dimensions are supported at this time; reduce the
         number of dimensions")
  }
  #   2. sds and means must have the same dimensions:
  sds <- as.matrix(sds)
  if (all(dim(means) != dim(sds))) {
    stop("'means' and 'sds' must have the same dimensions")
  }
  #   3. species and means must have the same "dimensions":
  if (length(species) != nrow(means)) {
    stop("The length of 'species' does not match the number of rows of 'means'
         and 'sds'")
  }
  #	4. NA's not allowed in means, sds & species:
  if (any(is.na(means)) | any(is.na(sds)) | any(is.na(species))) {
    stop("NA values are not allowed in 'means', 'sds' or 'species'")
  }
  #	5. Compute the species or populations upon which calculations will be done:
  if (is.null(samples)) {
    species_base <- species
    if (length(unique(species_base)) == 1){
      type <- "One population_One species"
    } else{
      type <- "One population_Multiple species"
    }
  } else {
    if (length(samples) != nrow(means)) {
      stop("The length of 'samples' does not match the number of rows of 'means'
         and 'sds'")
    }
    if (any(is.na(samples))) {
      stop("NA values are not allowed in 'samples'")
    }
    species_base <- paste(species, samples, sep = ".")
    if (length(unique(species)) == 1){
      type <- "Multiple populations_One species"
    } else{
      type <- "Multiple populations_Multiple species"
    }
  }
  
  #	6. Define trait ranges:
  if (is.null(trait_ranges)) {
    trait_ranges <- rep (5, dimensions)
  }
  if (class(trait_ranges) != "list") {
    trait_ranges_aux <- trait_ranges
    trait_ranges <- list()
    for (dimens in 1:dimensions) {
      max_aux <- max(means[, dimens] + trait_ranges_aux[dimens] * sds[, dimens])
      min_aux <- min(means[, dimens] - trait_ranges_aux[dimens] * sds[, dimens])
      trait_ranges[[dimens]] <- c(min_aux, max_aux)
    }
  }
  #	6. Create the grid of cells in which the density function is evaluated:
  if (is.null(n_divisions)) {
    n_divisions_choose<- c(1000, 200, 50, 25)
    n_divisions<- n_divisions_choose[dimensions]
  }
  grid_evaluate<-list()
  edge_length <- list()
  cell_volume<-1
  for (dimens in 1:dimensions){
    grid_evaluate[[dimens]] <- seq(from = trait_ranges[[dimens]][1],
                                   to = trait_ranges[[dimens]][2],
                                   length=n_divisions)
    edge_length[[dimens]] <- grid_evaluate[[dimens]][2] -
      grid_evaluate[[dimens]][1]
    cell_volume <- cell_volume * edge_length[[dimens]]
  }
  evaluation_grid <- expand.grid(grid_evaluate)
  if (is.null(colnames(means))){
    names(evaluation_grid) <- paste0("Trait.",1:dimensions)
  } else {
    names(evaluation_grid) <- colnames(means)
  }
  if (dimensions==1){
    evaluation_grid <- as.matrix(evaluation_grid)
  }
  # Creation of lists to store results:
  results <- list()
  # DATA: To store data and common features
  results$data <- list()
  results$data$evaluation_grid <- evaluation_grid
  results$data$cell_volume <- cell_volume
  results$data$edge_length <- edge_length
  results$data$species <- species
  results$data$means <- means
  results$data$sds <- sds
  if (is.null(samples)){
    results$data$populations <-  NA
  } else{
    results$data$populations <-  species_base
  }
  
  results$data$alpha <- alpha
  results$data$pop_means <- list()
  results$data$pop_sds <- list()
  results$data$pop_sigma <- list()
  results$data$dimensions <- dimensions
  results$data$type <- type
  results$data$method <- "mean"
  
  # TPDs: To store TPDs features of each species/population
  results$TPDs<-list()
  
  
  ########Multivariate normal density calculation
  for (spi in 1:length(unique(species_base))) {
    # Some information messages
    if (spi == 1) { message(paste0("-------Calculating densities for ", type, "-----------\n")) }
    #Data selection
    selected_rows <- which(species_base == unique(species_base)[spi])
    results$data$pop_means[[spi]] <- means[selected_rows, ]
    results$data$pop_sds[[spi]] <- sds[selected_rows, ]
    names(results$data$pop_means)[spi] <- names(results$data$pop_sds)[spi]<-
      unique(species_base)[spi]
    if (dimensions > 1) {
      results$data$pop_sigma[[spi]] <- diag(results$data$pop_sds[[spi]]^2)
      
      multNormAux <- mvtnorm::dmvnorm(x = evaluation_grid,
                                      mean = results$data$pop_means[[spi]],
                                      sigma = results$data$pop_sigma[[spi]])
      multNormAux <- multNormAux / sum(multNormAux)
      # Now, we extract the selected fraction of volume (alpha), if necessary
      extract_alpha <- function(x){
        # 1. Order the 'cells' according to their probabilities:
        alphaSpace_aux <- x[order(x, decreasing=T)]
        # 2. Find where does the accumulated sum of ordered probabilities becomes
        #   greater than the threshold (alpha):
        greater_prob <- alphaSpace_aux[which(cumsum(alphaSpace_aux ) > alpha) [1]]
        # 3. Substitute smaller probabilities with 0:
        x[x < greater_prob] <- 0
        # 5. Finally, reescale, so that the total density is finally 1:
        x <- x / sum(x)
        return(x)
      }
      if (alpha < 1){
        multNormAux <- extract_alpha(multNormAux)
      }
      notZeroIndex <- which(multNormAux != 0)
      notZeroProb <- multNormAux[notZeroIndex]
      results$TPDs[[spi]] <- cbind(notZeroIndex, notZeroProb)
      
      
    }
    if (dimensions == 1) stop("This function is intended for > 1 dimension")
  }
  names(results$TPDs) <- unique(species_base)
  class(results) <- "TPDsp"
  return(results)
}

### TPDc_large: function to estimate Trait Probability Density of Communities based on TPDs_large objects
TPDc_large <- function(TPDs, sampUnit){
  sampUnit <- as.matrix(sampUnit)
  # 1. species names:
  if (is.null(colnames(sampUnit)) | any(is.na(colnames(sampUnit)))) {
    stop("colnames(sampUnit), must contain the names of the species.
      NA values are not allowed")
  }
  # 2. communities names:
  if (is.null(rownames(sampUnit)) | any(is.na(rownames(sampUnit)))){
    stop("rownames(sampUnit), must contain the names of the sampling units.
      NA values are not allowed")
  }
  # 3. Data values will be inherithed from TPDs, which must be of class TPD
  if (class(TPDs) != "TPDsp") {
    stop("TPDs must be an object of class 'TPDsp', created with the function
      'TPDs'")
  }
  species <- samples <- abundances <- numeric()
  for (i in 1:nrow(sampUnit)){
    samples <- c(samples, rep(rownames(sampUnit)[i], ncol(sampUnit)))
    species <- c(species, colnames(sampUnit))
    abundances <- c(abundances, sampUnit[i, ])
  }
  nonZero <- which(abundances > 0)
  samples <- samples[nonZero]
  species <- species[nonZero]
  abundances <- abundances[nonZero]
  
  # Creation of lists to store results:
  results <- list()
  results$data <- TPDs$data
  results$data$sampUnit <- sampUnit
  type <- results$data$type
  # All the species or populations in 'species' must be in the species or
  #   populations of TPDs:
  if (type == "Multiple populations_One species" |
      type == "Multiple populations_Multiple species") {
    species_base <- paste(species, samples, sep = ".")
    if (!all(unique(species_base) %in% unique(results$data$populations))) {
      non_found_pops <- which(unique(species_base) %in%
                                unique(results$data$populations) == 0)
      stop("All the populations TPDs must be present in 'TPDs'. Not present:\n",
           paste(species_base[non_found_pops], collapse=" / "))
    }
  }
  if (type == "One population_One species" |
      type == "One population_Multiple species") {
    species_base <- species
    if (!all(unique(species_base) %in% unique(results$data$species))) {
      non_found_sps <- which(unique(species_base) %in%
                               unique(results$data$species) == 0)
      stop("All the species TPDs must be present in 'TPDs'. Not present:\n",
           paste(species_base[non_found_sps], collapse=" / "))
    }
  }
  # END OF INITIAL CHECKS
  # TPDc computation
  results$TPDc <- list()
  results$TPDc$species <- list()
  results$TPDc$abundances <- list()
  results$TPDc$speciesPerCell <- list()
  # results$TPDc$RTPDs <- list()
  results$TPDc$TPDc <- list()
  
  for (samp in 1:length(unique(samples))) {
    selected_rows <- which(samples == unique(samples)[samp])
    species_aux <- species_base[selected_rows]
    abundances_aux <- abundances[selected_rows] / sum(abundances[selected_rows])
    RTPDsAux <- rep(0, nrow(results$data$evaluation_grid))
    TPDs_aux <- TPDs$TPDs[names(TPDs$TPDs) %in% species_aux]
    cellsOcc <- numeric()
    for (sp in 1:length(TPDs_aux)) {
      selected_name <- which(names(TPDs_aux) == species_aux[sp])
      cellsToFill <- TPDs_aux[[selected_name]][,"notZeroIndex"]
      cellsOcc <- c(cellsOcc, cellsToFill)
      probsToFill <- TPDs_aux[[selected_name]][,"notZeroProb"] * abundances_aux[sp]
      RTPDsAux[cellsToFill] <- RTPDsAux[cellsToFill] + probsToFill
    }
    TPDc_aux <- RTPDsAux
    notZeroIndex <- which(TPDc_aux != 0)
    notZeroProb <- TPDc_aux[notZeroIndex]
    results$TPDc$TPDc[[samp]] <- cbind(notZeroIndex, notZeroProb)
    results$TPDc$species[[samp]] <- species_aux
    results$TPDc$abundances[[samp]] <- abundances_aux
    results$TPDc$speciesPerCell[[samp]] <- table(cellsOcc)
    names(results$TPDc$TPDc)[samp] <-
      names(results$TPDc$species)[samp] <- names(results$TPDc$abundances)[samp] <-
      names(results$TPDc$speciesPerCell)[samp] <- unique(samples)[samp]
  }
  class(results) <- "TPDcomm"
  return(results)
}

### TPDRichness: Function to estimate functional richness from TPD object
TPDRichness <- function(TPDc = NULL, TPDs = NULL){
  # INITIAL CHECKS:
  # 1. At least one of TPDc or TPDs must be supplied.
  if (is.null(TPDc) & is.null(TPDs)) {
    stop("At least one of 'TPDc' or 'TPDs' must be supplied")
  }
  if (!is.null(TPDc) & class(TPDc) != "TPDcomm"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  if (!is.null(TPDs) & class(TPDs) != "TPDsp"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  # Creation of lists to store results:
  results <- list()
  # 1. Functional Richness
  Calc_FRich <- function(x) {
    results_FR <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux <- TPD[[i]]
      TPD_aux[TPD_aux > 0] <- cell_volume
      results_FR[i] <- sum(TPD_aux)
    }
    names(results_FR) <- names_aux
    return(results_FR)
  }
  
  # IMPLEMNENTATION
  if (!is.null(TPDc)) {
    results$communities <- list()
    # message("Calculating FRichness of communities")
    results$communities$FRichness <- Calc_FRich(TPDc)
  }
  if (!is.null(TPDs)) {
    if (TPDs$data$type == "One population_One species" |
        TPDs$data$type == "One population_Multiple species") {
      results$species <- list()
      # message("Calculating FRichness of species")
      results$species$FRichness <- Calc_FRich(TPDs)
    } else {
      results$populations <- list()
      # message("Calculating FRichness of populations")
      results$populations$FRichness <- Calc_FRich(TPDs)
    }
    if (TPDs$data$method == "mean") {
      message("WARNING: When TPDs are calculated using the TPDsMean function, Evenness
              and Divergence are meaningless!!")
    }
  }
  return(results)
}

### TPDRichness_large: Function to estimate functional richness from TPD_large object
TPDRichness_large <- function(TPDc = NULL, TPDs = NULL){
  # INITIAL CHECKS:
  # 1. At least one of TPDc or TPDs must be supplied.
  if (is.null(TPDc) & is.null(TPDs)) {
    stop("At least one of 'TPDc' or 'TPDs' must be supplied")
  }
  if (!is.null(TPDc) & class(TPDc) != "TPDcomm"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  if (!is.null(TPDs) & class(TPDs) != "TPDsp"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  # Creation of lists to store results:
  results <- list()
  # 1. Functional Richness
  Calc_FRich <- function(x) {
    results_FR <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux <- TPD[[i]][,"notZeroProb"]
      TPD_aux[TPD_aux > 0] <- cell_volume
      results_FR[i] <- sum(TPD_aux)
    }
    names(results_FR) <- names_aux
    return(results_FR)
  }
  
  # IMPLEMNENTATION
  if (!is.null(TPDc)) {
    results$communities <- list()
    # message("Calculating FRichness of communities")
    results$communities$FRichness <- Calc_FRich(TPDc)
  }
  if (!is.null(TPDs)) {
    if (TPDs$data$type == "One population_One species" |
        TPDs$data$type == "One population_Multiple species") {
      results$species <- list()
      # message("Calculating FRichness of species")
      results$species$FRichness <- Calc_FRich(TPDs)
    } else {
      results$populations <- list()
      # message("Calculating FRichness of populations")
      results$populations$FRichness <- Calc_FRich(TPDs)
    }
    if (TPDs$data$method == "mean") {
      message("WARNING: When TPDs are calculated using the TPDsMean function, Evenness
              and Divergence are meaningless!!")
    }
  }
  return(results)
}

### REND_large: Functional Evenness, Richness and Divergence from TPD_large object
REND_large <- function(TPDc = NULL, TPDs = NULL){
  # INITIAL CHECKS:
  # 1. At least one of TPDc or TPDs must be supplied.
  if (is.null(TPDc) & is.null(TPDs)) {
    stop("At least one of 'TPDc' or 'TPDs' must be supplied")
  }
  if (!is.null(TPDc) & class(TPDc) != "TPDcomm"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  if (!is.null(TPDs) & class(TPDs) != "TPDsp"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  # Creation of lists to store results:
  results <- list()
  # 1. Functional Richness
  Calc_FRich <- function(x) {
    results_FR <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux <- TPD[[i]][,"notZeroProb"]
      TPD_aux[TPD_aux > 0] <- cell_volume
      results_FR[i] <- sum(TPD_aux)
    }
    names(results_FR) <- names_aux
    return(results_FR)
  }
  
  # 2. Functional Evenness
  Calc_FEve <- function(x) {
    results_FE <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux <- TPD[[i]][,"notZeroProb"]
      TPD_eve <- rep((1 / length(TPD_aux)), times = length(TPD_aux))
      results_FE[i] <- sum(pmin(TPD_aux, TPD_eve))
    }
    names(results_FE) <- names_aux
    return(results_FE)
  }
  # 3. Functional Divergence
  Calc_FDiv <- function(x) {
    results_FD <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      evaluation_grid<-x$data$evaluation_grid
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      evaluation_grid<-x$data$evaluation_grid
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) { 
      notZeroCells <- TPD[[i]][,"notZeroIndex"]
      functional_volume <- evaluation_grid[notZeroCells, , drop=F]
      # Functional volume has to be standardised so that distances are
      # independent of the scale of the axes
      for (j in 1:ncol(functional_volume)){
        functional_volume[, j] <-
          (functional_volume[, j] - min(functional_volume[, j])) /
          (max(functional_volume[, j]) - min(functional_volume[, j]))
      }
      TPD_aux <- TPD[[i]][,"notZeroProb"]
      # 1. Calculate the center of gravity
      COG <- colMeans(functional_volume, na.rm=T)
      # 2. Calculate the distance of each point in the functional volume to the
      #   COG:
      dist_COG <- function(x, COG) {
        result_aux<-stats::dist(rbind(x, COG))
        return(result_aux)
      }
      COGDist <- apply(functional_volume, 1, dist_COG, COG)
      # 3. Calculate the mean of the COGDist's
      meanCOGDist <- mean(COGDist)
      # 4. Calculate the sum of the abundance-weighted deviances for distaces
      #   from the COG (AWdistDeviances) and the absolute abundance-weighted
      #   deviances:
      distDeviances <- COGDist - meanCOGDist
      AWdistDeviances <- sum(TPD_aux * distDeviances)
      absdistDeviances <- abs(COGDist - meanCOGDist)
      AWabsdistDeviances <- sum(TPD_aux * absdistDeviances)
      #Finally, calculate FDiv:
      results_FD[i] <- (AWdistDeviances + meanCOGDist) /
        ( AWabsdistDeviances +  meanCOGDist)
    }
    names(results_FD) <- names_aux
    return(results_FD)
  }
  # IMPLEMNENTATION
  if (!is.null(TPDc)) {
    results$communities <- list()
    message("Calculating FRichness of communities")
    results$communities$FRichness <- Calc_FRich(TPDc)
    message("Calculating FEvenness of communities")
    results$communities$FEvenness <- Calc_FEve(TPDc)
    message("Calculating FDivergence of communities")
    results$communities$FDivergence <- Calc_FDiv(TPDc)
  }
  if (!is.null(TPDs)) {
    if (TPDs$data$type == "One population_One species" |
        TPDs$data$type == "One population_Multiple species") {
      results$species <- list()
      message("Calculating FRichness of species")
      results$species$FRichness <- Calc_FRich(TPDs)
      message("Calculating FEvenness of species")
      results$species$FEvenness <- Calc_FEve(TPDs)
      message("Calculating FDivergence of species")
      results$species$FDivergence <- Calc_FDiv(TPDs)
    } else {
      results$populations <- list()
      message("Calculating FRichness of populations")
      results$populations$FRichness <- Calc_FRich(TPDs)
      message("Calculating FEvenness of populations")
      results$populations$FEvenness <- Calc_FEve(TPDs)
      message("Calculating FDivergence of populations")
      results$populations$FDivergence <- Calc_FDiv(TPDs)
    }
    if (TPDs$data$method == "mean") {
      message("WARNING: When TPDs are calculated using the TPDsMean function, Evenness
              and Divergence are meaningless!!")
    }
  }
  return(results)
}

### redundancy_large:  Functional Redundancy of Communities from TPD_large object
redundancy_large <- function(TPDc = NULL) {
  if (class(TPDc) != "TPDcomm") {
    stop("TPDc must be an object of class TPDcomm generated with the TPDc
		    function")
  }
  x <- TPDc
  results <- list()
  results$redundancy <- results$richness <- numeric()
  for (i in 1:length(x$TPDc$TPDc)) {
    TPDc_aux <- x$TPDc$TPDc[[i]][, "notZeroProb"]
    M <- x$TPDc$speciesPerCell[[i]]
    results$redundancy[i] <- sum(M * TPDc_aux) - 1
    results$richness[i] <- sum(x$TPDc$abundances[[i]] >0)
  }
  results$redundancyRelative <- results$redundancy / (results$richness -1)
  names(results$redundancy) <- names(results$richness) <-
    names(results$redundancyRelative) <- names(x$TPDc$TPDc)
  return(results)
}










# 
# ####################
# #################### TRANSFORMING PROBABILITIUES INTO QUANTILES
# 
# 
# # x is a TPDc object created with the "large" version
# PlotPercentileTPD_large <- function(x, whichPlot = 1, ncolors=1000, colorsPick="YlOrBr", 
#                               xlab="Trait 1", ylab="Trait2",
#                               cont=c(50,90), xlim=NULL, ylim=NULL){
#   TPDList <- x$TPDc$TPDc
#   percentilesTPD <- list()
#   for(comm in 1:length(TPDList)){
#     percentile <- rep(NA, nrow(TPDList[[comm]]))
#     TPDList[[comm]]<- cbind(TPDList[[comm]], percentile)
#     orderTPD <- order(TPDList[[comm]][,"notZeroProb"], decreasing = T)
#     TPDList[[comm]] <- TPDList[[comm]][orderTPD,]
#     TPDList[[comm]][,"percentile"] <- cumsum(TPDList[[comm]][,"notZeroProb"]) *100
#     TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][,"notZeroIndex"]),]
#     percentilesTPD[[comm]] <- TPDList[[comm]]
#   }
#   names(percentilesTPD) <- names(TPDList)
#   spacePercentiles <- matrix(data=0, nrow = nrow(x$data$evaluation_grid), 
#                              ncol = length(TPDList), 
#                              dimnames = list(1:nrow(x$data$evaluation_grid),
#                                              names(TPDList)))
#   trait1Edges <- unique(x$data$evaluation_grid[,1])
#   trait2Edges <- unique(x$data$evaluation_grid[,2])
#   imageMat <- array(NA, c(length(trait1Edges), 
#                           length(trait2Edges),
#                           length(percentilesTPD)))
#   ncolors <- 1000
#   ColorRamp <- colorRampPalette(RColorBrewer::brewer.pal(9,colorsPick)[c(9:7, 5:3)])(ncolors)
#   
#   for(comm in whichPlot){
#     percentileSpace <- x$data$evaluation_grid
#     percentileSpace$percentile <- rep(NA, nrow(percentileSpace))
#     percentileSpace[percentilesTPD[[comm]][,"notZeroIndex"], "percentile"] <- 
#       percentilesTPD[[comm]][,"percentile"]
#     for(i in 1:length(trait2Edges)){
#       colAux <- subset(percentileSpace, percentileSpace[,2] == trait2Edges[i])  
#       imageMat[, i, comm] <- colAux$percentile 
#     } 
#     xmin <- 0
#     xmax <- 100
#     ColorLevels<-seq(from=xmin, to=xmax, length=ncolors)
#     ColorRamp_ex <- ColorRamp[round(1+(min(imageMat[, , comm], na.rm=T)-xmin)*ncolors/(xmax-xmin)) : 
#                                 round( (max(imageMat[, , comm], na.rm=T)-xmin)*ncolors/(xmax-xmin) )]
#     if(is.null(xlim)) xlim <- range(trait1Edges)
#     if(is.null(ylim)) ylim <- range(trait2Edges)
#     image(x=trait1Edges, y=trait2Edges, z=imageMat[, , comm], col=ColorRamp_ex,
#           xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=xlim, ylim=ylim,
#           main = names(TPDList)[comm])
#     box(which="plot")
#     axis(1,tcl=0.3,lwd=0.8)
#     axis(2, las=1, tcl=0.3,lwd=0.8)
#     contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , comm], levels=cont,
#             drawlabels = T, labcex = 0.8, lwd=0.5, lty=1, col="black", add=T)
#   }
# }
# 
# # x is a TPDc object created with the "large" version
# PlotPercentileTPD <- function(x, whichPlot = 1, ncolors=1000, colorsPick="YlOrBr", 
#                                     xlab="Trait 1", ylab="Trait2",
#                                     cont=c(0.5,0.9), xlim=NULL, ylim=NULL,
#                               thresholdPlot = 0.99, gradientColorsF=NULL){
#   TPDList <- x$TPDc$TPDc
#   percentilesTPD <- list()
#   for(comm in whichPlot){
#     percentile <- rep(NA, length(TPDList[[comm]]))
#     TPDList[[comm]]<- cbind(index = 1:length(TPDList[[comm]]),
#                             prob = TPDList[[comm]], percentile)
#     orderTPD <- order(TPDList[[comm]][,"prob"], decreasing = T)
#     TPDList[[comm]] <- TPDList[[comm]][orderTPD,]
#     TPDList[[comm]][,"percentile"] <- cumsum(TPDList[[comm]][,"prob"])
#     TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][,"index"]),]
#     percentilesTPD[[comm]] <- TPDList[[comm]]
#   }
#   names(percentilesTPD) <- names(TPDList)[whichPlot]
#   spacePercentiles <- matrix(data=0, nrow = nrow(x$data$evaluation_grid), 
#                              ncol = length(whichPlot), 
#                              dimnames = list(1:nrow(x$data$evaluation_grid),
#                                              names(TPDList)[whichPlot]))
#   trait1Edges <- unique(x$data$evaluation_grid[,1])
#   trait2Edges <- unique(x$data$evaluation_grid[,2])
#   imageMat <- array(NA, c(length(trait1Edges), 
#                           length(trait2Edges),
#                           length(percentilesTPD)))
#   ncolors <- 1000
#   ColorRamp <- rev(gradientColorsF(ncolors))
#   
#   for(comm in whichPlot){
#     percentileSpace <- x$data$evaluation_grid
#     percentileSpace$percentile <- rep(NA, nrow(percentileSpace))
#     percentileSpace[, "percentile"] <-percentilesTPD[[comm]][,"percentile"]
#     for(i in 1:length(trait2Edges)){
#       colAux <- subset(percentileSpace, percentileSpace[,2] == trait2Edges[i])  
#       imageMat[, i, comm] <- colAux$percentile 
#     } 
#     imageMat[imageMat >thresholdPlot] <- NA
#     xmin <- 0
#     xmax <- 1
#     ColorLevels<-seq(from=xmin, to=xmax, length=ncolors)
#     ColorRamp_ex <- ColorRamp[round(1+(min(imageMat[, , comm], na.rm=T)-xmin)*ncolors/(xmax-xmin)) : 
#                                 round( (max(imageMat[, , comm], na.rm=T)-xmin)*ncolors/(xmax-xmin) )]
#     if(is.null(xlim)) xlim <- range(trait1Edges)
#     if(is.null(ylim)) ylim <- range(trait2Edges)
#     image(x=trait1Edges, y=trait2Edges, z=imageMat[, , comm], col=ColorRamp_ex,
#           xaxs="r", yaxs="r", xlab=xlab, ylab=ylab, axes=F, xlim=xlim, ylim=ylim,
#           main = "", asp = 1)
#     box(which="plot")
#     axis(1,tcl=0.3,lwd=0.8)
#     axis(2, las=1, tcl=0.3,lwd=0.8)
#     contour(x=trait1Edges, y=trait2Edges, z=imageMat[, , comm], levels=cont,
#             drawlabels = T, labcex = 0.8, lwd=0.5, lty=1, col="black", add=T)
#   }
# }
# 
# 
# 
# ###################
# ###################
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 