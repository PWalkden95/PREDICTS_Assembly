require(gawdis)

rao_diversity_gaw <- function (comm, traits = NULL, phylodist = NULL, checkdata = TRUE, 
          ord = "metric", put.together = NULL, standardize = TRUE, 
          ...) 
{
  diver.internal <- function(community, distance) {
    if (any(is.na(distance))) {
      distance.na <- ifelse(is.na(distance), 0, 1)
      inter.na <- community %*% distance.na
      adjustment <- rowSums(sweep(community, 1, inter.na, 
                                  "*", check.margin = FALSE))
      distance[is.na(distance)] <- 0
      inter <- community %*% distance
      res <- rowSums(sweep(community, 1, inter, "*", 
                           check.margin = FALSE))
      res <- ifelse(adjustment > 0, res/adjustment, res)
    }
    else {
      inter <- community %*% distance
      res <- rowSums(sweep(community, 1, inter, "*", 
                           check.margin = FALSE))
    }
    return(res)
  }
  
  
  
  res <- list(call = match.call())

  list.warning <- list()

  comm <- as.matrix(comm)
  N <- nrow(comm)
  S <- ncol(comm)
  dist.1 <- 1 - diag(x = rep(1, S))
  if (!is.null(traits)) {
    if(class(traits) == "dist"){
      dist.functional <- sqrt(as.matrix(traits)) 
    } else {
    traits <- as.data.frame(traits)
    m <- ncol(traits)
    weights <- rep(1, m)
    make.names <- is.null(colnames(traits))
    colnames(traits) <- colnames(traits, do.NULL = FALSE, 
                                 prefix = "T")
    names(weights) <- colnames(traits)
    
    dist.functional <- sqrt(as.matrix(gawdis(x = traits)))
    }
    }

  comm <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/")
  SD <- diver.internal(comm, dist.1)
  res$Simpson <- SD
  if (!is.null(traits)) {
    FD <- diver.internal(comm, dist.functional)
    res$FunRao <- FD
    res$FunRedundancy <- SD - FD
  }
  if (!is.null(phylodist)) {
    PD <- diver.internal(comm, dist.phylogenetic)
    res$PhyRao <- PD
    res$PhyRedundancy <- SD - PD
  }
  return(res)
}
