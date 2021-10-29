

Calc_FRich <- function(data){
  cell_volume <- data[["data"]][["cell_volume"]]
  
  site_data <- data[-1]
  
FRich <- data.frame(SSBS = names(site_data),FRich = NA)
j <- 1
for(i in 1:length(site_data)){
  TPD_aux <- site_data[[i]][["TPDc"]][["RelativeAbundance"]]
  TPD_aux[TPD_aux > 0] <- cell_volume
  FRich[j,"FRich"] <- sum(TPD_aux)
  j <- j + 1
}
return(FRich)
}




################################################
####### Functional Evenness calc
Calc_FEve <- function(data) {
  
  site_data <- data[-1]
  j <- 1
  FEve <- data.frame(SSBS = names(site_data),FEve = NA)
  for (i in 1:length(site_data)) {
    TPD <- data[[i]][["TPDc"]][["RelativeAbundance"]]
    TPD_aux <- TPD[TPD > 0]
    TPD_eve <- rep((1/length(TPD_aux)), times = length(TPD_aux))
    FEve[j,"FEve"] <- sum(pmin(TPD_aux, TPD_eve))
    j <- j +1
  }
  return(FEve)
}


##############################################################
#############################################################
# Functional Divergence
Calc_FDiv <- function(data) {
  evaluation_grid <- data$data$evaluation_grid
  cell_volume <- data$data$cell_volume
  
  site_data <- data[-1]
  
  FDiv <- data.frame(SSBS = names(site_data),FDiv = NA)
  k <- 1
  for (i in 1:length(site_data)) {
    TPD <- site_data[[i]][["TPDc"]][["RelativeAbundance"]]
    functional_volume <- evaluation_grid[TPD > 0,
                                         , drop = F]
    for (j in 1:ncol(functional_volume)) {
      functional_volume[, j] <- (functional_volume[,
                                                   j] - min(functional_volume[, j]))/(max(functional_volume[,
                                                                                                            j]) - min(functional_volume[, j]))
    }
    TPD_aux <- TPD[TPD > 0]
    COG <- colMeans(functional_volume, na.rm = T)
    dist_COG <- function(x, COG) {
      result_aux <- stats::dist(rbind(x, COG))
      return(result_aux)
    }
    COGDist <- apply(functional_volume, 1, dist_COG,
                     COG)
    meanCOGDist <- mean(COGDist)
    distDeviances <- COGDist - meanCOGDist
    AWdistDeviances <- sum(TPD_aux * distDeviances)
    absdistDeviances <- abs(COGDist - meanCOGDist)
    AWabsdistDeviances <- sum(TPD_aux * absdistDeviances)
    FDiv[k,"FDiv"] <- (AWdistDeviances + meanCOGDist)/(AWabsdistDeviances +
                                                         meanCOGDist)
    k <- k + 1
  }
  return(FDiv)
}


TPD_FD_metrics <- function(data){

  FRich <- Calc_FRich(data)
  FEve <- Calc_FEve(data)
  FDiv <- Calc_FDiv(data)
  
  
  
FD_metrics <- FRich %>% dplyr::left_join(FEve, by = "SSBS") %>% dplyr::left_join(FDiv, by = "SSBS")

return(FD_metrics)
}



FD_metrics <- TPD_FD_metrics(PREDICTS_tpds)

############################################
############################################
obj_2_string <-function(x){
str <- deparse(substitute(x))
return(str)
}


Calc_dissim <- function(data,sites1,sites2){


  
results_samp <- list()
results_samp$dissim$dissimilarity <- NA
results_samp$dissim$P_shared <- NA
results_samp$dissim$P_non_shared <- NA

sites1_data <- TPD_plot_data(data, sites1)
sites2_data <- TPD_plot_data(data, sites2)


TPD_i <- sites1_data[["pl_dat"]][["prob"]]
TPD_j <- sites2_data[["pl_dat"]][["prob"]]



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
results_samp$dissim$dissimilarity <- results_samp$dissim$dissimilarity <- 1 - O_aux
results_samp$dissim$P_non_shared <- results_samp$dissim$P_non_shared     <- (2 * min(B_aux, C_aux))/(A_aux +
                                                                                         2 * min(B_aux, C_aux))
results_samp$dissim$P_shared <- results_samp$dissim$P_shared  <- 1 - results_samp$dissim$P_non_shared



return(results_samp)

}



test <- Calc_dissim(PREDICTS_tpds, sites1 = primary_for, sites2 = primary_non_for)



rm(list = ls())
require(tidyverse)
PREDICTS_tpds <- readRDS("Outputs/PREDICTS_sites_tpds.rds")
randomisations_tpds <- readRDS("Outputs/randomisations_TPD_morpho.rds")
PREDICTS_tpds_for <- readRDS("Outputs/PREDICTS_sites_for_tpds.rds")
PREDICTS <- readRDS("Outputs/refined_predicts.rds") %>% dplyr::distinct(SSBS, Predominant_habitat, Use_intensity, Biome, UN_subregion, Realm) %>%
  dplyr::mutate(Predominant_habitat = ifelse(grepl(Predominant_habitat, pattern = "secondary", ignore.case = TRUE), "Secondary vegetation",
                                             paste(Predominant_habitat)))
TPD_LU <- data.frame(SSBS = names(PREDICTS_tpds)) %>% dplyr::left_join(PREDICTS)
source("Functions/TPD_3D_plots.R")
land_use <- PREDICTS %>% dplyr::distinct(Predominant_habitat, .keep_all = FALSE) %>% dplyr::filter(Predominant_habitat != "Cannot decide") %>%
  pull() %>% as.character()
table(TPD_LU$Predominant_habitat, TPD_LU$UN_subregion)
table(TPD_LU$Predominant_habitat, TPD_LU$Biome)
table(TPD_LU$Predominant_habitat, TPD_LU$Realm)
primary_for <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[5], Realm == "Australasia") %>% pull(SSBS) %>% as.character()
primary_non_for <- TPD_LU %>% dplyr::filter(Predominant_habitat == land_use[6], Realm == "Australasia") %>% pull(SSBS) %>% as.character()
TPD_Diff_plot(data = PREDICTS_tpds, sites1 = primary_for, sites2 = secondary, T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body",
              method = "prob")
TPD_3d_plot(PREDICTS_tpds, sites = primary_non_for,  T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body", method = "prob")
TPD_3d_plot(PREDICTS_tpds, sites = primary_for,  T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body", method = "prob")
TPD_Diff_plot(data = PREDICTS_tpds, sites1 = primary_for, sites2 = primary_non_for, T1lab = "Locomotion", T2lab = "Foraging",T3lab = "Body",
              method = "prob")
site_dissim <- Calc_dissim(PREDICTS_tpds, sites1 = primary_for, sites2 = primary_non_for)
View(site_dissim)
test <- FD_metrics(PREDICTS_tpds[c("data",sites1,sites2)])
test <- FD_metrics(PREDICTS_tpds[c("data",primary_for,primary_non_for)])
View(PREDICTS_tpds)
