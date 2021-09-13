rm(list = ls())

require(tidyverse)


randomisations <- readRDS("Outputs/site_randomisations.rds")
dist_traits <- readRDS("Outputs/study_species_pool_dist_ninety.rds")
PREDICTS <- readRDS("Outputs/refined_predicts.rds")


source("Functions/PREDICTS_rao_func.R")


random_FD_calc <- function(list){
  
  site_name <- substr(colnames(list)[3],1,nchar(colnames(list)[3])-9)
  
  study <- PREDICTS %>% dplyr::filter(SSBS == site_name) %>% distinct(SS) %>% pull()
  
  dist_mat <- dist_traits[[study]]
  
  
  t_frame <- data.frame(RelativeAbundance = rep(list[["RelativeAbundance"]],1001),Birdlife_Name = NA, SSBS = NA)
  
  start_row <- 1
  end_row <- nrow(list)
  for(i in 2:ncol(list)){
    if(start_row ==1){
      t_frame[c(start_row:end_row), "SSBS"] <- site_name
      t_frame[c(start_row:end_row), "Birdlife_Name"] <- list[,i]
    } else {
      t_frame[c(start_row:end_row), "SSBS"] <- colnames(list)[i]
      t_frame[c(start_row:end_row), "Birdlife_Name"] <- list[,i]
    }
    
    start_row <- end_row +1
    end_row <- end_row + nrow(list)
  }
 
  comm_data <- t_frame %>% distinct(SSBS)
  comm_rao <- Rao_Q_Func_bias(data = t_frame,traits = dist_mat)
  
  ses <- (comm_rao$FunRao[1] - mean(comm_rao$FunRao[2:length(comm_rao$FunRao)]))/sd(comm_rao$FunRao[2:length(comm_rao$FunRao)])
  
  return(ses)
}

ses <- lapply(randomisations,random_FD_calc)

write_rds(ses, file = "Outputs/randomisations_raos_q.rds")
ses <- readRDS("Outputs/randomisations_raos_q.rds")

PREDICTS <- readRDS("Outputs/refined_predicts.rds") %>% dplyr::select(SSBS,SS,Predominant_habitat,Use_intensity, Biome) %>%
  distinct(SSBS,SS,Predominant_habitat,Use_intensity, .keep_all = TRUE)
drop_sites <- PREDICTS %>% dplyr::filter(Predominant_habitat == "Cannot decide") %>% pull(SSBS) %>% as.character()

ses <- ses[-which(names(ses)%in%drop_sites)]



ses_scores <- data.frame(SSBS = names(randomisations)) %>% dplyr::left_join(PREDICTS) %>% filter(Predominant_habitat != "Cannot decide")
ses_scores$ses <- unlist(ses) 

ses_scores <- ses_scores %>% dplyr::mutate(Predominant_habitat = ifelse(grepl(Predominant_habitat, pattern = "Primary"),
                                                                        "Primary", paste(Predominant_habitat)),
                                           Predominant_habitat = ifelse(grepl(Predominant_habitat, pattern = "secondary",ignore.case = TRUE),
                                                                        "Secondary vegetation", paste(Predominant_habitat)),
                                           LUI = paste(Predominant_habitat,Use_intensity, sep = "_"),
                                           LUI = factor(LUI),
                                           LUI = relevel(LUI, ref = "Primary_Minimal use")) 

sp_rich <- c()
for(sit in names(ses)){
  sp_rich <- c(sp_rich,nrow(randomisations[[sit]]))
}

ses_scores$species_richness <- sp_rich



require(lme4)


table(ses_scores$Predominant_habitat)

test <- lmer(ses ~ LUI + Biome + (1|SS), data = ses_scores) 

summary(test)
