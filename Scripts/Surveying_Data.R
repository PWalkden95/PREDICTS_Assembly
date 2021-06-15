rm(list = ls())

require(tidyverse)
require(magrittr)
require(gtools)

####### Load in PREDICTS, traits & Taxonomy data

Forage <- readRDS("../../Datasets/GBD/Trophic_Foraging_Niche.rds")
PREDICTS <- readRDS("../../Datasets/PREDICTS/diversity-2021-02-24-03-32-59.rds")
PREDICTS_Taxonomy <- readRDS("../../Datasets/PREDICTS/PREDICTS_AVES_Updated_Taxonomy.rds")
Jetz_Traits <- read.csv("../../Datasets/GBD/GBD_Jetz_averages_11_Nov_2020.csv")

#### Change colname for simplicity

colnames(Jetz_Traits)[1] <- "Jetz_Name" 

#### Filter PREDICTS to birds that are resolved to the species level, join updated taxonomy to be compatable to Jetz phylogeny, join traits

PREDICTS_Aves <- PREDICTS %>% dplyr::filter(Rank == "Species", Class == "Aves") %>% 
  dplyr::left_join(PREDICTS_Taxonomy[,c("Taxon", "Jetz_Name")], by = "Taxon") %>%
  dplyr::left_join(Forage[,c(34:46)], by = c("Jetz_Name" = "Taxon")) %>%
  dplyr::left_join(Jetz_Traits[,c(1,8:18)], by = "Jetz_Name")


### For the study we would  like to have abundance recorded sites with greater than 1 species & is focussed on birds or surveys the
### entire community to ensure that survey is capturing the bird community to the fullest extent - those focussed on specific species
### will not be suitable for analyses. 

PRED_Assem <- PREDICTS_Aves %>% dplyr:: filter(Diversity_metric_type == "Abundance", Measurement > 0) %>% dplyr::group_by(SS) %>%
  dplyr::mutate(spp_num = n_distinct(Jetz_Name)) %>% dplyr::filter(spp_num > 1, grepl(Study_name, pattern = "bird", ignore.case = TRUE)|Sampling_target == "Entire community") %>%
  droplevels()


## after manually checking that the study surveyed the complete community of birds (Some studies marked as surveying Entire communities)
## had occassionally targeted specfic species, unsuitable for commuity assembly analyses. Espc. when working with observational data

drop <- data.frame(study = c("VK1_2011__Zimmerman", "DI1_2004__Naidoo", "JB3_2016__Gardner", "JB3_2017__Kajtoch",
                             "MJ1_2009__Lehouck", "HP1_2010__Bicknell" , "JB3_2019__Leso", "SE2_2010__McCarthy", "SE1_2011__Arbelaez",
                             "HB1_2009__Parry", "DL1_2010__Sosa", "VK1_2007__StLaurent", "DI1_2011__Dawson"),
                   reason= c("songbirds", "songbirds", "group_abundance", "orchard birds",
                             "frugivores","targeted sampling", "targeted sampling", "camera_traps", "Mixed Bird Flocks", "not entire community",
                             "couldnt find paper", "passerines", "proportion of plots"))

## filter out these studies 

refine_data <- PRED_Assem %>% filter(!(Source_ID %in% drop$study))

#### Correct for sampling effort 

refine_data <- refine_data %>% 
  
  
  dplyr::group_by(SS) %>% ## group by study 
  
  dplyr::mutate(Max_Sampling_Effort = max(Sampling_effort)) %>%  # max samplin effort in each site
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(Rescaled_Sampling_Effort = Sampling_effort/Max_Sampling_Effort) %>%  ## rescale sampling effort so that the maximum effort with a study is 1
  
  dplyr::mutate(Effort_Corrected_Measurement = ifelse(Diversity_metric_is_effort_sensitive == TRUE,
                                                      Measurement/Rescaled_Sampling_Effort,
                                                      Measurement)) %>%
  
  droplevels()



## Lasky study is split into surveys conducted in the wet and dry season and they should be combined 


data <- refine_data
stud_combine <- "HP1_2010__Lasky"
#### enter PREDICTS data and Source_ID of studies to be combined 

combine_study <- function(data, stud_combine) {
  
  study <- data %>% filter(grepl(SS, pattern = stud_combine)) %>% droplevels()

cat(sprintf("Number of blocks in study: %i\n", length(unique(study$SS))))

cross <- as.data.frame(unclass(table(study$SS, study$Site_name)))

for(block in rownames(cross)){
  cat(sprintf("Number of sites in block %s: %i\n", block, length(which(cross[block,] > 0))))
}
  
combinations <- gtools::combinations(n = length(rownames(cross)), r = 2, v = rownames(cross))

for(i in 1:nrow(combinations)){
same <- which(c(cross[combinations[1],] > 0) == c(cross[combinations[2],] > 0))
diff <- as.numeric(which(c(cross[combinations[1],] > 0) != c(cross[combinations[2],] > 0)))

cat(sprintf("Block %s surveyed %i of the same sites as block %s\n", combinations[1], length(same), combinations[2]))

if(is_empty(diff)){
  cat("No difference in sites")
} else {
  cat(sprintf("sites different: %i \n", length(diff)))}
}

cat("Combining the same sites....")

combine_stud <- study %>% group_by(Site_name, Jetz_Name) %>%
    dplyr::mutate(Effort_Corrected_Measurement = sum(Effort_Corrected_Measurement)) %>% ungroup() %>%
  dplyr::distinct(Site_name, Jetz_Name, .keep_all = TRUE) 

cat("Resolving SSBS names...")

site_names <- dplyr::distinct(combine_stud, Site_name) %>% pull() %>% as.character()

j <- 1
for(site in site_names){
combine_stud <- combine_stud %>% dplyr::mutate(SSBS = ifelse(Site_name == site, paste(stud_combine," 1 ",j, sep = ""), paste(SSBS)))
j <- j + 1
}

combine_stud <- combine_stud %>% dplyr::mutate(SSBS = factor(SSBS),
                               SS = paste(stud_combine," 1"),
                               SS = factor(SS),
                               SSB = paste(stud_combine," 1 ", "1",sep = ""),
                               SSB = factor(SSB),
                               SSS = SSBS)


combine_stud <- data %>% filter(!(grepl(SS, pattern = stud_combine))) %>% rbind(combine_stud)

return(combine_stud)

}

refine_data <- combine_study(refine_data, "JB3_2019__daSilva")

test <- refine_data %>% filter(grepl(SS, pattern = "daSilva")) %>% droplevels()

table(test$Site_name, test$Study_number)


unique(refine_data$SS)

#### Visulise the sites geographically

wm<-map_data("world") %>% filter(region != "Antartica" ) %>% fortify()

## site coords

site_points <- refine_data %>% distinct(SSBS,Longitude,Latitude)

# generate and plot map

site_plot<-ggplot()+ coord_fixed()+
  geom_map(data =wm, map = wm,
           aes(group = group, map_id= region),
           fill = "darkgrey")+
  geom_point(data = fortify(site_points), aes(Longitude, Latitude),
             colour = "blue", size = 1)+
  theme_classic()

plot(site_plot)


diff<- which(FALSE)


## Save 

write_rds(refine_data, "Outputs/refined_predicts.rds")
  