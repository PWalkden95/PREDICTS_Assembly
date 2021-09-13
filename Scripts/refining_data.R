rm(list = ls())

require(tidyverse)
require(magrittr)
require(gtools)

######### Load in PREDICTS dtaaset with Birdlife taxonomy and traits

PREDICTS_Aves <- readRDS("../PREDICTS_Taxonomy/PREDICTS_BL_traits.rds")

### For the study we would  like to have abundance recorded sites with greater than 1 species & is focussed on birds or surveys the
### entire community to ensure that survey is capturing the bird community to the fullest extent - those focussed on specific species
### will not be suitable for analyses. 

keep <- c("FE1_2015__Marinaro","JB3_2013__Gove", "DL1_2011__Latta", "DB1_2010__Dures")

PRED_Assem <- PREDICTS_Aves %>% dplyr:: filter(Diversity_metric_type == "Abundance", Measurement > 0) %>% dplyr::group_by(SS) %>%
  dplyr::mutate(spp_num = n_distinct(Birdlife_Name)) %>% dplyr::filter(spp_num > 1, grepl(Study_name, pattern = "bird", ignore.case = TRUE)|Sampling_target == "Entire community"|Source_ID %in% keep) %>%
  droplevels()



## after manually checking that the study surveyed the complete community of birds (Some studies marked as surveying Entire communities)
## had occassionally targeted specfic species, unsuitable for commuity assembly analyses. Espc. when working with observational data

drop <- data.frame(study = c("VK1_2011__Zimmerman", "DI1_2004__Naidoo", "JB3_2016__Gardner",
                             "HP1_2010__Bicknell" , "SE2_2010__McCarthy", "SE1_2011__Arbelaez",
                             "HB1_2009__Parry", "VK1_2007__StLaurent", "DI1_2011__Dawson", "DL1_2013__Bartolommei", "LK1_2009__Hayward"),
                   reason= c("songbirds", "songbirds", "group_abundance",
                             "targeted sampling", "camera_traps", "Mixed Bird Flocks", "not entire community",
                             "passerines", "proportion of plots", "Owls and nocturnal birds", "Not whole bird community"))

## filter out these studies 

refine_data <- PRED_Assem %>% filter(!(Source_ID %in% drop$study)) %>% filter(SS != "GN1_2010__Hvenegaard 2")

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


#### enter PREDICTS data and Source_ID of studies to be combined 
source("Functions/Combine_studies.R")


studies <- as.character(unique(refine_data$SS))
reduce <- substr(studies,1,(nchar(studies)-2))
duplicate <- unique(reduce[which(duplicated(reduce))])

retain <- c("MJ1_2016__Jung", "XH1_2015__Rachakonda")  ### each of the studies are distinct sites 
duplicate <- duplicate[!(duplicate %in% retain)]  



for(dup in duplicate){
refine_data <- combine_study(refine_data, dup) %>% droplevels()
}

###################################################
##################################################
# Here there are still many sites whose use intensity is classified as cannot decide therefore after some double checking we can assign classification
refine_data <- refine_data %>% dplyr::muate(Use_intensity = ifelse(SS == "DL1_2009__Azpiroz 1" & Predominant_habitat == "Cropland",
                                                                   "Minimal use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "DL1_2010__Proenca 2" & Predominant_habitat == "Plantation forest",
                                                                   "Minimal use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "DL1_2011__Mallari 1",
                                                                   "Minimal use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "DL1_2012__Dallimer 1" & Predominant_habitat == "Secondary vegetation (indeterminate age)",
                                                                   "Minimal use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SSBS == c("TN1_2006__Soh 1 Cameron 1","TN1_2006__Soh 1 Fraser 6") ,
                                                                   "Minimal use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "AD1_2012__Yamaura 1" & Predominant_habitat == "Pasture",
                                                                   "Light use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "AD1_2012__Yamaura 1" & Predominant_habitat == "Young Secondary Vegetation",
                                                                   "Minimal use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "AD1_2012__Yamaura 1" & Predominant_habitat == "Plantation forest",
                                                                   "Minimal use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "GN1_2010__Hvenegaard 1" & Predominant_habitat == "Cropland",
                                                                   "Light use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SSBS %in% c("HP1_2008__Gomes 1  1","HP1_2008__Gomes 1  9"),
                                                                   "Minimal use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SSBS %in% c("HP1_2008__Gomes 1  5","HP1_2008__Gomes 1  6"),
                                                                   "Light use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "HP1_2008__Gomes 1" & Predominant_habitat == "Secondary vegetation (indeterminate age)",
                                                                   "Light use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "SH1_2014__Walker 2",
                                                                   "Minimal use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "GN1_2018__AlvarezAlvarez 1" & Predominant_habitat == "Cropland",
                                                                   "Minimal use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "GN1_2018__Karp 1" & Predominant_habitat == "Pasture",
                                                                   "Light use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "HW1_2007__Chapman 1" & Use_intensity == "Cannot decide" & Predominant_habitat == "Urban",
                                                                   "Minimal use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "HW1_2007__Chapman 1" & Use_intensity == "Cannot decide",
                                                                   "Light use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "HW1_2011__Cerezo 1" & Predominant_habitat == "Cropland",
                                                                   "Light use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "HW1_2011__Cerezo 1" & Predominant_habitat == "Pasture",
                                                                   "Intense use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "MJ1_2008__Munyekenye 1" & Predominant_habitat == "Primary forest",
                                                                   "Minimal use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "MJ1_2008__Munyekenye 1" & Predominant_habitat == "Secondary vegetation (indeterminate age)",
                                                                   "Minimal use",paste(Use_intensity)),
                                            Use_intensity = ifelse(SS == "MJ1_2013__Reynolds 1" & Predominant_habitat == "Mature secondary vegetation",
                                                                   "Intense use",paste(Use_intensity)))




##### Additionally some sites that are geographically proximate or are point counts within a continuous habitat need to be merged to faciliate easier
##### quantification of functional diversity metrics as a fuller account of the community is derived

source("Functions/refine_data_blocks.R")

### first load in some manual inspection of studies and sites determine whether and how merging should take place

merge_studs <- read.csv("Outputs/checked_studies.csv")

source("Functions/Site_merging.R")

refine_data <- site_merge(refine_data,merge_studs)




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
  scale_x_continuous(limits = c(-180,180),breaks = seq(-180,180,30))+
  scale_y_continuous(limits = c(-90,90), breaks = seq(-90,90,30))+
  theme_classic()

plot(site_plot)


## Save 

write_rds(refine_data, "Outputs/refined_predicts.rds")