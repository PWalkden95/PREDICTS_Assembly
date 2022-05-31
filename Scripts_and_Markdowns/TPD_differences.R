####################################################################################
####################################################################################
################# TPD Visulation script in addition to mapping dietary niche onto LUR hypervolumes


rm(list = ls())

## load in tidyverse

require(tidyverse)


### function to allow the opening of files in markdowns and non-markdowns

markdown_rds_open <- function(path) {
  open_file <- try(readRDS(path))
  
  if (class(open_file) == "try-error") {
    open_file <- readRDS(paste("../", path, sep = ""))
  }
  
  return(open_file)
  
}

#### realm sites


realms_land_use_sites <- function(realm, land_use) {
  return(
    TPD_LU %>% dplyr::filter(Predominant_habitat == land_use, Realm == realm) %>% dplyr::distinct(SSBS) %>% pull()
  )
}

#### realm land use loop



realm_land_use_loop <-
  function(FUN,
           objects = NULL,
           object_names = NULL) {
    ## if you need some objects to store date in beforehand
    
    
    if(!is.null(objects)){
    object_list <- rep(list(c()),length(objects))
    
    for( i in 1:length(objects)){
    
      if (objects[i] %in% c("vector", "dataframe")) {
        if(i == 1){
          next()
        }
        object_list[[i]] <- c()
      }
      
      if (objects[i] == "list") {
        object_list[[i]] <- list()
      }
      
    }
      
      names(object_list) <- object_names
    }
    
    
    
    
    
    ## for realms and land uses
    
    for (r in realms) {
      for (LU in realm_land_uses[[r]]) {
        value <- FUN(realm = r, land_use = LU)
        
        if (!is.null(objects)) {
          for (i in 1:length(objects)) {
            if (objects[i] == "list") {
              object_list[[i]][[r]][[LU]] <- value[[i]]
            }
            
            if (objects[i] == "vector") {
              object_list[[i]] <- c(object_list[[i]], value)
            }
            
            if (objects[i] == "dataframe") {
              object_list[[i]] <- rbind(object_list[[i]], value)
              
            }
          }
        }
      }
    }
    ## if you have objetcs you want to store data in the store in the object list
    
    
    if (!is.null(objects)) {
      if(length(objects) == 1 & objects %in% c("vector","dataframe")){
        return(object_list[[1]])
      } else {
      
      return(object_list)
      }
    } else {
      return(value)
    }
  }

##read in all relevant files

## obeserved hypervolumes

PREDICTS_tpds <-
  markdown_rds_open("Outputs/PREDICTS_sites_tpds.rds")

## null expectation hypervolumes

PREDICTS_randomisations <-
  readRDS("Outputs/randomisations_TPD_morpho.rds")

## predicts + collpasing secondary and primary vegetation

PREDICTS_full <- readRDS("Outputs/refined_predicts.rds")

PREDICTS <- PREDICTS_full %>%  ## PREDICTS data
  dplyr::distinct(SSBS, Predominant_habitat, Use_intensity, Realm, SS) %>% ## pull out land_use type, Subregion, realm etc
  dplyr::mutate(
    Predominant_habitat = ifelse(
      grepl(
        Predominant_habitat,
        pattern = "secondary",
        ignore.case = TRUE
      ),
      "Secondary vegetation",
      paste(Predominant_habitat)
    ),
    Predominant_habitat = ifelse(
      grepl(
        Predominant_habitat,
        pattern = "primary",
        ignore.case = TRUE
      ),
      "Primary vegetation",
      paste(Predominant_habitat)
    )
  ) %>%
  data.frame()

### individula species tpds

species_TPD <- readRDS("Outputs/species_tpds_morpho.rds")

### regional species pool

species_pools <-
  readRDS("Outputs/predicts_sites_species_pools.rds")

### Dietary niche data for bird species

Forage <-
  readRDS(file = "../PREDICTS_Taxonomy/PREDICTS_imputed_BL_traits_forage.rds")



### data frame of each site and it's land use classification

TPD_LU <-
  data.frame(SSBS = names(PREDICTS_tpds)) %>% dplyr::left_join(PREDICTS) %>% dplyr::filter(Predominant_habitat != "Cannot decide")


## look at the distribution of sites among land uses

table(TPD_LU$Predominant_habitat, TPD_LU$Realm)


## list them for ease of access later

realm_land_uses <- list()
realm_land_uses[["Afrotropic"]] <-
  c(
    "Primary vegetation",
    "Secondary vegetation",
    "Plantation forest",
    "Cropland",
    "Pasture",
    "Urban"
  )
realm_land_uses[["Australasia"]] <-
  c("Primary vegetation",
    "Secondary vegetation",
    "Pasture",
    "Urban")
realm_land_uses[["Indo-Malay"]] <-
  c("Primary vegetation",
    "Secondary vegetation",
    "Plantation forest",
    "Cropland")
realm_land_uses[["Nearctic"]] <-
  c("Primary vegetation", "Pasture", "Cropland", "Urban")
realm_land_uses[["Neotropic"]] <-
  c(
    "Primary vegetation",
    "Secondary vegetation",
    "Plantation forest",
    "Pasture",
    "Cropland",
    "Urban"
  )
realm_land_uses[["Palearctic"]] <-
  c(
    "Primary vegetation",
    "Secondary vegetation",
    "Plantation forest",
    "Cropland",
    "Urban"
  )


## load in TPD analysis functions

source("Functions/TPD_3D_plots.R")


## pull out all the realms

realms <- TPD_LU %>% distinct(Realm) %>% pull() %>% as.character()


################ create some directories to store data visualisations

dir.create("Outputs/TPD_3D_Plots")

for (r in realms) {
  dir.create(paste("Outputs/TPD_3D_Plots", r, sep = "/"))
  dir.create(paste("Outputs/TPD_3D_Plots", r, "land_uses" , sep = "/"))
  dir.create(paste("Outputs/TPD_3D_Plots", r, "land_use_differences" , sep = "/"))
  dir.create(paste("Outputs/TPD_3D_Plots", r, "randomised_holes" , sep = "/"))
  dir.create(paste("Outputs/TPD_3D_Plots", r, "foraging_guilds_mapped" , sep = "/"))
  dir.create(paste(
    "Outputs/TPD_3D_Plots",
    r,
    "foraging_guilds_mapped" ,
    "null",
    sep = "/"
  ))
  dir.create(paste(
    "Outputs/TPD_3D_Plots",
    r,
    "foraging_guilds_mapped" ,
    "observed",
    sep = "/"
  ))
}


## before the plots are produced I need to generate the scale for the legend which is not going to be linear but dependent on the quantiles of
## probability of occupancy of all the hypervolumes from the observed LUR systems




cell_probability_of_occupancy_function <- function(realm,land_use){
  
  land_use_sites <-
    TPD_LU %>% dplyr::filter(Predominant_habitat == land_use, Realm == realm) %>% dplyr::distinct(SSBS) %>% pull()

  
  sites_TPD_data <-
    TPD_plot_data(data = PREDICTS_tpds, site = land_use_sites)
  
  all_cell_probabilities <- sites_TPD_data[["pl_dat"]][["prob"]]
  
  
    return(all_cell_probabilities[which(all_cell_probabilities > 0)])
  

    
}


cell_probability_of_occupancy <- realm_land_use_loop(FUN = cell_probability_of_occupancy_function, objects = c("vector"),
                                                     object_names = c("cell_probability_of_occupancy"))


### with all the probability of occupancy values where the spread might be quite large but the vast amjority of values are clumped quite low a
### linear scale for the legend wouldn't work so I will split up the data in 5000 quantiles

# colours for TPD_plots

legend_col <- as.numeric(quantile(cell_probability_of_occupancy,seq(0.0002,1,0.0002)))


## save these quantiles as they will be used when plotting the TPD but also when generating the legend for figures

write_rds(legend_col, file = "Functions/TPD_colours.rds")


##### now with that done we can get to plotting and saving the TPDs

## create a function for our realm land use loop

plotting_TPD <- function(realm, land_use) {
  sites <- realms_land_use_sites(realm = realm, land_use =  land_use)
  
  
  TPD_3d_plot(
    data = PREDICTS_tpds,
    sites = sites,
    ylab = "Locomotion",
    xlab = "Foraging",
    zlab = "Body",
    method = "prob",
    save = TRUE,
    file = paste(
      "Outputs/TPD_3D_Plots/",
      realm,
      "/land_uses/",
      land_use,
      "_TPD_plot.png",
      sep = ""
    ),
    title = "",
    grid = FALSE,
    free_limits = FALSE
  )
}


realm_land_use_loop(plotting_TPD)

### now we have saved the the observed hypervolumes we want to be able to map dietary guilds onto the functional trait spaces


mapping_dietary_guild_function <- function(realm, land_use) {
  sites <- realms_land_use_sites(realm = realm, land_use = land_use)
  
  tpd_for_dat <-
    TPD_forage_mapping_data(
      data = PREDICTS_tpds,
      randata = PREDICTS_randomisations,
      for_data = Forage,
      sites = sites
    )
  
  return(tpd_for_dat)
  
}

### run the loop

mapping_lists <-
  realm_land_use_loop(
    FUN = mapping_dietary_guild_function,
    objects = c("list", "list"),
    object_names = c("TPD_for_mapping_data",
                     "TPD_for_mapping_data_random")
  )

### extract the lists

TPD_for_mapping_data <- mapping_lists[[1]]
TPD_for_mapping_data_random <- mapping_lists[[2]]

write_rds(file = "Outputs/TPD_forage_mapping.rds", TPD_for_mapping_data)
write_rds(file = "Outputs/TPD_forage_mapping_random.rds", TPD_for_mapping_data_random)


TPD_for_mapping_data <-
  markdown_rds_open("Outputs/TPD_forage_mapping.rds")
TPD_for_mapping_data_random <-
  markdown_rds_open("Outputs/TPD_forage_mapping_random.rds")


###########################  have the observed LUR hypervolumes with idetary nichees mapped visualised


observed_foraging_plot <- function(realm, land_use) {
  tpd_for_dat <- TPD_for_mapping_data[[realm]][[land_use]]
  
  
  TPD_forage_mapping_plot(
    data = PREDICTS_tpds,
    tpd_for_dat,
    T1lab = "",
    T2lab = "",
    T3lab = "",
    title = "",
    save = TRUE,
    s_file = paste(
      "Outputs/TPD_3D_Plots/",
      r,
      "/foraging_guilds_mapped/",
      "observed/",
      LU,
      "_TPD_plot.png",
      sep = ""
    ),
    animation = TRUE,
    a_file = paste(
      "Outputs/TPD_3D_Plots/",
      r,
      "/foraging_guilds_mapped/",
      "observed/",
      LU,
      "_TPD_animation",
      sep = ""
    )
  )
}


### then for the null LURs


null_foraging_plot <- function(realm, land_use) {
  tpd_for_dat_random <-
    TPD_for_mapping_data_random[[realm]][[land_use]]
  
  
  TPD_forage_mapping_plot(
    data = PREDICTS_tpds,
    tpd_for_dat_random,
    T1lab = "",
    T2lab = "",
    T3lab = "",
    title = "",
    save = TRUE,
    s_file = paste(
      "Outputs/TPD_3D_Plots/",
      r,
      "/foraging_guilds_mapped/",
      "observed/",
      LU,
      "_TPD_plot.png",
      sep = ""
    ),
    animation = TRUE,
    a_file = paste(
      "Outputs/TPD_3D_Plots/",
      r,
      "/foraging_guilds_mapped/",
      "observed/",
      LU,
      "_TPD_animation",
      sep = ""
    )
  )
}

realm_land_use_loop(FUN = observed_foraging_plot)
realm_land_use_loop(FUN = null_foraging_plot)






#####################################################
####################################################

proportional_occupancy_function <- function(realm, land_use) {
  sites <- realms_land_use_sites(realm = realm, land_use = land_use)
  
  for_mapping_data <- TPD_for_mapping_data[[realm]][[land_use]]
  for_mapping_data_random <-
    TPD_for_mapping_data_random[[realm]][[land_use]]
  
  
  proportional_occupancy_df <- c()
  
  for (g in c("In", "Gr", "Fr", "Ne", "Om")) {
    prop_occ <-
      proportional_occupancy(
        data = PREDICTS_tpds,
        randata = PREDICTS_randomisations,
        fordata = for_mapping_data,
        ranfordata = for_mapping_data_random,
        sites = sites,
        guild = g
      )
    
    
    df <- data.frame(
      Realm = realm,
      Land_use = land_use,
      Guild = g,
      proportional_occupancy = prop_occ$proportional_occupancy,
      Similarity = prop_occ$similarity
    )
    
    proportional_occupancy_df <-
      rbind(proportional_occupancy_df, df)
    
  }
  
  return(proportional_occupancy_df)
  
}


prop_occ_list <-
  realm_land_use_loop(
    FUN = proportional_occupancy_function,
    objects = c("dataframe"),
    object_names = c("proportional_occupancy_dataframe")
  )




#############################################
#############################################
## what proportion of cells are correctly allocated.

calculate_assignment_accuracy <- function(row) {
  assigned_trophic_niche <- row[, "Trophic_niche"]
  
  if (assigned_trophic_niche == "Unclassified") {
    return(0)
  }
  
  column_names <- colnames(assigned_cells)[6:15]
  
  value <-
    sum(row[, assigned_trophic_niche], na.rm = TRUE) / sum(row[, column_names], na.rm = TRUE)
  
  return(ifelse(is.infinite(value), 1, value))
}

#################################################
#################################################

assignment_accuracy_function <- function(realm, land_use) {
  assigned_cells <- TPD_for_mapping_data[[realm]][[land_use]]
  
  proportion_correct <- c()
  for (i in 1:nrow(assigned_cells)) {
    proportion_correct <-
      c(proportion_correct ,
        calculate_assignment_accuracy(assigned_cells[i,]))
    
  }
  
  proportion_correct_frame <-
    data.frame(
      Realm = realm,
      Land_use = land_use,
      dietary_niche_assignment_accuracy = mean(proportion_correct),
      unassigned_cells = sum(assigned_cells[, "Trophic_niche"] == "Unclassified") /
        nrow(assigned_cells)
    )
  
  return(proportion_correct_frame)
  
  
}


dietary_niche_accuracy <-
  realm_land_use_loop(FUN = assignment_accuracy_function,
                      objects = "dataframe",
                      "dietary_niche_assignment_accuracy")

dietary_niche_accuracy <- dietary_niche_accuracy[[1]]

mean(dietary_niche_accuracy$dietary_niche_assignment_accuracy)
sd(dietary_niche_accuracy$dietary_niche_assignment_accuracy)

mean(dietary_niche_accuracy$unassigned_cells)
sd(dietary_niche_accuracy$unassigned_cells)


########################### couple of extras



###############################################
###############################################


lu_combos <-
  matrix(gtools::combinations(
    n = length(land_use),
    r = 2,
    v = land_use,
    set = TRUE
  ),
  ncol = 2) %>%
  data.frame() %>% set_colnames(c("LU1", "LU2"))


for (r in realms) {
  for (i in 1:nrow(lu_combos)) {
    sites_1 <-
      TPD_LU %>% dplyr::filter(Realm == r, Predominant_habitat == lu_combos[i, "LU1"]) %>% dplyr::distinct(SSBS) %>% pull()
    sites_2 <-
      TPD_LU %>% dplyr::filter(Realm == r, Predominant_habitat == lu_combos[i, "LU2"]) %>% dplyr::distinct(SSBS) %>% pull()
    
    
    if (length(sites_1) > 0 & length(sites_2) > 0) {
      TPD_Diff_plot(
        data = PREDICTS_tpds,
        sites1 = sites_1,
        sites2 = sites_2,
        T1lab = "Locomotion",
        T2lab = "Foraging",
        T3lab = "Body",
        method = "prob",
        title = paste(r, "Difference", lu_combos[i, "LU1"], lu_combos[i, "LU2"], sep = ""),
        save = TRUE,
        file = paste(
          "Outputs/TPD_3D_Plots/",
          r,
          "/land_use_differences/",
          lu_combos[i, "LU1"],
          "_",
          lu_combos[i, "LU2"],
          "_difference_TPD_plot.png",
          sep = ""
        )
      )
      
    }
  }
}


######################################
###################################
## calculate dissimilarity between sites and the randomised communitites






for (r in realms) {
  for (LU in land_use) {
    ransites <-
      TPD_LU %>% dplyr::filter(Predominant_habitat == LU, Realm == r) %>% dplyr::distinct(SSBS) %>% pull()
    
    if (length(ransites) > 0) {
      TPD_ranDiff_plot(
        data = PREDICTS_tpds,
        randata = PREDICTS_randomisations,
        sites = ransites,
        threshold = 1,
        T1lab = "locomotion",
        T2lab = "foraging",
        T3lab = "body",
        title = "",
        method = "prob",
        save = TRUE,
        file = paste(
          "Outputs/TPD_3D_Plots/",
          r,
          "/randomised_holes/",
          LU,
          "_random_holes_TPD_plot.png",
          sep = ""
        ),
        observed = FALSE,
        grid = FALSE,
        free_limits = TRUE
      )
      
      
    }
  }
  
  
}
