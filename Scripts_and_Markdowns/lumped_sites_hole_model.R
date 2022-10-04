require(future)
require(future.apply)
require(car)
require(lmerTest)
require(merTools)



  
  
  
  
  pull_sites <- function(study){
  
  sites <- TPD_LU %>% filter(SS == study)
  
  land_uses <- sites %>% dplyr::distinct(Predominant_habitat) %>% pull()
  
  site_list <- list()
  
  for(LU in land_uses){
    sits <- sites %>% dplyr::filter(Predominant_habitat == LU) %>% dplyr::distinct(SSBS, .keep_all = FALSE) %>% pull() 
    
    site_list[[paste(study,LU,sep = "/")]] <- sits
    }
    return(site_list)  
  }
  
  
  studies <- unique(TPD_LU$SS)
  
  lumped_sites <- list()
  
  for(study in studies){
    
    stud_sites <- pull_sites(study = study)
    
    lumped_sites <- c(lumped_sites,stud_sites)
  }
  
  
  site_hole <- c()
  for(i in 1:length(lumped_sites)){
  
    
    
    data <- data.frame(SS = str_split(names(lumped_sites)[i],pattern = "/")[[1]][1],
                       Predominant_habitat = str_split(names(lumped_sites)[i],pattern = "/")[[1]][2])
    
    hole <- TPD_holes(data = PREDICTS_tpds, sites = lumped_sites[[i]],randata = PREDICTS_randomisations, threshold = 0.95)
      
    
    if(is_empty(hole)){
      data$internal_hole_vol <- NA
      data$internal_hole_prop <- NA
      
      data$external_hole_vol <- NA
      data$external_hole_prop <- NA
      
    }else{
      
      data$internal_hole_vol <- hole$internal$total_hole_volume
      data$internal_hole_prop <- hole$internal$proportion_holes_volume
      
      data$external_hole_vol <- hole$external$total_hole_volume
      data$external_hole_prop <- hole$external$proportion_holes_volume
      
    }
    
    fd_mets <- TPD_FD_metrics(data = PREDICTS_tpds, sites = lumped_sites[[i]])
    
    data[,c("FRich","FDiv","FEve","FRound")] <- as.numeric(fd_mets[nrow(fd_mets),])
    
    site_hole <- rbind(site_hole,data)
    
  }

  
  
site_hole <- site_hole %>% dplyr::left_join(PREDICTS[,c("SS","Realm")], by = "SS")  %>% dplyr::distinct(SS,Predominant_habitat, .keep_all = TRUE)
  


site_hole$Predominant_habitat <- relevel(factor(site_hole$Predominant_habitat), ref = "Primary forest")
site_hole$Realm <- relevel(factor(site_hole$Realm), ref = "Neotropic")




table(site_hole$Realm, site_hole$Predominant_habitat)




site_hole$logFRich <- log(site_hole$FRich)



site_mod <- lmer(sqrt(internal_hole_prop) ~ Predominant_habitat + logFRich + Realm + Predominant_habitat:Realm + (1|SS), data = site_hole)



summary(site_mod)
Anova(site_mod, type = "III")

ggResidpanel::resid_panel(site_mod)




realms <- c(rep("Afrotropic",7),
            rep("Australasia",7),
            rep("Nearctic",7),
            rep("Neotropic",7),
            rep("Indo-Malay",7),
            rep("Palearctic",7))

LUs <- rep(c("Primary forest", "Primary non-forest", "Secondary vegetation", "Plantation forest", "Pasture", "Cropland", "Urban"),6)



lala <- data.frame(Realm = realms, Predominant_habitat = LUs) 
lala$logFRich <- log(mean(site_hole$FRich))


#for(i in 1:nrow(lala)){
#  ph <- lala[i,"Predominant_habitat"]
#  r <- lala[i,"Realm"]
#lala[i,"logF_Rich"] <- log(mean(site_hole[which(site_hole$Predominant_habitat == ph & site_hole$Realm == r),"F_Rich"]))
#}

drop_rows <- c(30,33, 17)

lala <- lala[-drop_rows,]


pred_fun <- function(x){
  as.numeric(predict(x, newdata = lala, re.form = NA))
}

lala$predicted <- pred_fun(site_mod)^2
lala$upper <- NA
lala$lower <- NA


booted_mod <- bootMer(site_mod, FUN = function(x) pred_fun(x), nsim = 200)

estimates <- booted_mod$t0^2

bootstraps <- booted_mod$t^2

for(i in 1:ncol(bootstraps)){
  lala[i,"upper"] <- as.numeric(quantile(bootstraps[,i], 0.95))
  lala[i,"lower"] <- as.numeric(quantile(bootstraps[,i], 0.05))
  
}

###############################
###############################
# Now the plot

lala$Predominant_habitat <- factor(lala$Predominant_habitat, levels = c("Primary forest", "Primary non-forest", "Secondary vegetation", "Plantation forest", "Pasture", "Cropland", "Urban"))

levels(lala$Predominant_habitat) <- c("PriFor","PriNFor","SecVeg","PlnFor","Pas","Crp","Urb")


for(r in unique(realms)){
  
  p_d <- lala %>% dplyr::filter(Realm == r)
  
  pri_intercept <- p_d %>% dplyr::filter(Predominant_habitat == "PriFor") %>% pull(predicted)
  
  
  pd <- position_dodge(0.5)
  plot <- ggplot(data = p_d, aes(x = Predominant_habitat, y = predicted, colour = Predominant_habitat)) +
    geom_errorbar(aes(ymin=ifelse(lower< 0, 0, lower) , ymax= upper), colour="black", width=.1, position=pd, linetype = 1) +
    geom_point(position=pd,size=6)+
    xlab("Land use") +
    ylab("Holeyness") +
    scale_colour_hue(name="Use intensity",    # Legend label, use darker colors
                     breaks=c("Minimal", "Light","Intense","All"),
                     labels=c("Minimal", "Light","Intense","All"),
                     l=40) +                    # Use darker colors, lightness=40+
    expand_limits(y=0) +
    ggtitle(paste(r))+
    theme_classic() +
    geom_hline(yintercept= pri_intercept, linetype='dotted', col = 'red')+
    theme(legend.justification=c(1,0),
          legend.position=c(1,0.65),
          text = element_text(size = 20))+
    theme(axis.text.x = element_text(angle = 90))
  
  plot(plot)
}
