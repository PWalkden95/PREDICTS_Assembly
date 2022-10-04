

require(raster)
require(tidyverse)


my.colors <-
  colorRampPalette(c("white", "orange", "green", "darkgreen"))


colour_df <- data.frame(value = legend_col,colour = my.colors(length(legend_col)))


sequence <- seq(0,max(legend_col[-5001]), length = 5000)

legend_matrix <- data.frame(legend_axis = sequence,
                            colour_value = legend_col[-5001],
                            colour = my.colors(5000),
                            position = 1)



plot <- ggplot(data = legend_matrix, aes(x = legend_axis, y = position, fill = legend_axis))+
  geom_raster(show.legend = FALSE) +
  geom_vline(xintercept = 0) +
  scale_fill_gradientn(colours=legend_matrix$colour,
                       breaks=legend_matrix$colour_value) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank()) + 
  annotate("text", label = paste(as.character(signif(max(legend_matrix$colour_value), digits = 3))," - ",signif(max(legend_col),digits = 3)), x = max(legend_matrix$legend_axis), y = 0.4, size = 20,fontface = "bold") +
  geom_segment(aes(x=max(legend_matrix$legend_axis), xend=max(legend_matrix$legend_axis), y=0.48, yend=0.55), colour="black",lwd = 2) +
  geom_segment(aes(x=sequence[find_position(max(sequence)*0.1, sequence)], xend=sequence[find_position(max(sequence)*0.1, sequence)], y=0.48, yend=0.55), colour="black",lwd = 2) +
  geom_segment(aes(x=sequence[find_position(max(sequence)*0.2, sequence)], xend=sequence[find_position(max(sequence)*0.2, sequence)], y=0.48, yend=0.55), colour="black",lwd = 2) +
  geom_segment(aes(x=sequence[find_position(max(sequence)*0.3, sequence)], xend=sequence[find_position(max(sequence)*0.3, sequence)], y=0.48, yend=0.55), colour="black",lwd = 2) +
  geom_segment(aes(x=sequence[find_position(max(sequence)*0.4, sequence)], xend=sequence[find_position(max(sequence)*0.4, sequence)], y=0.48, yend=0.55), colour="black",lwd = 2) +
  geom_segment(aes(x=sequence[find_position(max(sequence)*0.5, sequence)], xend=sequence[find_position(max(sequence)*0.5, sequence)], y=0.48, yend=0.55), colour="black",lwd = 2) +
  geom_segment(aes(x=sequence[find_position(max(sequence)*0.6, sequence)], xend=sequence[find_position(max(sequence)*0.6, sequence)], y=0.48, yend=0.55), colour="black",lwd = 2) +
  geom_segment(aes(x=sequence[find_position(max(sequence)*0.7, sequence)], xend=sequence[find_position(max(sequence)*0.7, sequence)], y=0.48, yend=0.55), colour="black",lwd = 2) +
  geom_segment(aes(x=sequence[find_position(max(sequence)*0.8, sequence)], xend=sequence[find_position(max(sequence)*0.8, sequence)], y=0.48, yend=0.55), colour="black",lwd = 2) +
  geom_segment(aes(x=sequence[find_position(max(sequence)*0.9, sequence)], xend=sequence[find_position(max(sequence)*0.9, sequence)], y=0.48, yend=0.55), colour="black",lwd = 2) 




for(i in seq(0.1,0.9,0.1)){
  
  axis_position <- sequence[find_position(max(sequence)*i, sequence)]
  axis_value <- legend_matrix$colour_value[find_position(max(sequence)*i, sequence)]
  
  plot <- plot +
    annotate("text", label = paste(as.character(signif(axis_value, digits = 3))), x = axis_position, y = 0.4, size = 20,fontface = "bold") 
   
  
}

  
  



ggsave("C:/Users/patri/Desktop/test.png",plot,device = "png", height = 8, width = 55,dpi=300, limitsize = FALSE )



raster_value <- c()
for(i in legend_col){
  raster_value <- c(raster_value,find_position(i,legend_col))
  
}


raster_value[4000]

colour_df <- data.frame(legend_axis = sequence, value = raster_value) %>%
  dplyr::left_join(legend_matrix, by = legen)






create_raster <- raster(nrow = 100, ncol = 10000, xmn = 0, xmx = max(legend_col), ymn=0, ymx=0.5)

raster_matrix <- matrix(rep(seq))


res(create_raster) <- c(max(legend_col)/10000,0.5/100)

extent(create_raster) <- c(0,max(legend_col),0,0.5)

test <- create_raster









raster_plot <- ggplot()+
  geom_raster()



test

plot(test, xlim = c(0,max(legend_col)))
abline(v = 0)
axis(1, labels = c(0.2,0.4,0.6))

test




