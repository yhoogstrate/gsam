#!/usr/bin/env R


setwd("~/projects/gsam")


library(ggplot2)
library(tidyverse)

# rm(geomnet)
# rm(geom_circle)
# detach("package:geomnet")
# detach("package:geomnet", unload=TRUE)
# detach("function:geom_circle", unload=TRUE)
# remove.packages('geomnet')
# devtools::install_github("yhoogstrate/geomnet")
# library(geomnet)



#devtools::install_github("yhoogstrate/geomnet")
library(geomnet)


col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))




# > plt
#             AAB1     AAC1     AAD1     AAF1     AAF2     AAJ1
# AJAP1   6.792060 5.161832 7.320888 5.750563 8.084654 6.998988
# CHD5    7.048284 4.729151 5.568896 6.923165 9.552252 8.421078
# FBXO2   5.517842 5.274932 5.099118 4.902737 7.686516 6.200491
# GABRD   5.962070 6.069547 7.433840 5.167848 7.068670 5.738693
# HSPB7   5.158878 5.375160 6.101220 5.670945 5.787876 5.464138
# PLCH2   6.027242 5.120413 4.676221 5.810424 7.591514 6.172717
# PRKCZ   6.027242 6.598770 5.436837 6.295483 7.761095 7.659296
# TMEM88B 4.995062 4.243987 4.550191 4.243987 5.698900 4.621858


# > labels
#         direction_up direction_down
# AJAP1           TRUE          FALSE
# CHD5            TRUE          FALSE
# FBXO2           TRUE          FALSE
# GABRD           TRUE          FALSE
# HSPB7           TRUE          FALSE
# PLCH2           TRUE          FALSE
# PRKCZ           TRUE          FALSE
# TMEM88B         TRUE          FALSE


cor_cor_plot <- function(normalised_correlated_data, labels, method="ward.D2") {
  
  #normalised_correlated_data <- plt
  
  # remove duplicate entries:
  plt <- normalised_correlated_data %>%
    tibble::rownames_to_column('__hugo_symbol__') %>%
    dplyr::filter(!duplicated(`__hugo_symbol__`)) %>%
    tibble::column_to_rownames('__hugo_symbol__')
  
  # determine correlation
  plt <- plt %>%
    as.matrix %>% 
    t() %>%
    cor()
  
  # find order by taking correlation of the correlation
  h <- hclust( as.dist(1 - cor(plt)) , method = method ) # Geniale manier om te clusteren!!!
  o <- h$labels[h$order] %>% rev()
  
  ph <- ggdendro::ggdendrogram(h, rotate = TRUE, theme_dendro = FALSE) +
    ggdendro::theme_dendro()
  
  # re-order to cor-cor clustering order and transform from data.frame into matrix
  plt <- plt %>%
    as.data.frame %>%
    dplyr::select(o) %>%
    t() %>%
    as.data.frame %>%
    dplyr::select(o) %>%
    t() %>%
    as.matrix
  
  
  # to test:
  # corrplot::corrplot(plt)
  

    
  # add x and y co-ordinates later on to melted table
  o.join <- data.frame(name = o, i = 1:length(o))
  
  plt.expanded2 <- reshape2::melt(plt) %>%
    dplyr::rename(y = Var1) %>%
    dplyr::rename(x = Var2) %>%
    dplyr::mutate(x = as.factor(x)) %>%
    dplyr::mutate(y = as.factor(y)) %>%
    dplyr::left_join(o.join %>% dplyr::rename(x.order = i), by=c('x' = 'name'))%>%
    dplyr::left_join(o.join %>% dplyr::mutate(i = nrow(.) - i + 1  ) %>% dplyr::rename(y.order = i), by=c('y' = 'name'))
  
  rm(o.join)
  
  
  
  p1 <- ggplot(plt.expanded2,
         aes( x = x.order, y = y.order, 
             radius = ((abs(value) * 0.7) + 0.3) / 2 - 0.05 ,  # [0.3 , 0.8] + 0.2 smoothened from lwd/border
             fill=value,
             col=value,
             label=x)
         ) +
    geom_tile( col="gray", fill="white", lwd=0.3) +
    scale_fill_gradientn( colours = col2(200), na.value = "grey50", limits = c(-1,1) , guide="none") + # guide = "colourbar",
    scale_color_gradientn( colours = col2(200), na.value = "grey50", limits = c(-1,1) , guide="none" ) +
    geomnet::geom_circle(radius.fixed = T) + # not to be confused w/ geomnet::geom_circle
    scale_x_discrete(labels = NULL, breaks = NULL) +
    theme(legend.position = 'bottom',
          axis.text.y = element_text(size = 6, angle = 0, hjust = 1, vjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, color="gray80"),
          
          text = element_text(size=13),
          
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank()
          ) +
    labs(y = NULL, x=NULL, main=NULL) +
    ggplot2::coord_fixed() +
    scale_y_continuous(name=NULL, breaks = length(o):1, labels = o)
  
  
  
  plt <- data.frame(gid = o, i = 1:length(o)) %>%
    dplyr::left_join(labels %>% tibble::rownames_to_column('gid'), by=c('gid' = 'gid')) %>%
    reshape2::melt(id.vars = c('gid','i'))  
  
  
  p2 <- ggplot(plt , aes(x = i , y = variable , fill=value, label=gid)) +
    geom_tile(col='white',lwd=0.6) +
    #scale_x_discrete(position = "bottom")  + 
    scale_x_discrete(labels = NULL, breaks = NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.text.y = element_text(size = 7, angle = 0, hjust = 1, vjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank()
          ) +
    guides(fill="none") +
    ggplot2::coord_fixed(ratio = 2.75) +
    labs(x=NULL, y=NULL) + 
    scale_fill_manual(values=c('TRUE'='red','FALSE'='gray95')) 
    #scale_fill_manual(values=c('TRUE'='gray40','FALSE'='gray95')) 

  
  #(p2 / p1) / ph
  
  #(p2 + plot_layout(guides = 'collect') ) /
  #(p1 + (ph))
  
  
  layout <- '
A#
BC'
  wrap_plots(A = p2, B = p1, C = (ph + plot_spacer () )  , design = layout)
  
  #  (p2 + plot_spacer()  + plot_layout(widths = c(1, 0.5)) ) /
  #  (p1 + ph + plot_layout(widths = c(1, 0.5)))
    
    
  
  # (p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
  #  (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt")))
  
  
}



