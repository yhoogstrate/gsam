#!/usr/bin/env R


setwd("~/projects/gsam")


library(ggplot2)
library(tidyverse)
#library(ggforce)

rm(geomnet)
rm(geom_circle)
detach("package:geomnet")
detach("package:geomnet", unload=TRUE)
detach("function:geom_circle", unload=TRUE)
remove.packages('geomnet')
devtools::install_github("yhoogstrate/geomnet")
library(geomnet)



col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))





cor_cor_plot <- function(normalised_correlated_data, labels) {
  
  normalised_correlated_data <- plt
  
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
  h <- hclust( as.dist(1 - cor(plt)) ) # Geniale manier om te clusteren!!!
  o <- h$labels[h$order]
  rm(h)
  
  
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
          axis.text.y = element_text(size = 9, angle = 0, hjust = 1, vjust = 0.5),
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
  
  
  p2 <- ggplot(plt, aes(x = i , y = variable, fill=value, label=gid)) +
    geom_tile(col='white',lwd=0.6) +
    #scale_x_discrete(position = "bottom")  + 
    scale_x_discrete(labels = NULL, breaks = NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank()
          ) +
    guides(fill="none") +
    ggplot2::coord_fixed() +
    labs(x=NULL, y=NULL) + 
    scale_fill_manual(values=c('TRUE'='gray40','FALSE'='gray95')) 
    #xlim(1, length(o) )

  
  p2 / p1
  
  
}


cor_cor_plot(plt, labels)




# plt.expanded2 <- plt.expanded2 %>%
#   dplyr::mutate(r = ((abs(value) * 0.7) + 0.3) / 2 - 0.05) %>%
#   dplyr::mutate(r = r / 15)



#saveRDS(plt.expanded2 , file="tmp/tmp.Rds")
plt.expanded2  <- readRDS(file="tmp/tmp.Rds")

# , 
# fill=value,
# col=value,
# label=x

plt <- plt.expanded2 %>% dplyr::slice_head(n=2)

ggplot(plt, aes( x = x.order, y = y.order, radius = r)) +
  geom_circle(radius.fixed = T) # not to be confused w/ geomnet::geom_circle




