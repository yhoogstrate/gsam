#!/usr/bin/env R

# make scalebar

# as from code 
col2 <- grDevices::colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                      "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                      "#4393C3", "#2166AC", "#053061"))


plt <- data.frame(c = -20:20/20) |> 
  dplyr::mutate(x = 1:n())

ggplot(plt, aes(x=x,y=c, col=c)) +
  geom_point() +
  ggplot2::scale_fill_gradientn(colours = col2(200), na.value = "grey50", limits = c(-1, 1)) + # guide = "colourbar",
  ggplot2::scale_color_gradientn(colours = col2(200), na.value = "grey50", limits = c(-1, 1)) 

ggsave("output/figures/recursive_cor_plot__w__legend.pdf")
