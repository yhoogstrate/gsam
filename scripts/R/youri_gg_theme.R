#!/usr/bin/env R

youri_gg_theme <- theme(
  text = element_text(family = 'Helvetica'),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = 'bottom',
  plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
  panel.background = element_rect(fill = 'white', colour = 'white'),
  axis.title = element_text(face = "bold",size = rel(1)),
  axis.title.y = element_text(angle=90,vjust =2),
  axis.title.x = element_text(vjust = -0.2),
  axis.text = element_text(),
  axis.line = element_line(colour="black"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_line(colour = 'grey20', linetype = 'dotted'),
  panel.grid.minor.y = element_line(colour = 'grey50', linetype = 'dotted')
) + theme(
  panel.grid.major.x = element_line(colour = 'grey20', linetype = 'dotted',size=0.25),
  panel.grid.minor.x = element_line(colour = 'grey50', linetype = 'dotted')
)
