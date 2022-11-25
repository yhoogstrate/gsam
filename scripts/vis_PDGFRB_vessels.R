#!/usr/bin/env R

data <- readRDS('data/gsam/IF/Experiments revisions/COL1A1_PDGFRB_CD31_GFAP/rerun analysis/measurement_CD31_Collagen_PDGFRB.RDS')


## full intensities per tile ----

parse_csv <- function (fn) {
  dat <- read.csv(fn) |> 
    dplyr::mutate(fn2 = fn)
  
  return(dat)
}


files.CD31 <- data.frame(fn=Sys.glob("data/2022 11 14 Mean modes/*/*Mode_Endothelaial.csv")) |> 
  dplyr::mutate(data = pbapply::pblapply(fn ,parse_csv)) |> 
  tidyr::unnest(data) |> 
  dplyr::rename_with( ~ paste0(.x, ".CD31")) |> 
  dplyr::mutate(basename = gsub("^.+/([^/]+)","\\1",fn.CD31)) |>
  dplyr::mutate(regionname = gsub(".tif_.+$",".tif",basename))
stopifnot(files.CD31$fn.CD31 == files.CD31$fn2.CD31)

files.PDGFRB <- data.frame(fn=Sys.glob("data/2022 11 14 Mean modes/*/*Mode_Pericytes.csv")) |> 
  dplyr::mutate(data = pbapply::pblapply(fn ,parse_csv)) |> 
  tidyr::unnest(data) |> 
  dplyr::rename_with( ~ paste0(.x, ".PDGFRB"))  |> 
  dplyr::mutate(basename = gsub("^.+/([^/]+)","\\1",fn.PDGFRB)) |>
  dplyr::mutate(regionname = gsub(".tif_.+$",".tif",basename))
stopifnot(files.PDGFRB$fn.PDGFRB == files.PDGFRB$fn2.PDGFRB)

stopifnot(files.CD31$regionname %in% files.PDGFRB$regionname)
stopifnot(files.PDGFRB$regionname %in% files.CD31$regionname)



plt <- files.CD31 |> 
  dplyr::left_join(files.PDGFRB, by=c('regionname'='regionname'), suffix=c('','')) |> 
  dplyr::mutate(figname = gsub("region.+$","",regionname)) |> 
  dplyr::mutate(figname = gsub("^Levi batch ","b",figname)) |> 
  dplyr::mutate(figname = gsub(".czi - Scene #","s",figname,fixed=T)) |> 
  dplyr::mutate(figname = gsub(" .+$","",figname)) |> 
  dplyr::filter(figname %in% c("b3s6","b3s7","b3s8","b4s3","b4s4","b4s5","b6s1","b6s5","b6s7","b6s8") == F ) |>  # exclude normals 
  dplyr::left_join(
    openxlsx::read.xlsx('data/2022 11 14 Mean modes/sample_annotation_COL.xlsx'),
    by=c('figname' = 'sample'),
    suffix==c('','') )




stopifnot(length(unique(plt$figname)) == 26)




sids <- plt |> 
  dplyr::pull(sid) |> 
  unique() |> 
  sort()

df <- data.frame()
for(ssid in sids) {
  plt.single <- plt |> 
    dplyr::filter(sid == ssid) |> 
    dplyr::filter(Mean.CD31 > 0 & Mean.PDGFRB > 0)
  
  ggplot(plt.single, aes(x=Mean.CD31, y=Mean.PDGFRB)) +
    geom_point() +
    ggpubr::stat_cor(method = "spearman", aes(label = ..r.label..)) +
    labs(subtitle=unique(plt.single$sid))
  
  
  df <- rbind(df, data.frame(
    sid = unique(plt.single$sid),
    R = cor(plt.single$Mean.CD31, plt.single$Mean.PDGFRB, method="spearman"), # sometimes excessive outliers affecting squared errors
    type = "tumor"
  ) |> 
    dplyr::mutate(    resection = ifelse(grepl("1",sid),"primary","recurrence"))
  
  )
}

m = median(df$R)



df <- rbind(
  df |> dplyr::mutate(type="data"),
  df |> dplyr::mutate(R= 0, type="offset") # nice for lines
)



ggplot(df, aes(x=sid, y=R, group=sid)) +
  #facet_grid(cols = vars( resection), scales = "free", space="free") +
  geom_hline(yintercept = m, lwd = 0.6, lty = "dotted", col = "red") +
  annotate(geom="text", x="AAY2", y=0.99, label=paste0("Median: ",round(m,2)),hjust=0)  +
  geom_line() +
  geom_point(data = df |>  dplyr::filter(type == "data")) +
  labs(x=NULL, y = "Correlation PDGFRB vs. CD31") +
  scale_y_continuous(breaks=c(-1,0,1), limits=c(-1,1)) +
  theme_bw()  +
  theme(
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.text.x = element_text(angle = 90, vjust = 0.45,hjust=1) ,
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25)
  )

ggsave("output/figures/2022_figure_S-PDGFRB_x_CD31_cor.pdf", width=8.3 / 4,height=8.3/4, scale=2)







