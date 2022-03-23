#!/usr/bin/env R

# load data ----

source("scripts/R/job_gg_theme.R")
source("scripts/R/youri_gg_theme.R")

source("scripts/R/gsam_metadata.R")


# plt ----


ids.ecm <- c('AAJ','AAY','AOA','GAH')
ids.neuron <- c('AZF','AZG','CDA','CDD','JAK','HAE')

plt <- data.frame(pid = c(ids.neuron, ids.ecm)) %>% 
  dplyr::mutate(staining = ifelse(pid %in% ids.neuron,"neuron","ECM")) %>% 
  dplyr::mutate(tumour.marker = ifelse(F,"","GFAP")) %>% 
  dplyr::left_join(
    gsam.rna.metadata %>% 
      # endothelial.component, C4.signature, C5.signature
      dplyr::select(sid, pid, resection, neuron.component, oligodendrocyte.component, extracellular.matrix.component, blacklist.pca) %>% 
      dplyr::filter(blacklist.pca == F)
    , by=c('pid'='pid')
  ) %>% 
  dplyr::mutate(neuron.component = -1 * neuron.component) # pca inverted this component


rm(ids.neuron, ids.ecm)


plt <- plt %>%
  tidyr::pivot_longer(cols=c('neuron.component','extracellular.matrix.component')) %>% 
  dplyr::rename(signature.score = value) %>% 
  dplyr::rename(signature = name)



ggplot(plt, aes(x=reorder(pid, sid),y = signature.score)) +
  geom_path( arrow = arrow(ends = "last", type = "closed", angle=15, length = unit(0.12, "inches"))  ) +
  facet_grid(rows = vars(signature), cols=vars(staining), space = "free_x", scales="free") +
  theme_bw()


