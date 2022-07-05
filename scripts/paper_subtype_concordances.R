#!/usr/bin/env R


source('scripts/R/palette.R')
source("scripts/R/gsam_metadata.R")



# giga.RNA.quantity.(ng)"                      "giga.nM.on.qPCR"                             "tumour.percentage.dna"                      
# "gliovis.svm_call"                            "gliovis.knn_call"                            "gliovis.gsea_call"                          
# "gliovis.equal_call"                          "gliovis.majority_call"                       "svm.Classical"                              
# "svm.Mesenchymal"                             "svm.Proneural"                               "knn.Classical"                              
# "knn.Mesenchymal"                             "knn.Proneural"                               "ssGSEA.Proneural.score"                     
# "ssGSEA.Classical.score"                      "ssGSEA.Mesenchymal.score"                    "ssGSEA.Proneural"              
# "ssGSEA.Classical"                            "ssGSEA.Mesenchymal"                          "NMF.123456.PCA.SVM.class"                   
# "NMF.123456.PCA.SVM.Classical.p"              "NMF.123456.PCA.SVM.Proneural.p"              "NMF.123456.PCA.SVM.Mesenchymal.p"           
# "NMF:123456.1"                                "NMF:123456.2"                                "NMF:123456.3"    

plt <- gsam.rna.metadata %>% 
  dplyr::mutate(gliovis.match = ifelse(gliovis.majority_call == NMF.123456.PCA.SVM.class,"Concordant","Discordant")) %>% 
  dplyr::filter(!is.na(NMF.123456.PCA.SVM.class) & !is.na(gliovis.majority_call) ) %>% 
  
  dplyr::select(sid, gliovis.match, gliovis.majority_call , NMF.123456.PCA.SVM.class,
                "gliovis.svm_call",
                "gliovis.knn_call",
                "gliovis.gsea_call" 
                ) %>% 
  dplyr::mutate(order = rank( NMF.123456.PCA.SVM.class)) %>%  # , gliovis.majority_call, gliovis.svm_call, gliovis.gsea_call, gliovis.knn_call
  tidyr::pivot_longer(cols = c(
    'gliovis.majority_call' , 'NMF.123456.PCA.SVM.class',
    
    "gliovis.svm_call",
    "gliovis.knn_call",
    "gliovis.gsea_call"
  )) %>% 
  dplyr::mutate(facet_y = case_when(
    name == "NMF.123456.PCA.SVM.class" ~ 1,
    name == "gliovis.majority_call" ~ 2,
    TRUE ~ 3
  )) %>% 
  dplyr::mutate(name = ifelse(name == "NMF.123456.PCA.SVM.class", "GITS SVM", name))%>% 
  dplyr::mutate(name = ifelse(name == "gliovis.majority_call", "GlioVis (majority call)", name)) %>% 
  dplyr::mutate(name = ifelse(name == "gliovis.svm_call", "GlioVis sub (SVM)", name)) %>% 
  dplyr::mutate(name = ifelse(name == "gliovis.knn_call", "GlioVis sub (KNN)", name)) %>% 
  dplyr::mutate(name = ifelse(name == "gliovis.gsea_call", "GlioVis sub (GSEA)", name))


ggplot(plt, aes(x = reorder(sid, order), y=name, fill=value) ) +
  geom_tile(col="black") +
  facet_grid(cols = vars(gliovis.match), rows=vars(facet_y), scales = "free", space="free") +
  theme_bw() +
  scale_fill_manual(values = subtype_colors)

ggsave("output/figures/subtype_concordance.svg",width=12,height=1.6)
ggsave("output/figures/subtype_concordance.pdf",width=12,height=1.6)



