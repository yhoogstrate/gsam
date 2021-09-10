#!/usr/bin/env R


# load data ----

source("scripts/R/cnv_matrix.R")


# cnv ----

locs <- cnv_matrix_genes %>%
  dplyr::filter( gene %in% c('PTCH1') ) %>%
  tibble::rownames_to_column('loc') %>%
  dplyr::pull('loc')


plt.r1 <- cnv_matrix %>%
  tibble::rownames_to_column('loc') %>%
  tibble::as_tibble(.name_repair="unique") %>%
  dplyr::select(contains('1') | contains('loc')) %>%
  as.data.frame() %>%
  dplyr::filter(loc %in% locs) %>%
  dplyr::mutate(loc=NULL) %>%
  colMeans()



plt.r2 <- cnv_matrix %>%
  tibble::rownames_to_column('loc') %>%
  tibble::as_tibble(.name_repair="unique") %>%
  dplyr::select(contains('2') | contains('loc')) %>%
  as.data.frame() %>%
  dplyr::filter(loc %in% locs) %>%
  dplyr::mutate(loc=NULL) %>%
  colMeans()


plt <- rbind(data.frame(sid = names(plt.r1), lr = plt.r1, res="R1"),
  data.frame(sid = names(plt.r2), lr = plt.r2, res="R2"))


ggplot(plt, aes(x=res, y= lr)) +
  geom_point() +
  geom_violin(draw_quantiles = c(0.5), col="black") +
  geom_jitter( position=position_jitter(0.2), size=0.9) +
  #labs(x = NULL, col = NULL, y = "Tumour cell percentage" ) +
  geom_signif(
    comparisons = list(c("R1","R2")),
    test="wilcox.test"
    #col="black"
  ) + 
  job_gg_theme 



                      