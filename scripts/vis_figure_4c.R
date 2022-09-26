#!/usr/bin/env R


# load data ----


source('scripts/load_results.out.R')
source('scripts/load_G-SAM_expression_data.R')
source('scripts/load_G-SAM_metadata.R')


## select genes / vars ----


tmp.labels <- results.out |> 
  dplyr::filter(!is.na(gid)) |> 
  dplyr::filter(hugo_symbol %in% c('CREB5', 'COA1', 'AHCYL1') == F) |>

  dplyr::mutate(primary.marker.genes = ifelse(primary.marker.genes %in% c('PN subtype', 'CL subtype', 'MES subtype'), NA, primary.marker.genes) ) |>
  
  dplyr::mutate(`pericyte` = hugo_symbol %in% c("RGS5","PDGFRB","CD248","HEYL","CFH")) |> 

  dplyr::filter(!is.na(primary.marker.genes) | C1.2022 | pericyte) |>

  dplyr::select(.data$gid, 
                .data$hugo_symbol, 
                
                .data$C1.2022,
                .data$primary.marker.genes,
                .data$pericyte
                ) |>
  
  dplyr::bind_rows(c(gid="NMF:H[1] (CL)", hugo_symbol = "NMF:H[1] (CL)", primary.marker.genes = 'CL subtype')) |>
  dplyr::bind_rows(c(gid="NMF:H[2] (PN)",  hugo_symbol = "NMF:H[2] (PN)",  primary.marker.genes = 'PN subtype')) |>
  dplyr::bind_rows(c(gid="NMF:H[3] (MES)",  hugo_symbol = "NMF:H[3] (MES)",  primary.marker.genes = 'MES subtype')) |>
  dplyr::bind_rows(c(gid='tumor-% DNA',    hugo_symbol = 'tumor-% DNA',    primary.marker.genes = 'tumor-% DNA')) |>
  
  dplyr::mutate(pericyte = ifelse(is.na(pericyte), F, pericyte)) |> 
  dplyr::mutate(`C1.2022` = ifelse(is.na(C1.2022), F, C1.2022)) |> 
  
  dplyr::mutate(decoy = T) |> 
  tidyr::pivot_wider(names_from = `primary.marker.genes`, values_from = decoy, values_fill = FALSE) |> 
  dplyr::mutate(`NA` = NULL) |> 
  as.data.frame() |> 
  
  dplyr::rename(`NMF:H[1] (CL)` = `CL subtype`) |> 
  dplyr::rename(`NMF:H[2] (PN)` = `PN subtype`) |> 
  dplyr::rename(`NMF:H[3] (MES)` = `MES subtype`) |> 
  dplyr::rename(`OD` = `oligodendrocyte`) |> 
  dplyr::rename(`NE` = `neuron`) |> 
  dplyr::rename(`AC` = `astrocyte`) |> 
  dplyr::rename(`EN` = `endothelial`) |>
  dplyr::rename(`PE` = `pericyte`) |> 
  dplyr::rename(`TIL/TC` = `TIL / T-cell`) |> 
  dplyr::rename(`TAM/MG` = `microglia/TAM`) |> 
  dplyr::rename(`C1 (ECM)` = `C1.2022`) |> 
  dplyr::rename(`T% (DNA)` = `tumor-% DNA`) |> 
  dplyr::rename(`chr7 (T)` = `chr7 gained (tumor)`)


## select samples ----


tmp.metadata <- gsam.rna.metadata |> 
  dplyr::filter(blacklist.pca == F) |> 
  dplyr::filter(pat.with.IDH == F) |> 
  dplyr::filter(sid %in% c('BAI2', 'CAO1-replicate', 'FAB2', 'GAS2-replicate') == F ) |> 
  dplyr::filter(tumour.percentage.dna >= 15)



## select pairs ----


tmp.metadata.paired <- tmp.metadata %>%
  dplyr::filter(pid %in% 
                  (tmp.metadata %>%
                     dplyr::group_by(pid) %>%
                     dplyr::tally() %>%
                     dplyr::filter(n == 2) %>% 
                     dplyr::ungroup() %>%
                     dplyr::filter(!duplicated(pid)) %>%
                     dplyr::pull(pid))
  ) %>%
  dplyr::mutate(pid = as.factor(as.character(pid))) # re-factor


tmp.expression.paired <- gsam.gene.expression.all.vst |> 
  as.data.frame() |> 
  dplyr::select(tmp.metadata.paired$sid)
stopifnot(colnames(tmp.expression.paired) == tmp.metadata.paired$sid)



## select data ----


plt <- tmp.labels |>
  dplyr::filter(gid != 'tumor-% DNA') |>
  dplyr::filter(hugo_symbol %in% c('NMF:H[1] (CL)', 'NMF:H[2] (PN)', 'NMF:H[3] (MES)') == F) |> # lines below add the NMF matrices
  dplyr::select('gid','hugo_symbol') |>
  dplyr::left_join(tmp.expression.paired |> 
                     tibble::rownames_to_column('gid'),
                     by=c('gid'='gid')) |> 
  dplyr::mutate(hugo_symbol = NULL) |>
  tibble::column_to_rownames('gid') |>
  t() |> 
  as.data.frame() |>
  tibble::rownames_to_column('sid') |>
  dplyr::left_join(
    
    gsam.rna.metadata |>
      dplyr::select(sid, `NMF:150:1`, `NMF:150:2`, `NMF:150:3`, tumour.percentage.dna) |> 
      dplyr::rename(`NMF:H[1] (CL)` = `NMF:150:1`) |> 
      dplyr::rename(`NMF:H[2] (PN)` = `NMF:150:2`) |> 
      dplyr::rename(`NMF:H[3] (MES)` = `NMF:150:3`) |> 
      dplyr::rename(`tumor-% DNA` = tumour.percentage.dna)

    , by=c('sid' = 'sid')
  ) |>
  tibble::column_to_rownames('sid') |>
  t() |>
  as.data.frame()



# Add NMF 1 - 3
# Add ECM signature
# Add EPIC macrophage score [not necessary]

odd <- 1:ncol(plt) |> purrr::keep(~ . %% 2 == 1)
even <- 1:ncol(plt) |> purrr::keep(~ . %% 2 == 0)



plt.r1 <- plt[,odd] 
plt.r2 <- plt[,even]
stopifnot ( gsub("^(...).*$","\\1",colnames(plt.r1)) == gsub("^(...).*$","\\1",colnames(plt.r2)) )
rm(plt)


plt <- log2(plt.r1 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) /
            plt.r2 %>% `colnames<-`(gsub("^(...).*$","\\1",colnames(.))) ) |> 
  tibble::rownames_to_column('gid') |> 
  dplyr::left_join(
    tmp.labels |>
      dplyr::select(hugo_symbol, gid), 
    by = c('gid' = 'gid'),
    suffix = c('', '')) |> 
  dplyr::mutate(gid=NULL) |> 
  tibble::column_to_rownames('hugo_symbol')


stopifnot(tmp.labels$hugo_symbol == rownames(plt))




tmp.labels <- tmp.labels |> 
  dplyr::mutate(gid = NULL) |>
  tibble::column_to_rownames('hugo_symbol')
  
tmp.labels <- tmp.labels |> 
  dplyr::select(`T% (DNA)`,
                `chr7 (T)`,
                `EN`,
                OD,
                `NMF:H[2] (PN)`,
                NE,
                `NMF:H[1] (CL)`,
                AC,
                `C1 (ECM)`,
                `NMF:H[3] (MES)`,
                `PE`,
                `TIL/TC`,
                `TAM/MG`
                ) |>  # order
  dplyr::mutate_all(function(arg) { return (ifelse(arg, as.logical(NA), arg)) }) # gray labels


# plot & export ----


recursiveCorPlot::recursiveCorPlot(plt, tmp.labels, 6 * 0.8, 1 , caption=paste0("G-SAM: n=",ncol(plt)," pairs"))
ggsave("output/figures/2022_figure_4c___logFc_gene_per_pair.pdf",width=8.3 / 2,height=8.3/2, scale=2.1)



