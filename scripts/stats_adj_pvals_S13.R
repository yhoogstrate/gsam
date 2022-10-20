

data.frame(pval =c(
  0.74, 0.9, 0.64,
  0.16, 0.51, 0.2,
  0.6, 0.75, 0.5,
  0.48, 0.12, 0.13,
  0.3, 0.66, 0.41
)) |> 
  dplyr::mutate(padj = p.adjust(pval, method="fdr")) |> 
  dplyr::mutate(padj.f = format.pval(padj, digits=3))


