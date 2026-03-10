#!/usr/bin/env R

# data ----

#source('scripts/load_G-SAM_metadata.R')
#meta.p <- gsam.patient.metadata |> 
#  dplyr::select(studyID, survivalDays, survival.events)
#write.csv(meta.p , "/tmp/meta.p.csv")

meta.p <- read.csv("/tmp/meta.p.csv")


data = readxl::read_xlsx('zip/vst_primary_cat_exp2_categorization_only.xlsx') |> 
  tibble::column_to_rownames('...1') |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column('sample_id') |> 
  dplyr::mutate(patient_id = gsub("^(...).+$","\\1",sample_id)) |> 
  dplyr::mutate(resection = gsub("^...(.).*?$","\\1",sample_id)) |> 
  dplyr::mutate(Category = factor(Category, levels=c("Other", "Double_High")))



# primaries ----



km <- data |> 
  dplyr::left_join(meta.p, by=c('patient_id'='studyID'))



library(survival)
library(survminer)

fit_overall <- surv_fit(Surv(survivalDays, survival.events) ~ 1, data = km)
fit_group <- surv_fit(Surv(survivalDays, survival.events) ~ Category, data = km)


## ggsurvplot export bugfix - https://github.com/kassambara/survminer/issues/152#issuecomment-938941051
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}



km_plot.curve <- survminer::ggsurvplot(
  fit_group,                    # The fitted survival object
  data = km,                  # Specify the original data
  risk.table = TRUE,            # Add a risk table below the plot
  pval = TRUE,                  # Add p-value (from log-rank test) to the plot
  conf.int = TRUE,              # Add confidence intervals
  xlab = "Time in Days",        # X-axis label
  legend.title = "Caetegory", # Legend title
  #fontsize = 0.25,
  ggtheme = theme_light()       # Apply a ggplot theme,
)


km_plot.risk_table <- survminer::ggsurvtable(fit_group,
                                             data = km, 
                                             tables.theme = theme_light(), 
                                             ggtheme = theme_light(),
                                             fontsize = 3.14,
                                             #font.family = theme_cellpress_font_family,
                                             xlab = "OS from primary (Days)")



# hacky stuff needed for neat gg compatible exports, and nearly fit the panels
patchwork:::`/.ggplot`(km_plot.curve$plot, km_plot.risk_table$risk.table) + patchwork::plot_layout(heights = c(4,1))


ggsave("output/figures/analysis_survival_prognostic__km.pdf", width=5.5,height=3.5)


