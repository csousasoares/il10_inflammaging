## Load Packages and Set Seed --------------------------------------------------

set.seed(123)

library(tidyverse)
library(ComplexHeatmap)

## Load data from the Vadim Gladyshev Lab obtained from the
## UK Biobank (thousands of samples across ages)

vadim_proteome <- read.csv("input_data\\vadim_proteome.csv", 
                           header = T, sep = ";")

il_10_human <- c("CCL3",
                 "CCL4",
                 "CCL5",
                 "IL10",
                 "GZMK",
                 "GZMB",
                 "CD38",
                 "EOMES",
                 "LAG3",
                 "TOX2",
                 "PDCD1",
                 "TIGIT",
                 "PRF1",
                 "CDKN1A",
                 "CDKN2A",
                 "CCNA2") ## Get the Human IL-10 signature

vadim_il10sig <- vadim_proteome %>% 
  dplyr::filter(protein %in% il_10_human)

vadim_il10sig$log10qval <- -log10(vadim_il10sig$q.value)

vadim_il10sig <- vadim_il10sig %>% mutate(log10qval_corr = case_when(
  log10qval == Inf ~ 193, ## Values with Infinite log10qvalue clipped at 193.
  log10qval < 193 ~ log10qval
))

vadim_il10sig_corr <- vadim_il10sig %>% 
  select(protein, correlation.coefficient..r.)

row.names(vadim_il10sig_corr) <- vadim_il10sig_corr$protein

vadim_il10sig_corr <- vadim_il10sig_corr %>% 
  select(correlation.coefficient..r.)


qvals <- vadim_il10sig$q.value ## Get q-values for plotting.

ha = rowAnnotation(log10qval = anno_barplot(qvals), width = unit(2, "cm"))

col_fun = colorRamp2::colorRamp2(c(-0.4, 0, 0.4), 
                                 c("deepskyblue", "white", "red"))
col_fun(seq(-0.4, 0.4))

sig <- matrix(
  data = qvals,
  nrow = 13,
  ncol = 1
)

heatmap_plasma_proteome <- Heatmap(
        vadim_il10sig_corr,
        col = col_fun,
        row_title = "Plasma Proteome vs Chronological Age (UK Biobank)",
        column_labels = "",
        row_dend_gp = gpar(lwd = 0.5),
        heatmap_legend_param = list(
            title = "Correlation Coefficient", at = c(-0.4, 0, 0.4), 
            labels = c("-0.4", "0", "0.4"),
            border = "black",
            title_position = "leftcenter-rot",
            legend_height = unit(4, "cm"),
            grid_width = unit(0.5, "cm"),
            legend_gp = gpar(lwd = 0.5)),
        border_gp = gpar(col = "black", lty = 1, lwd = 0.5),
        rect_gp = gpar(col = "white", lwd = 2),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(sig[i, j] < 0.0001) {
            grid.text("****", x, y)
          } else if(sig[i, j] < 0.001) {
            grid.text("***", x, y)
          } else if(sig[i, j] < 0.01) {
            grid.text("**", x, y)
          } else if(sig[i, j] < 0.05) {
            grid.text("*", x, y)
          }
        }
        ) ## Heatmap with asterisks for q-values.
        
heatmap_plasma_proteome ## Save image manually



