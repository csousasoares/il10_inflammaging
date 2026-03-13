## Load Packages and Set Seed --------------------------------------------------

set.seed(123)

library(tidyverse)
library(clusterProfiler)
library(scales)

## Load Data -------------------------------------------------------------------

## Load .csv with info regarding changes in genes for individual 
## immune cell types in aging. Data from "Multimodal profiling reveals 
## tissue-directed signatures of human immune cells altered with age"
## doi: 10.1038/s41590-025-02241-4

aging_change <- read.csv("input_data\\aging_changes.csv", 
                         sep = ",", header = T)

unique(aging_change$assay)
unique(aging_change$tissue_group) ## Check tissues and genes

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


## Inspect individual organ info:

il_10_aging <- aging_change %>% dplyr::filter(id %in% il_10_human)

il_10_aging_bm <- il_10_aging %>% dplyr::filter(tissue_group == "bma") 

il_10_aging_blood <- il_10_aging %>% dplyr::filter(tissue_group == "blo") 

il_10_aging_spleen <- il_10_aging %>% dplyr::filter(tissue_group == "spl") 

head(il_10_aging_bm)
head(il_10_aging_blood)
head(il_10_aging_spleen)


## GSEA BM ---------------------------------------------------------------------

il_10_gs <- data.frame(
  gs_name = "IL10_sig",
  gene_symbol = il_10_human
)

colnames(aging_change)

aging_change <- aging_change %>% dplyr::filter(aveexpr > 2.5) 

## Filtered genes with very low expression.

aging_change_bm <- aging_change %>% 
  dplyr::filter(tissue_group == "bma")

aging_change_bm <- aging_change_bm %>% 
  dplyr::mutate(direction = case_when(
  logfc < 0 ~ -1,
  logfc >= 0 ~ 1
)) %>% dplyr::mutate(log10pval = -log10(p_value)) %>% 
  dplyr::mutate(metric = logfc)

aging_change_bm <- aging_change_bm %>% 
  group_by(assay) %>% dplyr::arrange(desc(metric))

head(aging_change_bm)

aging_change_bm <- group_split(aging_change_bm)


results_bm <- lapply(aging_change_bm, function(df) {
  df <- as.data.frame(df)
  
  # Proper numeric named vector
  ranks <- as.numeric(df[["metric"]])
  names(ranks) <- as.character(df[["id"]])
  ranks <- sort(ranks, decreasing = TRUE)
  
  gsea <- GSEA(ranks,
               TERM2GENE = il_10_gs,
               verbose = FALSE,
               seed = TRUE,
               pvalueCutoff = 1,
               minGSSize = 0,
               pAdjustMethod = "BH",
               by = "fgsea")
  
  # Convert to data frame and keep NES etc.
  out <- as.data.frame(gsea)[, c("ID", "NES", "pvalue", "p.adjust")]
  out$assay <- unique(df$assay)
  # Add an identifier for this split (e.g. assay type)
  
  return(out)
})

bm_gsea_all <- do.call(rbind.data.frame, results_bm)
bm_gsea_all$tissue <- "bone_marrow"


## GSEA Blood ------------------------------------------------------------------

aging_change_blood <- aging_change %>% 
  dplyr::filter(tissue_group == "blo")

aging_change_blood <- aging_change_blood %>% 
  dplyr::mutate(direction = case_when(
  logfc < 0 ~ -1,
  logfc >= 0 ~ 1
)) %>% dplyr::mutate(log10pval = -log10(p_value)) %>% 
  dplyr::mutate(metric = logfc)

aging_change_blood <- aging_change_blood %>% 
  group_by(assay) %>% dplyr::arrange(desc(metric))

head(aging_change_blood)

aging_change_blood <- group_split(aging_change_blood)


results_blood <- lapply(aging_change_blood, function(df) {
  df <- as.data.frame(df)
  
  # Proper numeric named vector
  ranks <- as.numeric(df[["metric"]])
  names(ranks) <- as.character(df[["id"]])
  ranks <- sort(ranks, decreasing = TRUE)
  
  gsea <- GSEA(ranks,
               TERM2GENE = il_10_gs,
               verbose = FALSE,
               seed = TRUE,
               pvalueCutoff = 1,
               minGSSize = 0,
               pAdjustMethod = "BH",
               by = "fgsea")
  
  # Convert to data frame and keep NES etc.
  out <- as.data.frame(gsea)[, c("ID", "NES", "pvalue", "p.adjust")]
  out$assay <- unique(df$assay)
  # Add an identifier for this split (e.g. assay type)
  
  return(out)
})

blood_gsea_all <- do.call(rbind.data.frame, results_blood)
blood_gsea_all$tissue <- "blood"


## GSEA Spleen -----------------------------------------------------------------

aging_change_spleen <- aging_change %>% dplyr::filter(tissue_group == "spl")
aging_change_spleen <- aging_change_spleen %>% 
  dplyr::mutate(direction = case_when(
  logfc < 0 ~ -1,
  logfc >= 0 ~ 1
)) %>% dplyr::mutate(log10pval = -log10(p_value)) %>% 
  dplyr::mutate(metric = logfc)

aging_change_spleen <- aging_change_spleen %>% 
  group_by(assay) %>% dplyr::arrange(desc(metric))

head(aging_change_spleen)

aging_change_spleen <- group_split(aging_change_spleen)


results_spleen <- lapply(aging_change_spleen, function(df) {
  df <- as.data.frame(df)
  
  # Proper numeric named vector
  ranks <- as.numeric(df[["metric"]])
  names(ranks) <- as.character(df[["id"]])
  ranks <- sort(ranks, decreasing = TRUE)
  
  gsea <- GSEA(ranks,
               TERM2GENE = il_10_gs,
               verbose = FALSE,
               seed = TRUE,
               pvalueCutoff = 1,
               minGSSize = 0,
               pAdjustMethod = "BH",
               by = "fgsea")
  
  # Convert to data frame and keep NES etc.
  out <- as.data.frame(gsea)[, c("ID", "NES", "pvalue", "p.adjust")]
  out$assay <- unique(df$assay)
  # Add an identifier for this split (e.g. assay type)
  
  return(out)
})

spleen_gsea_all <- do.call(rbind.data.frame, results_spleen)

spleen_gsea_all$tissue <- "spleen"


## Merge GSEA Results ----------------------------------------------------------


all_cells <- do.call("rbind", 
                     list(blood_gsea_all, bm_gsea_all, spleen_gsea_all)) %>% 
  dplyr::filter(p.adjust < 0.1) %>% 
  dplyr::mutate(log10padj = -log10(p.adjust))

all_cells <- all_cells %>% 
  dplyr::group_by(assay) %>% dplyr::mutate(n = n())


ggplot(all_cells, aes(x = tissue, y = forcats::fct_reorder(assay, NES))) +
  geom_point(shape = 21, aes(fill = NES, size = p.adjust)) +
  scale_fill_gradient(low = "blue", high = "red", limits = c(-2.5,2.5)) +
  scale_size_continuous(range = c(10, 5), transform = "log10") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "GSEA IL-10 Signature - Human (padjust < 0.1)",
       y = "Cell Types", x = "Tissues")

ggsave(
  "il10_sig_16_genes_gsea_immune_cells_across_tissues_padj_01_new_gradient.pdf",
  path = "Plots",
  width = 5,
  height = 7,
  unit = "in"
)

head(all_cells, 10)

writeLines(capture.output(sessionInfo()), "sessionInfo_human_data.txt")
