### FOR GITHUB
# ============================================================
# Differential expression data loading and harmonization
# ============================================================

library(dplyr)
library(readxl)

# ------------------------------------------------------------
# Helper function
# ------------------------------------------------------------
prepare_deg_table <- function(data, dataset_name) {
  data %>%
    dplyr::select(gene_name, log2FC, pvalAdj) %>%
    dplyr::mutate(Dataset = dataset_name)
}

# ------------------------------------------------------------
# Load input data
# ------------------------------------------------------------

CD4_pMT10_plus_vs_minus_Zn <- read.csv(
  "CD4_pMT10_plus_vs_minus_Zn_GSEA.csv",
  header = TRUE,
  sep = ","
)

CD8_pMT10_plus_vs_minus_Zn <- read.csv(
  "CD8_pMT10_plus_vs_minus_Zn_GSEA.csv",
  header = TRUE,
  sep = ","
)

CD4_pMT10IFNgKO_plus_vs_minus_Zn <- read.csv(
  "CD4_pMT10IFNgKO_plus vs minus_Zn_GSEA_July2025.csv",
  header = TRUE,
  sep = ","
)

CD8_pMT10IFNgKO_plus_vs_minus_Zn <- read.csv(
  "CD8_pMT10IFNgKO_plus vs minus_Zn_GSEA_July2025.csv",
  header = TRUE,
  sep = ","
)


# ------------------------------------------------------------
# Standardize DEG tables
# ------------------------------------------------------------

DEGs_CD4_pMT10_plus_vs_minus_Zn <- prepare_deg_table(
  CD4_pMT10_plus_vs_minus_Zn,
  "CD4_pMT10_plus_vs_minus_Zn"
)

DEGs_CD8_pMT10_plus_vs_minus_Zn <- prepare_deg_table(
  CD8_pMT10_plus_vs_minus_Zn,
  "CD8_pMT10_plus_vs_minus_Zn"
)


DEGs_CD4_pMT10IFNgKO_plus_vs_minus_Zn <- prepare_deg_table(
  CD4_pMT10IFNgKO_plus_vs_minus_Zn,
  "CD4_pMT10IFNgKO_plus_vs_minus_Zn"
)

DEGs_CD8_pMT10IFNgKO_plus_vs_minus_Zn <- prepare_deg_table(
  CD8_pMT10IFNgKO_plus_vs_minus_Zn,
  "CD8_pMT10IFNgKO_plus_vs_minus_Zn"
)


# Respective to Figure 4e and S4b 
# ============================================================
# GSEA for IL-10-related signature
# ============================================================

library(dplyr)
library(tibble)
library(fgsea)
library(ggplot2)
library(scales)

# ------------------------------------------------------------
# IL-10 signature
# ------------------------------------------------------------
selected_genes <- c(
  "Ccl3", "Ccl4", "Ccl5", "Il10", "Gzmk", "Gzmb", "Cd38",
  "Eomes", "Lag3", "Tox2", "Pdcd1", "Tigit", "Prf1",
  "Cdkn1a", "Cdkn2a", "Ccna2"
)

combined_interventions <- bind_rows(
  DEGs_CD4_pMT10_plus_vs_minus_Zn,
  DEGs_CD4_pMT10IFNgKO_plus_vs_minus_Zn,
  DEGs_CD8_pMT10_plus_vs_minus_Zn,
  DEGs_CD8_pMT10IFNgKO_plus_vs_minus_Zn
)


combined_aging_core <- bind_rows(
  aging_sig_Mm_AT,
  reprog_Mm_Hs_AT,
  Interventions_Max_lifespan_AT,
  Interventions_Median_lifespan_AT,
  Aging_Liver_AT
) %>%
  mutate(
    Dataset = case_when(
      Dataset == "aging_sig_Mm_AT" ~ "Aging Sig",
      Dataset == "reprog_Mm_Hs_AT" ~ "Reprogramming Sig",
      Dataset == "Interventions_Max_lifespan_AT" ~ "Interv Max Lifespan Sig",
      Dataset == "Interventions_Median_lifespan_AT" ~ "Interv Median Lifespan Sig",
      Dataset == "Aging_Liver_AT" ~ "Aging Liver Sig",
      TRUE ~ Dataset
    )
  )


# ------------------------------------------------------------
# GSEA helper function
# ------------------------------------------------------------
run_signature_gsea <- function(combined_data, selected_genes) {
  gsea_results_list <- list()
  plots_list <- list()
  
  for (dataset_name in unique(combined_data$Dataset)) {
    message("Running GSEA for: ", dataset_name)
    
    dataset_data <- combined_data %>%
      filter(Dataset == dataset_name)
    
    # Positive enrichment
    dataset_positive <- dataset_data %>%
      mutate(log2FC = as.numeric(log2FC)) %>%
      filter(log2FC > 0, !is.na(log2FC), is.finite(log2FC)) %>%
      arrange(log2FC)
    
    if (nrow(dataset_positive) > 0) {
      ranks_pos <- deframe(dataset_positive %>% select(gene_name, log2FC))
      
      gsea_result_pos <- tryCatch(
        {
          fgseaMultilevel(
            pathways = list(selected_genes = selected_genes),
            stats = ranks_pos,
            minSize = 5,
            maxSize = 500,
            scoreType = "pos"
          ) %>%
            mutate(Direction = "Positive", Dataset = dataset_name)
        },
        error = function(e) {
          message("GSEA positive failed for: ", dataset_name, " | ", e$message)
          return(NULL)
        }
      )
      
      if (!is.null(gsea_result_pos)) {
        plot_pos <- plotEnrichment(selected_genes, ranks_pos) +
          labs(title = paste("GSEA:", dataset_name, "- Positive"))
        
        gsea_results_list[[paste0(dataset_name, "_positive")]] <- gsea_result_pos
        plots_list[[paste0(dataset_name, "_plot_pos")]] <- plot_pos
      }
    }
    
    # Negative enrichment
    dataset_negative <- dataset_data %>%
      mutate(log2FC = as.numeric(log2FC)) %>%
      filter(log2FC < 0, !is.na(log2FC), is.finite(log2FC)) %>%
      arrange(desc(log2FC))
    
    if (nrow(dataset_negative) > 0) {
      ranks_neg <- deframe(dataset_negative %>% select(gene_name, log2FC))
      
      gsea_result_neg <- tryCatch(
        {
          fgseaMultilevel(
            pathways = list(selected_genes = selected_genes),
            stats = ranks_neg,
            minSize = 5,
            maxSize = 500,
            scoreType = "neg"
          ) %>%
            mutate(Direction = "Negative", Dataset = dataset_name)
        },
        error = function(e) {
          message("GSEA negative failed for: ", dataset_name, " | ", e$message)
          return(NULL)
        }
      )
      
      if (!is.null(gsea_result_neg)) {
        plot_neg <- plotEnrichment(selected_genes, ranks_neg) +
          labs(title = paste("GSEA:", dataset_name, "- Negative"))
        
        gsea_results_list[[paste0(dataset_name, "_negative")]] <- gsea_result_neg
        plots_list[[paste0(dataset_name, "_plot_neg")]] <- plot_neg
      }
    }
  }
  
  gsea_summary <- bind_rows(gsea_results_list)
  
  if (nrow(gsea_summary) == 0) {
    warning("No GSEA results were generated.")
    return(list(summary = NULL, plots = plots_list))
  }
  
  # Standardize adjusted p-value column name
  names(gsea_summary)[names(gsea_summary) == "padj"] <- "pvalAdj"
  
  # Add summary metrics for plotting
  gsea_summary <- gsea_summary %>%
    mutate(
      Significant = pvalAdj < 0.1,
      absNES = abs(NES),
      color_value = -log10(pvalAdj),
      label = Dataset,
      is_signif = ifelse(Significant, "Significant", "Not significant")
    )
  
  return(list(summary = gsea_summary, plots = plots_list))
}

# ------------------------------------------------------------
# Example: run one analysis set
# ------------------------------------------------------------
# Choose the dataset collection relevant for the panel you want to plot.
# Example shown here uses IFNgKO-related immune comparisons.
gsea_output <- run_signature_gsea(
  combined_data = combined_interventions,
  selected_genes = selected_genes
)

gsea_summary <- gsea_output$summary
plots_list <- gsea_output$plots

# Use only one ordering block depending on the analysis being plotted.

# Example ordering for IFNgKO-related immune comparisons
if (!is.null(gsea_summary)) {
  gsea_summary$label <- factor(
    gsea_summary$label,
    levels = c(
      "CD8_pMT10IFNgKO_plus_vs_minus_Zn",
      "CD8_pMT10_plus_vs_minus_Zn",
      "CD4_pMT10IFNgKO_plus_vs_minus_Zn",
      "CD4_pMT10_plus_vs_minus_Zn"
    )
  )
}

# Alternative ordering for aging-related signatures:
# gsea_summary$label <- factor(
#   gsea_summary$label,
#   levels = c(
#     "Aging Sig",
#     "Reprogramming Sig",
#     "Interv Max Lifespan Sig",
#     "Interv Median Lifespan Sig",
#     "Aging Liver Sig"
#   )
# )

# ------------------------------------------------------------
# Plot : continuous color scale by -log10(adjusted p-value)
# ------------------------------------------------------------
if (!is.null(gsea_summary)) {
  p_gsea_padj <- ggplot(gsea_summary, aes(y = label, x = 0)) +
    geom_segment(
      aes(xend = NES, yend = label, color = color_value),
      size = 1,
      lineend = "round"
    ) +
    geom_point(
      aes(x = NES, color = color_value),
      size = 3
    ) +
    scale_color_gradient(
      name = "-log10(padj)",
      low = "lightblue",
      high = "darkblue"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    xlim(-2.5, 2.5) +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      title = "GSEA for IL-10 Signature"
    ) +
    theme_bw() +
    theme(
      panel.border = element_rect(color = "black", size = 1.5),
      axis.text = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold")
    )
}




# Respective to Figure 4c,d and S5c,d
# ============================================================
# IL-10 signature and age-related tissue expression
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------
# IL-10 signature
# ------------------------------------------------------------
selected_genes <- c(
  "Ccl3", "Ccl4", "Ccl5", "Il10", "Gzmk", "Gzmb", "Cd38",
  "Eomes", "Lag3", "Tox2", "Pdcd1", "Tigit", "Prf1",
  "Cdkn1a", "Cdkn2a", "Ccna2"
)
selected_genes1 <- c("Ilra")
selected_genes2 <- c("Ilrb")

# ------------------------------------------------------------
# Load bulk tissue aging counts
# Data extracted from: Ageing Hallmarks Exhibit Organ-Specific Temporal Signatures. doi:10.1038/s41586-020-2499-y.
# ------------------------------------------------------------
raw_counts_long <- read.csv("Bulk_Tissue_Aging_Counts.csv", header = TRUE)

# Log2-transform counts - normalization step 
raw_counts_long <- raw_counts_long %>%
  mutate(Count = log2(Count + 1))

# ============================================================
# PART A. Correlation between IL-10 signature and Il10 expression
# ============================================================

# ------------------------------------------------------------
# Prepare IL-10 signature data
# ------------------------------------------------------------
il10_signature_long <- raw_counts_long %>%
  filter(Gene %in% selected_genes)

# Order used in plots/facets
tissue_order_all <- c(
  "BAT", "MAT", "Brain", "Lung", "Kidney", "Limb_Muscle",
  "Skin", "WBC", "Heart", "Bone", "Small_Intestine", "SCAT", "Liver"
)

# Summarize mean expression by tissue, dataset, and age
il10_signature_grouped <- il10_signature_long %>%
  group_by(Tissue, Dataset, Age) %>%
  summarise(mean_count = mean(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Dataset, values_from = mean_count) %>%
  drop_na()

# Main tissues for display
main_tissues <- c("Marrow", "Spleen", "Pancreas", "GAT")

il10_signature_main <- il10_signature_grouped %>%
  filter(Tissue %in% main_tissues)

il10_signature_main$Tissue <- factor(il10_signature_main$Tissue, levels = main_tissues)

# Correlation statistics by tissue
cor_results_il10_signature <- il10_signature_main %>%
  group_by(Tissue) %>%
  summarise(
    cor = cor(Il10, `Il10 Signature`, method = "pearson"),
    p_value = cor.test(Il10, `Il10 Signature`, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "r = ", round(cor, 2),
      "\n", "p = ", format.pval(p_value, digits = 2)
    ),
    Tissue = factor(Tissue, levels = main_tissues)
  )

p_il10_signature_correlation <- ggplot(
  il10_signature_main,
  aes(x = `Il10 Signature`, y = Il10)
) +
  geom_point(aes(color = as.factor(Age)), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1) +
  facet_wrap(~ Tissue, scales = "free", ncol = 4) +
  scale_color_discrete(name = "Age (months)") +
  labs(
    x = "IL-10 Signature (mean expression)",
    y = "Il10",
    title = "Correlation between IL-10 Signature and Il10 Expression"
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", size = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "right"
  ) +
  geom_text(
    data = cor_results_il10_signature,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1,
    vjust = 1.1,
    inherit.aes = FALSE,
    size = 5
  )

print(p_il10_signature_correlation)

ggsave(
  filename = "Il10_vs_Il10_Signature_Correlation_Main.pdf",
  plot = p_il10_signature_correlation,
  width = 14,
  height = 10,
  units = "in"
)

ggsave(
  filename = "Il10_vs_Il10_Signature_Correlation_Supp.pdf",
  plot = p_il10_signature_correlation,
  width = 14,
  height = 10,
  units = "in"
)


# ============================================================
# PART B. Age trends for IL-10 signature / Il10rb / Il10
# ============================================================

# ------------------------------------------------------------
# Prepare datasets
#
signature_data <- raw_counts_long %>%
  filter(Gene %in% selected_genes) %>%
  mutate(Dataset = "Il10 Signature")

il10rb_data <- raw_counts_long %>%
  filter(Gene %in% selected_genes1) %>%
  mutate(Dataset = "Il10ra")

il10_data <- raw_counts_long %>%
  filter(Gene %in% selected_genes2) %>%
  mutate(Dataset = "Il10rb")

combined_signature_gene_data <- bind_rows(
  signature_data,
  il10ra_data,
  il10rb_data
)

# ------------------------------------------------------------
# Correlation of age with IL-10 signature signal
# ------------------------------------------------------------
cor_values_signature_age <- signature_data %>%
  group_by(Tissue) %>%
  summarise(
    cor = cor.test(Age, Count, method = "pearson")$estimate,
    p_value = cor.test(Age, Count, method = "pearson")$p.value,
    .groups = "drop"
  )

# Tissue lists
supp_tissues <- c(
  "Heart", "GAT", "Spleen", "SCAT", "Limb_Muscle",
  "Lung", "Liver", "Bone", "Small Intestine",
  "BAT", "Kidney", "WBC", "Brain"
)

main_tissues <- c("Marrow", "Spleen", "Pancreas", "GAT")

# Main tissue subset
signature_data_main <- signature_data %>%
  filter(Tissue %in% main_tissues)

cor_values_signature_main <- cor_values_signature_age %>%
  filter(Tissue %in% main_tissues)

# Supplementary tissue subset
signature_data_supp <- signature_data %>%
  filter(!Tissue %in% main_tissues)

cor_values_signature_supp <- cor_values_signature_age %>%
  filter(!Tissue %in% main_tissues)

# Set factor order for main tissues
signature_data_main$Tissue <- factor(signature_data_main$Tissue, levels = main_tissues)
cor_values_signature_main$Tissue <- factor(cor_values_signature_main$Tissue, levels = main_tissues)

# ------------------------------------------------------------
# Plot: IL-10 signature across main tissues with age
# ------------------------------------------------------------
p_bulk_tissues_signature <- ggplot(
  signature_data_main,
  aes(x = Age, y = Count, color = Tissue)
) +
  geom_smooth(method = "loess", se = TRUE, size = 3, color = "darkblue") +
  facet_wrap(~ Tissue, scales = "free_y", ncol = 4) +
  labs(
    title = "IL-10 Signature in Aging Mouse Tissues",
    x = "Age (Months)",
    y = "Normalized counts - Log2"
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", size = 1.5),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  geom_text(
    data = cor_values_signature_main,
    aes(
      x = Inf,
      y = Inf,
      label = paste(
        "r = ", round(cor, 2),
        "\n p = ", format.pval(p_value)
      )
    ),
    hjust = 1.1,
    vjust = 1.1,
    color = "black",
    size = 5,
    inherit.aes = FALSE
  )

print(p_bulk_tissues_signature)

ggsave(
  filename = "Il10_Signature_Across_Tissues_Age.pdf",
  plot = p_bulk_tissues_signature,
  width = 14,
  height = 10,
  units = "in"
)

ggsave(
  filename = "Il10_Across_Tissues_Age_Main.pdf",
  plot = p_bulk_tissues_signature,
  width = 14,
  height = 10,
  units = "in"
)

ggsave(
  filename = "Il10_Across_Tissues_Age_Supp.pdf",
  plot = p_bulk_tissues_signature,
  width = 14,
  height = 10,
  units = "in"
)

# ============================================================
# PART C. Combined smooth trends by dataset in infiltrated tissues
# ============================================================

combined_main_tissues <- combined_signature_gene_data %>%
  filter(Tissue %in% main_tissues)

p_combined_bulk_tissues <- ggplot(
  combined_main_tissues,
  aes(x = Age, y = Count, color = Dataset)
) +
  geom_smooth(method = "loess", se = TRUE, size = 2.5) +
  facet_wrap(~ Tissue, scales = "free_y") +
  labs(
    title = "IL-10 Signature and Related Genes across Tissues",
    x = "Age (Months)",
    y = "Normalized counts - Log2",
    color = "Signature"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "top"
  )

print(p_combined_bulk_tissues)

# Optional check
unique(combined_signature_gene_data$Tissue)

# ============================================================
# PART D. Multi-gene age trends and tissue-specific correlations
# ============================================================

genes_of_interest <- c("Il10ra", "Il10rb", "Il10")

gene_age_data <- combined_signature_gene_data %>%
  filter(Gene %in% genes_of_interest)

# Correlation statistics by tissue and gene
cor_stats_genes <- gene_age_data %>%
  group_by(Tissue, Gene) %>%
  summarise(
    cor = cor.test(Age, Count)$estimate,
    p = cor.test(Age, Count)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("R = ", round(cor, 2), ", p = ", signif(p, 2))
  )

# Label positions for plot annotations
label_positions <- gene_age_data %>%
  group_by(Tissue) %>%
  summarise(
    x = min(Age) + 0.5,
    y_max = max(Count),
    y_min = min(Count),
    y_range = y_max - y_min,
    .groups = "drop"
  )

cor_stats_genes <- left_join(cor_stats_genes, label_positions, by = "Tissue") %>%
  mutate(
    y = case_when(
      Gene == "Il10" ~ y_max - 0.05 * y_range,
      Gene == "Il10ra" ~ y_max - 0.15 * y_range,
      Gene == "Il10rb" ~ y_max - 0.25 * y_range,
      TRUE ~ y_max
    )
  )

p_gene_age_trends <- ggplot(gene_age_data, aes(x = Age, y = Count, color = Gene)) +
  geom_smooth(method = "loess", se = TRUE, size = 2) +
  facet_wrap(~ Tissue, scales = "free_y") +
  geom_text(
    data = cor_stats_genes,
    aes(x = x, y = y, label = label, color = Gene),
    hjust = 0,
    vjust = 1,
    size = 4,
    inherit.aes = FALSE
  ) +
  labs(
    title = "Il10 and Il10ra/b Age-Related Expression in Infiltrated Tissues",
    x = "Age (Months)",
    y = "Normalized counts - Log2",
    color = "Gene"
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", size = 1.2),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "top"
  )

print(p_gene_age_trends)

ggsave(
  filename = "Il10_and_Il10ra_b_age_related_expression_in_infiltrated_tissues.pdf",
  plot = p_gene_age_trends,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600
)

# ============================================================
# PART E. Pearson correlation heatmap across selected tissues
# ============================================================

cor_stats_heatmap <- cor_stats_genes %>%
  mutate(
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(round(cor, 2), sig)
  ) %>%
  filter(Tissue %in% c("Marrow", "Spleen", "GAT", "Pancreas")) %>%
  mutate(Tissue = factor(Tissue, levels = c("Pancreas", "GAT", "Spleen", "Marrow")))

p_pearson_heatmap <- ggplot(cor_stats_heatmap, aes(x = Gene, y = Tissue, fill = cor)) +
  geom_tile(color = "black") +
  geom_text(aes(label = label), size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson (R)"
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 17) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  coord_fixed(ratio = 1)

print(p_pearson_heatmap)

ggsave(
  filename = "Pearson_Il10_4_Tissues.pdf",
  plot = p_pearson_heatmap,
  width = 4,
  height = 4,
  units = "in",
  dpi = 600
)



# Respective to Figure S1c
# ------------------------------------------------------------
# Gene Set Enrichment Analysis 
# ------------------------------------------------------------

combined_data <- bind_rows(
  DEGs_CD4_pMT10_plus_vs_minus_Zn,
  DEGs_CD4_pMT10IFNgKO_plus_vs_minus_Zn,
  DEGs_CD8_pMT10_plus_vs_minus_Zn,
  DEGs_CD8_pMT10IFNgKO_plus_vs_minus_Zn)


library(msigdbr)
c2cp <- read.gmt("h.all.v2023.1.Hs.entrez.gmt.txt")
hs_gsea <- msigdbr(species = "Mus musculus") 
hs_gsea %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat)

hs_gsea_c2 <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# GSEA wrapper function
run_gsea <- function(dataset_name, hs_gsea_c2) {
  df_sub <- combined_data %>%
    dplyr::filter(Dataset == dataset_name) %>%
    dplyr::select(gene_name, log2FC) %>%
    dplyr::distinct()
  
  # Create ranked gene list
  gsea_data <- df_sub$log2FC
  names(gsea_data) <- as.character(df_sub$gene_name)
  gsea_data <- sort(gsea_data, decreasing = TRUE)
  
  # Run GSEA with error handling
  gsea_res <- tryCatch({
    GSEA(gsea_data, TERM2GENE = hs_gsea_c2, verbose = FALSE, eps = 0)
  }, error = function(e) {
    message(paste("GSEA failed for", dataset_name, ":", e$message))
    return(NULL)
  })
  
  # Convert to tibble if result is valid
  if (!is.null(gsea_res) && "result" %in% slotNames(gsea_res)) {
    gsea_tbl <- as_tibble(gsea_res@result)
    
    if (nrow(gsea_tbl) > 0) {
      gsea_tbl <- gsea_tbl %>%
        mutate(phenotype = dataset_name, ID = sub("^HALLMARK_", "", ID))
      return(gsea_tbl)
    }
  }
  return(NULL)
}

# Run GSEA for each dataset
dataset_list <- unique(combined_data$Dataset)

gsea_results <- lapply(dataset_list, run_gsea, hs_gsea_c2 = hs_gsea_c2)
gsea_results <- Filter(Negate(is.null), gsea_results)

# Combine results
if (length(gsea_results) == 0) stop("No GSEA results were returned.")
gsea_combined <- bind_rows(gsea_results)

gsea_combined <- gsea_combined %>%
  dplyr::arrange(desc(NES)) %>%
  mutate(phenotype = factor(phenotype, levels = dataset_list))

top_GSEA <- gsea_combined %>%
  group_by(phenotype) %>%
  slice_max(order_by = NES, n = 50) %>%
  bind_rows(
    gsea_combined %>%
      group_by(phenotype) %>%
      slice_min(order_by = NES, n = 50)
  ) %>%
  ungroup()

data_combined <- ggplot(top_GSEA, aes(x = phenotype, y = ID)) +
  geom_point(aes(size = -log10(p.adjust), fill = NES), shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", size = 1),
    axis.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "right"
  ) +
  ggtitle("GSEA Results Il10 Sig for 4 DEG-Conditions")

print(data_combined)





# Respective to Figure S1f
# ============================================================
# Transcriptome directionality: CD4 vs CD8
# ============================================================

library(dplyr)
library(ggplot2)
library(ggrepel)

degs_cd8 <- DEGs_CD8_pMT10_plus_vs_minus_Zn
degs_cd4 <- DEGs_CD4_pMT10_plus_vs_minus_Zn

degs_cd8_selected <- degs_cd8 %>%
  dplyr::select(gene_name, log2FC, pvalAdj) %>%
  mutate(Dataset = "CD8_pMT10_plus_vs_minus_Zn")

degs_cd4_selected <- degs_cd4 %>%
  dplyr::select(gene_name, log2FC, pvalAdj) %>%
  mutate(Dataset = "CD4_pMT10_plus_vs_minus_Zn")

combined_directionality_data <- bind_rows(
  degs_cd8_selected,
  degs_cd4_selected
)

cd8_data <- combined_directionality_data %>%
  filter(Dataset == "CD8_pMT10_plus_vs_minus_Zn")

cd4_data <- combined_directionality_data %>%
  filter(Dataset == "CD4_pMT10_plus_vs_minus_Zn")

merged_directionality_data <- merge(
  cd8_data,
  cd4_data,
  by = "gene_name",
  suffixes = c("_CD8", "_CD4")
)

merged_directionality_data <- merged_directionality_data %>%
  mutate(
    same_direction = ifelse(
      (log2FC_CD8 > 0 & log2FC_CD4 > 0) |
        (log2FC_CD8 < 0 & log2FC_CD4 < 0),
      "Same Direction",
      "Different Direction"
    )
  )

same_direction_genes <- merged_directionality_data %>%
  filter(same_direction == "Same Direction")

top_labeled_genes <- merged_directionality_data %>%
  filter(same_direction == "Same Direction") %>%
  arrange(desc(log2FC_CD8)) %>%
  head(25) %>%
  bind_rows(
    merged_directionality_data %>%
      filter(same_direction == "Same Direction") %>%
      arrange(desc(log2FC_CD4)) %>%
      head(25)
  ) %>%
  distinct(gene_name, .keep_all = TRUE)

# Optional check for a specific gene
merged_directionality_data[merged_directionality_data$gene_name == "Il10rb", ]

p_transcriptome_directionality <- ggplot(
  merged_directionality_data,
  aes(x = log2FC_CD8, y = log2FC_CD4, color = same_direction)
) +
  geom_point(alpha = 0.6) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black",
    size = 1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    size = 1
  ) +
  geom_text_repel(
    data = top_labeled_genes,
    aes(label = gene_name),
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.5,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Same Direction" = "blue",
      "Different Direction" = "orange"
    )
  ) +
  labs(
    title = "Transcriptome Directionality in TCD4 and TCD8",
    x = "log2FC (TCD8)",
    y = "log2FC (TCD4)",
    color = "Direction"
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", size = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

print(p_transcriptome_directionality)

ggsave(
  filename = "Transcriptome_Directionality_TCD4_TCD8.pdf",
  plot = p_transcriptome_directionality,
  device = pdf,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

