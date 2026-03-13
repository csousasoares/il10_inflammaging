## Load Packages and Set Seed --------------------------------------------------

library(Seurat)
library(sctransform)
library(tidyverse)
library(patchwork)
library(sctransform)
library(presto)
library(org.Mm.eg.db)
library(ggrepel)
library(UCell)
library(org.Mm.eg.db)
library(stringr)

set.seed(123)


## gWAT Analysis ---------------------------------------------------------------

gat <- readRDS(file = "seurat_objects/gonadal_fat_pad.rds") ## Load object

DimPlot(gat, reduction = "umap", label = T, repel = T)

gat$cell_type <- stringr::str_wrap(gat$cell_type, 
                                   width = 25) ## Define before setting Idents

Idents(gat) <- "cell_type" 

## Set cell_type as identity for DimPlot, subset, etc.

gat_dimplot <- DimPlot(gat, 
                       reduction = "umap", 
                       label = F, 
                       repel = T, 
                       pt.size = 0.2,
                       alpha = 0.5) + 
  guides(color = guide_legend(override.aes = list(size=2), ncol=1)) +
  theme(legend.position = "left",
        legend.text = element_text(size = 10)) ## Create DimPlot

gat_dimplot


gat_dimplot_raster <- gat_dimplot %>% 
  ggrastr::rasterise(layers = "Point",
                     dpi = 600) ## Rasterize points

gat_dimplot_raster

il_10 <- c("Ccl3", "Ccl4", "Ccl5", "Il10", "Gzmk", "Gzmb", "Cd38", 
                    "Eomes", "Lag3", "Tox2", "Pdcd1", "Tigit", "Prf1", 
                    "Cdkn1a", "Cdkn2a", "Ccna2") ## Get IL-10 Signature


annots <- AnnotationDbi::select(org.Mm.eg.db, keys=il_10, 
                                columns="ENSEMBL", 
                                keytype="SYMBOL") ## Get Ensembl IDs

il_10_ensembl <- annots$ENSEMBL
il_10_ensembl

signatures <- list(
  il_10_sig = il_10_ensembl) ## Load IL-10 Sig. to UCell

gat <- AddModuleScore_UCell(gat, features=signatures, 
                            name=NULL) ## Calculate IL-10 Sig. Score

gat_feature_plot <- FeaturePlot(gat, 
                                reduction = "umap", 
                                features = "il_10_sig", 
                                split.by = "age",
                                pt.size = 0.3,
                                combine = T,
                                repel = T,
                                alpha = 0.65,
                                min.cutoff = 0, 
                                max.cutoff = 0.2,
                                order = T) +
  theme(legend.position = "right") &
  labs(x = "UMAP1", y = "UMAP2") &
  scale_color_gradient(low = "azure2", high = "red")

gat_feature_plot_1 <- gat_feature_plot[[1]] %>% 
  ggrastr::rasterise(layers = "Point", dpi = 600)

gat_feature_plot_2 <- gat_feature_plot[[2]] %>% 
  ggrastr::rasterise(layers = "Point", dpi = 600)

gat_feature_plot_3 <- gat_feature_plot[[3]] %>% 
  ggrastr::rasterise(layers = "Point", dpi = 600)
gat_feature_plot_3


wrap_plots(gat_dimplot_raster, gat_feature_plot_1, gat_feature_plot_2, 
           gat_feature_plot_3, widths = c(1,1,1,1))


ggsave(
  "gWAT_umap_il10_sig.pdf",
  path = "Plots",
  height = 5,
  width = 17,
  unit = "in"
)


## gWAT Analysis Custom Annotations --------------------------------------------


gat$cell_type <- stringr::str_wrap(gat$cell_type, 
                                   width = 600) ## Define before setting Idents


unique(gat$cell_type)

class(gat@meta.data$cell_type) ## Is a factor, should be character to edit

gat@meta.data$cell_type2 <- as.character(gat@meta.data$cell_type) 

## Apparently the vector has to be a character, not a factor

t_cells <- c("CD8-positive, alpha-beta T cell", 
             "CD4-positive, alpha-beta T cell", 
             "T cell")

gat@meta.data$cell_type2[gat@meta.data$cell_type %in% t_cells] <- "T cell"
gat@meta.data$cell_type2[gat@meta.data$cell_type == "mesenchymal stem cell of adipose tissue"] <- "mesenchymal stem cell"

Idents(gat) <- "cell_type2"

gat_dimplot <- DimPlot(gat, 
                       reduction = "umap", 
                       label = F, 
                       repel = T, 
                       pt.size = 0.2,
                       alpha = 0.5) + 
  guides(color = guide_legend(override.aes = list(size=2), ncol=1)) +
  theme(legend.position = "left",
        legend.text = element_text(size = 10)) ## Create DimPlot

gat_dimplot

gat_dimplot_raster <- gat_dimplot %>% 
  ggrastr::rasterise(layers = "Point",
                     dpi = 600) ## Rasterize points

gat_dimplot_raster


wrap_plots(gat_dimplot_raster, gat_feature_plot_1, gat_feature_plot_2, 
           gat_feature_plot_3, widths = c(1,1,1,1))


ggsave(
  "gWAT_umap_il10_sig_new_annotations.pdf",
  path = "Plots",
  height = 5,
  width = 17,
  unit = "in"
)


## Bone Marrow Analysis --------------------------------------------------------

bm <- readRDS(file = "seurat_objects/bm.rds") ## Load object

DimPlot(bm, reduction = "umap", label = T, repel = T) ## Check data

bm$cell_type <- stringr::str_wrap(bm$cell_type, 
                                   width = 25) ## Define before setting Idents

Idents(bm) <- "cell_type" 

## Set cell_type as identity for DimPlot, subset, etc.


bm_dimplot <- DimPlot(bm, 
                       reduction = "umap", 
                       label = F, 
                       repel = T, 
                       pt.size = 0.2,
                       alpha = 0.5) + 
  guides(color = guide_legend(override.aes = list(size=2), ncol=1)) +
  theme(legend.position = "left",
        legend.text = element_text(size = 10)) ## Create DimPlot

bm_dimplot


bm_dimplot_raster <- bm_dimplot %>% 
  ggrastr::rasterise(layers = "Point",
                     dpi = 600) ## Rasterize points

bm_dimplot_raster


il_10 <- c("Ccl3", "Ccl4", "Ccl5", "Il10", "Gzmk", "Gzmb", "Cd38", 
           "Eomes", "Lag3", "Tox2", "Pdcd1", "Tigit", "Prf1", 
           "Cdkn1a", "Cdkn2a", "Ccna2") ## Get IL-10 Signature


annots <- AnnotationDbi::select(org.Mm.eg.db, keys=il_10, 
                                columns="ENSEMBL", 
                                keytype="SYMBOL") ## Get Ensembl IDs

il_10_ensembl <- annots$ENSEMBL
il_10_ensembl

signatures <- list(
  il_10_sig = il_10_ensembl) ## Load IL-10 Sig. to UCell

bm <- AddModuleScore_UCell(bm, features=signatures, 
                            name=NULL) ## Calculate IL-10 Sig. Score

bm_feature_plot <- FeaturePlot(bm, 
                                reduction = "umap", 
                                features = "il_10_sig", 
                                split.by = "age",
                                pt.size = 0.3,
                                combine = T,
                                repel = T,
                                alpha = 0.65,
                                min.cutoff = 0, 
                                max.cutoff = 0.22,
                                order = T) +
  theme(legend.position = "right") &
  labs(x = "UMAP1", y = "UMAP2") &
  scale_color_gradient(low = "azure2", high = "red")

bm_feature_plot_1 <- bm_feature_plot[[1]] %>% 
  ggrastr::rasterise(layers = "Point", dpi = 600)

bm_feature_plot_2 <- bm_feature_plot[[2]] %>% 
  ggrastr::rasterise(layers = "Point", dpi = 600)

bm_feature_plot_3 <- bm_feature_plot[[3]] %>% 
  ggrastr::rasterise(layers = "Point", dpi = 600)
bm_feature_plot_3


wrap_plots(bm_dimplot_raster, bm_feature_plot_1, bm_feature_plot_2, 
           bm_feature_plot_3, widths = c(1,1,1,1))


ggsave(
  "bone_marrow_umap_il10_sig.pdf",
  path = "Plots",
  height = 5,
  width = 17,
  unit = "in"
)


## Spleen Analysis -------------------------------------------------------------


spleen <- readRDS(file = "seurat_objects/spleen.rds") ## Load object

DimPlot(spleen, reduction = "umap", label = T, repel = T) ## Check data

spleen$cell_type <- stringr::str_wrap(spleen$cell_type, 
                                  width = 25) ## Define before setting Idents

Idents(spleen) <- "cell_type" 

## Set cell_type as identity for DimPlot, subset, etc.



spleen_dimplot <- DimPlot(spleen, 
                      reduction = "umap", 
                      label = F, 
                      repel = T, 
                      pt.size = 0.3,
                      alpha = 0.5) + 
  guides(color = guide_legend(override.aes = list(size=2), ncol=1)) +
  theme(legend.position = "left",
        legend.text = element_text(size = 10)) ## Create DimPlot

spleen_dimplot


spleen_dimplot_raster <- spleen_dimplot %>% 
  ggrastr::rasterise(layers = "Point",
                     dpi = 600) ## Rasterize points

spleen_dimplot_raster


il_10 <- c("Ccl3", "Ccl4", "Ccl5", "Il10", "Gzmk", "Gzmb", "Cd38", 
           "Eomes", "Lag3", "Tox2", "Pdcd1", "Tigit", "Prf1", 
           "Cdkn1a", "Cdkn2a", "Ccna2") ## Get IL-10 Signature


annots <- AnnotationDbi::select(org.Mm.eg.db, keys=il_10, 
                                columns="ENSEMBL", 
                                keytype="SYMBOL") ## Get Ensembl IDs

il_10_ensembl <- annots$ENSEMBL
il_10_ensembl

signatures <- list(
  il_10_sig = il_10_ensembl) ## Load IL-10 Sig. to UCell

spleen <- AddModuleScore_UCell(spleen, features=signatures, 
                           name=NULL) ## Calculate IL-10 Sig. Score

spleen_feature_plot <- FeaturePlot(spleen, 
                               reduction = "umap", 
                               features = "il_10_sig", 
                               split.by = "age",
                               pt.size = 0.4,
                               combine = T,
                               repel = T,
                               alpha = 0.65,
                               min.cutoff = 0, 
                               max.cutoff = 0.25,
                               order = T) +
  theme(legend.position = "right") &
  labs(x = "UMAP1", y = "UMAP2") &
  scale_color_gradient(low = "azure2", high = "red")

spleen_feature_plot_1 <- spleen_feature_plot[[1]] %>% 
  ggrastr::rasterise(layers = "Point", dpi = 600)

spleen_feature_plot_2 <- spleen_feature_plot[[2]] %>% 
  ggrastr::rasterise(layers = "Point", dpi = 600)

spleen_feature_plot_3 <- spleen_feature_plot[[3]] %>% 
  ggrastr::rasterise(layers = "Point", dpi = 600)
spleen_feature_plot_3


wrap_plots(spleen_dimplot_raster, spleen_feature_plot_1, spleen_feature_plot_2, 
           spleen_feature_plot_3, widths = c(1,1,1,1))


ggsave(
  "spleen_umap_il10_sig.pdf",
  path = "Plots",
  height = 5,
  width = 17,
  unit = "in"
)



## Spleen Analysis Custom Annotations ------------------------------------------


spleen$cell_type <- stringr::str_wrap(spleen$cell_type, 
                                   width = 800) ## Define before setting Idents


unique(spleen$cell_type)

class(spleen@meta.data$cell_type) ## Is a factor, should be character to edit

spleen@meta.data$cell_type2 <- as.character(spleen@meta.data$cell_type) 

## Apparently the vector has to be a character, not a factor

t_cells <- c("CD8-positive, alpha-beta T cell", 
             "CD4-positive, alpha-beta T cell", 
             "T cell")

spleen@meta.data$cell_type2[spleen@meta.data$cell_type %in% t_cells] <- "T cell"

Idents(spleen) <- "cell_type2"

spleen_dimplot <- DimPlot(spleen, 
                       reduction = "umap", 
                       label = F, 
                       repel = T, 
                       pt.size = 0.2,
                       alpha = 0.5) + 
  guides(color = guide_legend(override.aes = list(size=2), ncol=1)) +
  theme(legend.position = "left",
        legend.text = element_text(size = 10)) ## Create DimPlot

spleen_dimplot

spleen_dimplot_raster <- spleen_dimplot %>% 
  ggrastr::rasterise(layers = "Point",
                     dpi = 600) ## Rasterize points

spleen_dimplot_raster


wrap_plots(spleen_dimplot_raster, spleen_feature_plot_1, spleen_feature_plot_2, 
           spleen_feature_plot_3, widths = c(1,1,1,1))


ggsave(
  "spleen_umap_il10_sig_new_annotations.pdf",
  path = "Plots",
  height = 5,
  width = 17,
  unit = "in"
)


## Pancreas FACS/Smart-Seq2 Analysis -------------------------------------------


pancreas <- readRDS(file = "seurat_objects/pancreas.rds") ## Load object

DimPlot(pancreas, reduction = "umap", label = T, repel = T) ## Check data

pancreas$cell_type <- stringr::str_wrap(pancreas$cell_type, 
                                      width = 25) ## Define before setting Idents

Idents(pancreas) <- "cell_type" 

## Set cell_type as identity for DimPlot, subset, etc.


pancreas_dimplot <- DimPlot(pancreas, 
                          reduction = "umap", 
                          label = F, 
                          repel = T, 
                          pt.size = 0.3,
                          alpha = 0.5) + 
  guides(color = guide_legend(override.aes = list(size=2), ncol=1)) +
  theme(legend.position = "left",
        legend.text = element_text(size = 10)) ## Create DimPlot

pancreas_dimplot


pancreas_dimplot_raster <- pancreas_dimplot %>% 
  ggrastr::rasterise(layers = "Point",
                     dpi = 600) ## Rasterize points

pancreas_dimplot_raster


il_10 <- c("Ccl3", "Ccl4", "Ccl5", "Il10", "Gzmk", "Gzmb", "Cd38", 
           "Eomes", "Lag3", "Tox2", "Pdcd1", "Tigit", "Prf1", 
           "Cdkn1a", "Cdkn2a", "Ccna2") ## Get IL-10 Signature


annots <- AnnotationDbi::select(org.Mm.eg.db, keys=il_10, 
                                columns="ENSEMBL", 
                                keytype="SYMBOL") ## Get Ensembl IDs

il_10_ensembl <- annots$ENSEMBL
il_10_ensembl

signatures <- list(
  il_10_sig = il_10_ensembl) ## Load IL-10 Sig. to UCell

pancreas <- AddModuleScore_UCell(pancreas, features=signatures, 
                               name=NULL) ## Calculate IL-10 Sig. Score

pancreas_feature_plot <- FeaturePlot(pancreas, 
                                   reduction = "umap", 
                                   features = "il_10_sig", 
                                   split.by = "age",
                                   pt.size = 0.4,
                                   combine = T,
                                   repel = T,
                                   alpha = 0.65,
                                   min.cutoff = 0, 
                                   max.cutoff = 0.2,
                                   order = T) +
  theme(legend.position = "right") &
  labs(x = "UMAP1", y = "UMAP2") &
  scale_color_gradient(low = "azure2", high = "red")

pancreas_feature_plot_1 <- pancreas_feature_plot[[1]] %>% 
  ggrastr::rasterise(layers = "Point", dpi = 600)

pancreas_feature_plot_2 <- pancreas_feature_plot[[2]] %>% 
  ggrastr::rasterise(layers = "Point", dpi = 600)

pancreas_feature_plot_3 <- pancreas_feature_plot[[3]] %>% 
  ggrastr::rasterise(layers = "Point", dpi = 600)
pancreas_feature_plot_3


wrap_plots(pancreas_dimplot_raster, pancreas_feature_plot_1, pancreas_feature_plot_2, 
           pancreas_feature_plot_3, widths = c(1,1,1,1))


ggsave(
  "pancreas_umap_il10_sig.pdf",
  path = "Plots",
  height = 5,
  width = 17,
  unit = "in"
)

writeLines(capture.output(sessionInfo()), "sessionInfo_TMS_UCell.txt")
