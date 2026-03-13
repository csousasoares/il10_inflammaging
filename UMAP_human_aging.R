## Load Packages and Set Seed --------------------------------------------------

set.seed(123)

library(Seurat)
library(ggplot2)
library(org.Hs.eg.db)
library(UCell)
library(patchwork)
library(scales)
library(ggsc)
library(Cairo)
library(stringr)

## Load and Inspect Blood T-Cell Seurat Object ---------------------------------

blood <- readRDS("input_data\\blood_human.rds") ## Load T-cells in the Blood.

Idents(blood) <- "cell_type"

unique(blood@meta.data[["MMoCHi_annotation"]])

mochi_all_dimplot <- DimPlot(blood, 
                             reduction = "umap", 
                             group.by = "MMoCHi_annotation",
                             raster = F,
                             alpha = 0.5,
                             raster.dpi = c(2048,2048),
                             pt.size = 0.01) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 

mochi_all_dimplot


mochi_all_dimplot <- ggrastr::rasterise(mochi_all_dimplot,
                                        layers='Point',
                                        dpi = 600)
mochi_all_dimplot

CD8_tcells <- c("cd8_temra", "cd8_naive", "cd8_tem", 
                "cd8_trm", "cd8_cm", "cd8_mait")
CD4_tcells <- c("cd4_cm", "cd4_naive", "cd4_tem", 
                "cd4_treg", "cd4_temra", "cd4_trm")
Gamma_Delta <- c("cd3_gd")

blood@meta.data$cell_type_broad[blood@meta.data$MMoCHi_annotation %in% CD4_tcells] <- "CD4+ T Cell"
blood@meta.data$cell_type_broad[blood@meta.data$MMoCHi_annotation %in% CD8_tcells] <- "CD8+ T Cell"
blood@meta.data$cell_type_broad[blood@meta.data$MMoCHi_annotation %in% Gamma_Delta] <- "Gamma-Delta T Cell"

unique(blood@meta.data$cell_type_broad)

mochi_cd4_8_gd_dimplot <- DimPlot(blood, 
                                  reduction = "umap", 
                                  group.by = "cell_type_broad",
                                  raster = F,
                                  alpha = 0.5,
                                  pt.size = 0.01,
                                  raster.dpi = c(2048,2048)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 

mochi_cd4_8_gd_dimplot

mochi_cd4_8_gd_dimplot <- ggrastr::rasterise(mochi_cd4_8_gd_dimplot,
                                             layers='Point',
                                             dpi = 600)
mochi_cd4_8_gd_dimplot


genes <- c("CD8A", "CD4")

annots <- AnnotationDbi::select(org.Hs.eg.db, keys=genes, 
                                columns="ENSEMBL", keytype="SYMBOL")

annots$ENSEMBL
intersect(row.names(blood), annots$ENSEMBL)


cd8_featureplot <- FeaturePlot(blood, ## Confirm CD4 and CD8 identity
                               features = c("ENSG00000153563"),
                               reduction = "umap",
                               pt.size = 0.01,
                               combine = T,
                               alpha = 0.5,
                               repel = T,
                               order = T,
                               raster = F,
                               raster.dpi = c(2048, 2048)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") &
  scale_color_viridis_b(option = "C", limits = c(0, 1), 
                        n.breaks = 8, direction = 1)


cd8_featureplot <- ggrastr::rasterise(cd8_featureplot,
                                      layers='Point',
                                      dpi = 600)
cd8_featureplot

cd4_featureplot <- FeaturePlot(blood, ## Confirm CD4 and CD8 identity
                               features = c("ENSG00000010610"),
                               reduction = "umap",
                               pt.size = 0.01,
                               combine = T,
                               alpha = 0.6,
                               raster = F,
                               raster.dpi = c(2048, 2048),
                               order = T) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") &
  scale_color_viridis_b(option = "C", limits = c(0, 0.3), 
                        n.breaks = 8, direction = 1)

cd4_featureplot <- ggrastr::rasterise(cd4_featureplot,
                                      layers='Point',
                                      dpi = 600)
cd4_featureplot


wrap_plots(mochi_cd4_8_gd_dimplot, 
           mochi_all_dimplot, 
           cd8_featureplot, 
           cd4_featureplot,
           widths = c(1,1,1,1))


ggsave(
  "T-Cells in Blood.pdf",
  path = "output_data\\umap_human_aging",
  width = 20,
  height = 6,
  units = "in"
)

## Blood T-Cell Aging Analysis -------------------------------------------------

## Divide ages into two groups as in the original paper:

unique(blood@meta.data$donor_age)

old <- c("56 years",
         "67 years",
         "70-74 years",
         "68 years",
         "64 years",
         "50-54 years",
         "40-44 years",
         "48 years",
         "55-60 years",
         "75 years",
         "73 years") ## Over 40y

young <- c("20-24 years",
           "23 years",
           "25 years",
           "35-39 years",
           "33 years",
           "27 years",
           "25-30 years",
           "20 years") ## Under 40y


blood@meta.data$age_group[blood@meta.data$donor_age %in% old] <- "older"
blood@meta.data$age_group[blood@meta.data$donor_age %in% young] <- "younger"
blood@meta.data$age_group <- factor(blood@meta.data$age_group, levels = c("younger", "older"))

## Now this will allow to work with <40y vs >40y cohort.

DimPlot(blood, 
        reduction = "umap", 
        group.by = "age_group",
        raster = T,
        pt.size = 4,
        raster.dpi = c(1024,1024)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 


## Using the annotation from the original paper:

Idents(blood) <- "MMoCHi_annotation"
Idents(blood)

## CD8 T-cells:

CD8_cells_mochi <- subset(blood, idents = CD8_tcells)
CD8_cells_mochi

# Wrap cell type names at ~15 characters (adjust width to your preference):

CD8_cells_mochi$MMoCHi_annotation <- stringr::str_wrap(
  CD8_cells_mochi$MMoCHi_annotation, width = 15)


cd8_dimplot <- DimPlot(CD8_cells_mochi, 
                       reduction = "umap", 
                       group.by = "MMoCHi_annotation",
                       raster = F,
                       alpha = 0.5,
                       raster.dpi = c(2048,2048),
                       pt.size = 0.01) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "left") & xlim(-10, 20) & 
  ylim(-10, 20)

cd8_dimplot <- ggrastr::rasterise(cd8_dimplot, 
                                  layers='Point'
                                  ,dpi = 600)
cd8_dimplot

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

annots <- AnnotationDbi::select(org.Hs.eg.db, 
                                keys = il_10_human, 
                                columns="ENSEMBL", 
                                keytype="SYMBOL") ## Convert genes.

il_10_ensembl <- annots$ENSEMBL
il_10_ensembl
il_10_ensembl <- intersect(il_10_ensembl, row.names(blood))

signatures <- list(
  il_10_sig = il_10_ensembl)

CD8_cells_mochi <- AddModuleScore_UCell(CD8_cells_mochi, 
                                        features=signatures, name=NULL)

cd8_feature_plot <- FeaturePlot(CD8_cells_mochi, 
                                features = c("il_10_sig"),
                                reduction = "umap",
                                pt.size = 0.01,
                                alpha = 0.5,
                                cols = c("azure2", "red"),
                                split.by = "age_group",
                                combine = T,
                                repel = T,
                                keep.scale = "feature",
                                order = T,
                                raster = F,
                                raster.dpi = c(512, 512)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") & 
  scale_color_viridis_b(option = "B",
                        limits = c(0,0.3),
                        n.breaks = 10) & xlim(-10, 20) & ylim(-10, 20)

cd8_feature_plot_1 <- ggrastr::rasterise(cd8_feature_plot[[1]], 
                                         layers = "Point", 
                                         dpi = 600)
cd8_feature_plot_1


cd8_feature_plot_2 <- ggrastr::rasterise(cd8_feature_plot[[2]], 
                                         layers = "Point", 
                                         dpi = 600)
cd8_feature_plot_2


wrap_plots(cd8_dimplot, 
           cd8_feature_plot_1,
           cd8_feature_plot_2,
           widths = c(1, 1, 1)) 

ggsave(
  "CD8 T-Cells Blood IL-10 Signature MMoCHi Annot.pdf",
  path = "output_data\\umap_human_aging/",
  width = 15,
  height = 5.25,
  units = "in"
)

## CD4 T-cells:

CD4_cells_mochi <- subset(blood, idents = CD4_tcells)

CD4_cells_mochi$MMoCHi_annotation <- stringr::str_wrap(
  CD4_cells_mochi$MMoCHi_annotation, width = 15)


cd4_dimplot <- DimPlot(CD4_cells_mochi, 
                       reduction = "umap", 
                       group.by = "MMoCHi_annotation",
                       raster = F,
                       alpha = 0.5,
                       raster.dpi = c(2048,2048),
                       pt.size = 0.01) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "left") & xlim(-10, 20) & 
  ylim(-10, 20)

cd4_dimplot <- ggrastr::rasterise(cd4_dimplot, 
                                  layers = "Point", 
                                  dpi = 600)

signatures <- list(
  il_10_sig = il_10_ensembl)

CD4_cells_mochi <- AddModuleScore_UCell(CD4_cells_mochi, 
                                        features = signatures, 
                                        name = NULL)

cd4_feature_plot <- FeaturePlot(CD4_cells_mochi, 
                                features = c("il_10_sig"),
                                reduction = "umap",
                                pt.size = 0.01,
                                alpha = 0.5,
                                cols = c("azure2", "red"),
                                split.by = "age_group",
                                combine = T,
                                repel = T,
                                keep.scale = "feature",
                                order = T,
                                raster = F,
                                raster.dpi = c(512, 512)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") & 
  scale_color_viridis_b(option = "B",
                        limits = c(0,0.25),
                        n.breaks = 10) & xlim(-10, 20) & 
  ylim(-10, 20)

cd4_feature_plot_1 <- ggrastr::rasterise(cd4_feature_plot[[1]], 
                                         layers = "Point", 
                                         dpi = 600)
cd4_feature_plot_1

cd4_feature_plot_2 <- ggrastr::rasterise(cd4_feature_plot[[2]], 
                                         layers = "Point", 
                                         dpi = 600)
cd4_feature_plot_2


wrap_plots(cd4_dimplot, 
           cd4_feature_plot_1, 
           cd4_feature_plot_2,
           widths = c(1, 1, 1))

ggsave(
  "CD4 T-Cells Blood IL-10 Signature MMoCHi Annot.pdf",
  path = "output_data\\umap_human_aging/",
  width = 15,
  height = 5.25,
  units = "in"
)

## Blood w/ New Color Scheme ---------------------------------------------------

## All T cells in the blood:

cols <- c(
  "cd3_gd" = "gray",
  "cd8_trm" = "orange3",
  "cd4_cm" = "deepskyblue2",
  "cd8_naive" = "deeppink",
  "cd4_tem" = "aquamarine",
  "cd4_trm" = "skyblue4",
  "cd8_tem" = "salmon",
  "cd4_treg" = "cornflowerblue",
  "cd8_temra" = "red",
  "cd4_naive" = "cyan2",
  "cd8_mait" = "orange",
  "cd4_temra" = "slateblue2",
  "cd8_cm" = "yellow3"
)


all_dimplot <- DimPlot(blood, 
                       reduction = "umap", 
                       group.by = "MMoCHi_annotation",
                       raster = F,
                       label = T,
                       label.size = 3,
                       alpha = 0.4,
                       pt.size = 0.01,
                       raster.dpi = c(2048,2038)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8, lineheight = 1),
        legend.position = "left",
        legend.key.spacing.y = unit(6, "pt")) &
  xlim(c(-10,25)) & ylim(c(-10, 20)) & 
  scale_color_manual(values = cols)

all_dimplot <- ggrastr::rasterise(all_dimplot, 
                                  layers = "Point", 
                                  dpi = 600)
all_dimplot

signatures <- list(
  il_10_sig = il_10_ensembl)

blood <- AddModuleScore_UCell(blood, 
                              features = signatures, 
                              name = NULL)

blood_feature_plot <- FeaturePlot(blood, 
                                  features = c("il_10_sig"),
                                  reduction = "umap",
                                  pt.size = 0.01,
                                  alpha = 0.4,
                                  split.by = "age_group",
                                  combine = T,
                                  repel = T,
                                  keep.scale = "feature",
                                  order = T,
                                  raster = F,
                                  raster.dpi = c(512, 512)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") & 
  scale_color_viridis_b(option = "B",
                        limits = c(0,0.3),
                        n.breaks = 10) & xlim(c(-10,25)) & 
  ylim(c(-10, 20))

blood_feature_plot

blood_feature_plot_1 <- ggrastr::rasterise(blood_feature_plot[[1]], 
                                           layers = "Point", 
                                           dpi = 600)
blood_feature_plot_1

blood_feature_plot_2 <- ggrastr::rasterise(blood_feature_plot[[2]], 
                                           layers = "Point", 
                                           dpi = 600)
blood_feature_plot_2


blood_umap <- wrap_plots(all_dimplot, 
                         blood_feature_plot_1,
                         blood_feature_plot_2,
                         widths = c(1, 1, 1)) 

blood_umap

saveRDS(blood_umap, "blood_umap.rds") 

ggsave(
  "All T-Cells Blood IL-10 Signature MMoCHi Annot.pdf",
  path = "output_data\\umap_human_aging",
  width = 20,
  height = 7,
  units = "in"
)


cols <- c(
  "cd3_gd" = "gray",
  "cd8_trm" = "orange3",
  "cd4_cm" = "deepskyblue2",
  "cd8_naive" = "deeppink",
  "cd4_tem" = "aquamarine",
  "cd4_trm" = "skyblue4",
  "cd8_tem" = "salmon",
  "cd4_treg" = "cornflowerblue",
  "cd8_temra" = "red",
  "cd4_naive" = "cyan2",
  "cd8_mait" = "orange",
  "cd4_temra" = "slateblue2",
  "cd8_cm" = "yellow3"
)


Idents(blood) <- "MMoCHi_annotation"
Idents(blood)

CD8_cells_mochi <- subset(blood, idents = CD8_tcells)
CD8_cells_mochi

CD8_cells_mochi$MMoCHi_annotation <- stringr::str_wrap(
  CD8_cells_mochi$MMoCHi_annotation, 
  width = 15)


cd8_dimplot <- DimPlot(CD8_cells_mochi, 
                       reduction = "umap", 
                       group.by = "MMoCHi_annotation",
                       raster = F,
                       pt.size = 0.05,
                       label = T,
                       alpha = 0.7,
                       raster.dpi = c(2048,2038)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8, lineheight = 1),
        legend.position = "left",
        legend.key.spacing.y = unit(6, "pt")) + 
  scale_color_manual(values = cols)

cd8_dimplot


signatures <- list(
  il_10_sig = il_10_ensembl)

CD8_cells_mochi <- AddModuleScore_UCell(CD8_cells_mochi, 
                                        features = signatures, name=NULL)

cd8_feature_plot <- FeaturePlot(CD8_cells_mochi, 
                                features = c("il_10_sig"),
                                reduction = "umap",
                                pt.size = 0.05,
                                alpha = 0.7,
                                cols = c("azure2", "red"),
                                split.by = "age_group",
                                combine = T,
                                repel = T,
                                keep.scale = "feature",
                                order = T,
                                raster = F,
                                raster.dpi = c(512, 512)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") & 
  scale_color_viridis_b(option = "B",
                        limits = c(0,0.3),
                        n.breaks = 10)


cd8_dimplot <- ggrastr::rasterise(cd8_dimplot[[1]], 
                                  layers = "Point", dpi = 600)
cd8_dimplot

cd8_feature_plot_1 <- ggrastr::rasterise(cd8_feature_plot[[1]], 
                                         layers = "Point", 
                                         dpi = 600)
cd8_feature_plot_1


cd8_feature_plot_2 <- ggrastr::rasterise(cd8_feature_plot[[2]], 
                                         layers = "Point", 
                                         dpi = 600)
cd8_feature_plot_2

wrap_plots(cd8_dimplot,
           cd8_feature_plot_1,
           cd8_feature_plot_2,
           widths = c(1, 1, 1))

ggsave(
  "CD8 T-Cells Blood IL-10 Signature MMoCHi Annot.pdf",
  path = "output_data\\umap_human_aging/",
  width = 15,
  height = 4,
  units = "in"
)

CD4_cells_mochi <- subset(blood, idents = CD4_tcells)

CD4_cells_mochi$MMoCHi_annotation <- stringr::str_wrap(
  CD4_cells_mochi$MMoCHi_annotation, width = 15)


cd4_dimplot <- DimPlot(CD4_cells_mochi, 
                       reduction = "umap", 
                       group.by = "MMoCHi_annotation",
                       raster = F,
                       label = T,
                       label.size = 3,
                       pt.size = 0.2,
                       raster.dpi = c(2048,2038)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8, lineheight = 1),
        legend.position = "left",
        legend.key.spacing.y = unit(6, "pt")) &
  xlim(c(-10,25)) & ylim(c(-10, 20)) & 
  scale_color_manual(values = cols)

cd4_dimplot


signatures <- list(
  il_10_sig = il_10_ensembl)

CD4_cells_mochi <- AddModuleScore_UCell(CD4_cells_mochi, 
                                        features = signatures, name = NULL)

cd4_feature_plot <- FeaturePlot(CD4_cells_mochi, 
                                features = c("il_10_sig"),
                                reduction = "umap",
                                pt.size = 0.2,
                                alpha = 0.8,
                                cols = c("azure2", "red"),
                                split.by = "age_group",
                                combine = T,
                                repel = T,
                                keep.scale = "feature",
                                order = T,
                                raster = F,
                                raster.dpi = c(512, 512)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") & 
  scale_color_viridis_b(option = "B",
                        limits = c(0,0.25),
                        n.breaks = 10) & xlim(c(-10,25)) & 
  ylim(c(-10, 20))

cd4_feature_plot




cd4_dimplot <- ggrastr::rasterise(cd4_dimplot[[1]], 
                                  layers = "Point", dpi = 600)
cd4_dimplot

cd4_feature_plot_1 <- ggrastr::rasterise(cd4_feature_plot[[1]], 
                                         layers = "Point", dpi = 600)
cd4_feature_plot_1


cd4_feature_plot_2 <- ggrastr::rasterise(cd4_feature_plot[[2]], 
                                         layers = "Point", dpi = 600)
cd4_feature_plot_2

wrap_plots(cd4_dimplot,
           cd4_feature_plot_1,
           cd4_feature_plot_2,
           widths = c(1, 1, 1))



ggsave(
  "CD4 T-Cells Blood IL-10 Signature MMoCHi Annot.pdf",
  path = "output_data\\umap_human_aging/",
  width = 15,
  height = 4,
  units = "in"
)


all_dimplot <- DimPlot(blood, 
                       reduction = "umap", 
                       group.by = "MMoCHi_annotation",
                       raster = F,
                       label = F,
                       label.size = 3,
                       pt.size = 0.02,
                       raster.dpi = c(2048,2038)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8, lineheight = 1),
        legend.position = "left",
        legend.key.spacing.y = unit(6, "pt")) &
  xlim(c(-10,25)) & ylim(c(-10, 20))

all_dimplot


signatures <- list(
  il_10_sig = il_10_ensembl)

blood <- AddModuleScore_UCell(blood, features = signatures, name=NULL)

blood_feature_plot <- FeaturePlot(blood, 
                                  features = c("il_10_sig"),
                                  reduction = "umap",
                                  pt.size = 0.02,
                                  alpha = 0.8,
                                  cols = c("azure2", "red"),
                                  split.by = "age_group",
                                  combine = T,
                                  repel = T,
                                  keep.scale = "feature",
                                  order = T,
                                  raster = F,
                                  raster.dpi = c(512, 512)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") & 
  scale_color_viridis_b(option = "B",
                        limits = c(0,0.25),
                        n.breaks = 10) & xlim(c(-10,25)) & 
  ylim(c(-10, 20))

blood_feature_plot


wrap_plots(all_dimplot, 
           blood_feature_plot, 
           widths = c(0.8, 2))  ##Save as 11x4 Inches in Landscape Mode


ggsave(
  "All T-Cells Blood IL-10 Signature MMoCHi Annot.pdf",
  path = "output_data\\umap_human_aging/",
  width = 11,
  height = 3.5,
  units = "in"
)

## Load and Inspect BM T-Cell Seurat Object ------------------------------------

bm <- readRDS("input_data\\bone_marrow_human.rds") ## load T-cells in the BM

Idents(bm) <- "cell_type"

DimPlot(bm, 
        reduction = "umap", 
        group.by = "cell_type",
        raster = T,
        pt.size = 2,
        raster.dpi = c(1024,1024)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 

DimPlot(bm, 
        reduction = "umap", 
        group.by = "MMoCHi_annotation",
        raster = T,
        pt.size = 2,
        raster.dpi = c(1024,1024)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 


mochi_all_dimplot <- DimPlot(bm, 
                             reduction = "umap", 
                             group.by = "MMoCHi_annotation",
                             raster = F,
                             alpha = 0.5,
                             raster.dpi = c(2048,2048),
                             pt.size = 0.01) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 

mochi_all_dimplot


mochi_all_dimplot <- ggrastr::rasterise(mochi_all_dimplot,
                                        layers='Point',
                                        dpi = 600)
mochi_all_dimplot

CD8_tcells <- c("cd8_temra", "cd8_naive", "cd8_tem", 
                "cd8_trm", "cd8_cm", "cd8_mait")

CD4_tcells <- c("cd4_cm", "cd4_naive", "cd4_tem", 
                "cd4_treg", "cd4_temra", "cd4_trm")

Gamma_Delta <- c("cd3_gd")

bm@meta.data$cell_type_broad[bm@meta.data$MMoCHi_annotation %in% CD4_tcells] <- "CD4+ T Cell"
bm@meta.data$cell_type_broad[bm@meta.data$MMoCHi_annotation %in% CD8_tcells] <- "CD8+ T Cell"
bm@meta.data$cell_type_broad[bm@meta.data$MMoCHi_annotation %in% Gamma_Delta] <- "Gamma-Delta T Cell"

unique(bm@meta.data$cell_type_broad)

mochi_cd4_8_gd_dimplot <- DimPlot(bm, 
                                  reduction = "umap", 
                                  group.by = "cell_type_broad",
                                  raster = F,
                                  alpha = 0.5,
                                  pt.size = 0.01,
                                  raster.dpi = c(2048,2048)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 

mochi_cd4_8_gd_dimplot

mochi_cd4_8_gd_dimplot <- ggrastr::rasterise(mochi_cd4_8_gd_dimplot,
                                             layers='Point',
                                             dpi = 600)
mochi_cd4_8_gd_dimplot


genes <- c("CD8A", "CD4")

annots <- AnnotationDbi::select(org.Hs.eg.db, 
                                keys = genes, 
                                columns = "ENSEMBL", 
                                keytype = "SYMBOL")

annots$ENSEMBL
intersect(row.names(bm), annots$ENSEMBL)


cd8_featureplot <- FeaturePlot(bm, ## Confirm CD4 and CD8 identity
                               features = c("ENSG00000153563"),
                               reduction = "umap",
                               pt.size = 0.01,
                               combine = T,
                               alpha = 0.5,
                               repel = T,
                               order = T,
                               raster = F,
                               raster.dpi = c(2048, 2048)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") &
  scale_color_viridis_b(option = "C", limits = c(0, 1), 
                        n.breaks = 8, direction = 1)


cd8_featureplot <- ggrastr::rasterise(cd8_featureplot,
                                      layers='Point',
                                      dpi = 600)
cd8_featureplot

cd4_featureplot <- FeaturePlot(bm, ## Confirm CD4 and CD8 identity
                               features = c("ENSG00000010610"),
                               reduction = "umap",
                               pt.size = 0.01,
                               combine = T,
                               alpha = 0.6,
                               raster = F,
                               raster.dpi = c(2048, 2048),
                               order = T) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") &
  scale_color_viridis_b(option = "C", limits = c(0, 0.3), 
                        n.breaks = 8, direction = 1)

cd4_featureplot <- ggrastr::rasterise(cd4_featureplot,layers='Point',dpi = 600)
cd4_featureplot


wrap_plots(mochi_cd4_8_gd_dimplot, 
           mochi_all_dimplot, 
           cd8_featureplot, 
           cd4_featureplot,
           widths = c(1,1,1,1))


ggsave(
  "T-Cells in Bone Marrow.pdf",
  path = "output_data\\umap_human_aging",
  width = 20,
  height = 6,
  units = "in"
)

## BM T-Cell Aging Analysis ----------------------------------------------------

unique(bm@meta.data$donor_age)

old <- c("56 years",
         "67 years",
         "70-74 years",
         "68 years",
         "64 years",
         "50-54 years",
         "40-44 years",
         "48 years",
         "55-60 years",
         "75 years",
         "73 years")

young <- c("20-24 years",
           "23 years",
           "25 years",
           "35-39 years",
           "33 years",
           "27 years",
           "25-30 years",
           "20 years")


bm@meta.data$age_group[bm@meta.data$donor_age %in% old] <- "older"
bm@meta.data$age_group[bm@meta.data$donor_age %in% young] <- "younger"
bm@meta.data$age_group <- factor(bm@meta.data$age_group, levels = c("younger", "older"))

## Now this will allow to work with <40y vs >40y cohort


unique(bm@meta.data$age_group)
unique(bm@meta.data$cell_type)

head(bm@meta.data[["age_group"]], 20)
head(bm@meta.data[["donor_age"]], 20)

DimPlot(bm, 
        reduction = "umap", 
        group.by = "age_group",
        raster = T,
        pt.size = 4,
        raster.dpi = c(1024,1024)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 

DimPlot(bm, 
        reduction = "umap", 
        group.by = "cell_type",
        raster = T,
        label = F,
        pt.size = 4,
        raster.dpi = c(1024,1024)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 


Idents(bm) <- "MMoCHi_annotation"
Idents(bm)

CD8_cells_mochi <- subset(bm, idents = CD8_tcells)
CD8_cells_mochi

CD8_cells_mochi$MMoCHi_annotation <- stringr::str_wrap(
  CD8_cells_mochi$MMoCHi_annotation, width = 15)


cols <- c(
  "cd3_gd" = "gray",
  "cd8_trm" = "orange3",
  "cd4_cm" = "deepskyblue2",
  "cd8_naive" = "deeppink",
  "cd4_tem" = "aquamarine",
  "cd4_trm" = "skyblue4",
  "cd8_tem" = "salmon",
  "cd4_treg" = "cornflowerblue",
  "cd8_temra" = "red",
  "cd4_naive" = "cyan2",
  "cd8_mait" = "orange",
  "cd4_temra" = "slateblue2",
  "cd8_cm" = "yellow3"
)

cd8_dimplot <- DimPlot(CD8_cells_mochi, 
                       reduction = "umap", 
                       group.by = "MMoCHi_annotation",
                       raster = F,
                       alpha = 0.5,
                       label = T,
                       raster.dpi = c(2048,2048),
                       pt.size = 0.01) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "left") & xlim(-10, 20) & ylim(-10, 20) & 
  scale_color_manual(values = cols)

cd8_dimplot <- ggrastr::rasterise(cd8_dimplot, layers='Point',dpi = 600)
cd8_dimplot

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

annots <- AnnotationDbi::select(org.Hs.eg.db, 
                                keys = il_10_human, 
                                columns = "ENSEMBL", 
                                keytype = "SYMBOL")

il_10_ensembl <- annots$ENSEMBL
il_10_ensembl
il_10_ensembl <- intersect(il_10_ensembl, row.names(bm))

signatures <- list(
  il_10_sig = il_10_ensembl)

CD8_cells_mochi <- AddModuleScore_UCell(CD8_cells_mochi, 
                                        features = signatures, 
                                        name = NULL)

cd8_feature_plot <- FeaturePlot(CD8_cells_mochi, 
                                features = c("il_10_sig"),
                                reduction = "umap",
                                pt.size = 0.01,
                                alpha = 0.5,
                                cols = c("azure2", "red"),
                                split.by = "age_group",
                                combine = T,
                                repel = T,
                                keep.scale = "feature",
                                order = T,
                                raster = F,
                                raster.dpi = c(512, 512)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") & 
  scale_color_viridis_b(option = "B",
                        limits = c(0,0.4),
                        n.breaks = 10) & 
  xlim(-10, 20) & ylim(-10, 20)

cd8_feature_plot_1 <- ggrastr::rasterise(cd8_feature_plot[[1]], 
                                         layers = "Point", 
                                         dpi = 600)
cd8_feature_plot_1


cd8_feature_plot_2 <- ggrastr::rasterise(cd8_feature_plot[[2]], 
                                         layers = "Point", 
                                         dpi = 600)
cd8_feature_plot_2


wrap_plots(cd8_dimplot, 
           cd8_feature_plot_1,
           cd8_feature_plot_2,
           widths = c(1, 1, 1))

ggsave(
  "CD8 T-Cells Bone Marrow IL-10 Signature MMoCHi Annot.pdf",
  path = "output_data\\umap_human_aging",
  width = 15,
  height = 4,
  units = "in"
)


CD4_cells_mochi <- subset(bm, idents = CD4_tcells)

CD4_cells_mochi$MMoCHi_annotation <- stringr::str_wrap(
  CD4_cells_mochi$MMoCHi_annotation, width = 15)


cd4_dimplot <- DimPlot(CD4_cells_mochi, 
                       reduction = "umap", 
                       group.by = "MMoCHi_annotation",
                       raster = F,
                       alpha = 0.5,
                       label = T,
                       raster.dpi = c(2048,2048),
                       pt.size = 0.01) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "left") & xlim(-10, 20) & ylim(-10, 20) & 
  scale_color_manual(values = cols)

cd4_dimplot <- ggrastr::rasterise(cd4_dimplot, layers = "Point", dpi = 600)
cd4_dimplot

signatures <- list(
  il_10_sig = il_10_ensembl)

CD4_cells_mochi <- AddModuleScore_UCell(CD4_cells_mochi, 
                                        features=signatures, name=NULL)

cd4_feature_plot <- FeaturePlot(CD4_cells_mochi, 
                                features = c("il_10_sig"),
                                reduction = "umap",
                                pt.size = 0.01,
                                alpha = 0.5,
                                cols = c("azure2", "red"),
                                split.by = "age_group",
                                combine = T,
                                repel = T,
                                keep.scale = "feature",
                                order = T,
                                raster = F,
                                raster.dpi = c(512, 512)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") & 
  scale_color_viridis_b(option = "B",
                        limits = c(0,0.3),
                        n.breaks = 10) & xlim(-10, 20) & ylim(-10, 20)

cd4_feature_plot_1 <- ggrastr::rasterise(cd4_feature_plot[[1]], 
                                         layers = "Point", dpi = 600)
cd4_feature_plot_1

cd4_feature_plot_2 <- ggrastr::rasterise(cd4_feature_plot[[2]], 
                                         layers = "Point", dpi = 600)
cd4_feature_plot_2


wrap_plots(cd4_dimplot, 
           cd4_feature_plot_1, 
           cd4_feature_plot_2,
           widths = c(1, 1, 1)) 


ggsave(
  "CD4 T-Cells Bone Marrow IL-10 Signature MMoCHi Annot.pdf",
  path = "immune_aging_16_genes\\Plots\\Bone Marrow",
  width = 15,
  height = 4,
  units = "in"
)


cols <- c(
  "cd3_gd" = "gray",
  "cd8_trm" = "orange3",
  "cd4_cm" = "deepskyblue2",
  "cd8_naive" = "deeppink",
  "cd4_tem" = "aquamarine",
  "cd4_trm" = "skyblue4",
  "cd8_tem" = "salmon",
  "cd4_treg" = "cornflowerblue",
  "cd8_temra" = "red",
  "cd4_naive" = "cyan2",
  "cd8_mait" = "orange",
  "cd4_temra" = "slateblue2",
  "cd8_cm" = "yellow3"
)


all_dimplot <- DimPlot(bm, 
                       reduction = "umap", 
                       group.by = "MMoCHi_annotation",
                       raster = F,
                       label = F,
                       label.size = 3,
                       pt.size = 0.01,
                       raster.dpi = c(2048,2038)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8, lineheight = 1),
        legend.position = "left",
        legend.key.spacing.y = unit(6, "pt")) &
  xlim(c(-10,25)) & ylim(c(-10, 20)) & scale_color_manual(values = cols)

all_dimplot <- ggrastr::rasterise(all_dimplot, layers = "Point", dpi = 600)
all_dimplot

signatures <- list(
  il_10_sig = il_10_ensembl)

bm <- AddModuleScore_UCell(bm, features=signatures, name=NULL)

bm_feature_plot <- FeaturePlot(bm, 
                               features = c("il_10_sig"),
                               reduction = "umap",
                               pt.size = 0.01,
                               alpha = 0.5,
                               split.by = "age_group",
                               combine = T,
                               repel = T,
                               keep.scale = "feature",
                               order = T,
                               raster = F,
                               raster.dpi = c(512, 512)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") & 
  scale_color_viridis_b(option = "B",
                        limits = c(0,0.4),
                        n.breaks = 10) & xlim(c(-10,25)) & ylim(c(-10, 20))

bm_feature_plot

bm_feature_plot_1 <- ggrastr::rasterise(bm_feature_plot[[1]], 
                                        layers = "Point", dpi = 600)
bm_feature_plot_1

bm_feature_plot_2 <- ggrastr::rasterise(bm_feature_plot[[2]], 
                                        layers = "Point", dpi = 600)
bm_feature_plot_2




bm_umap <- wrap_plots(all_dimplot, 
                      bm_feature_plot_1,
                      bm_feature_plot_2,
                      widths = c(1, 1, 1))

bm_umap

saveRDS(bm_umap, "bm_umap.rds") 


ggsave(
  "All T-Cells Bone Marrow IL-10 Gene.pdf",
  path = "output_data\\umap_human_aging",
  width = 20,
  height = 7,
  units = "in"
)

rm(bm)

## Load and Inspect Spleen T-Cell Seurat Object --------------------------------

spleen <- readRDS("input_data\\spleen_human.rds") ## Load T-cells in the spleen

Idents(spleen) <- "cell_type"

DimPlot(spleen, 
        reduction = "umap", 
        group.by = "cell_type",
        raster = T,
        pt.size = 2,
        raster.dpi = c(1024,1024)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 

DimPlot(spleen, 
        reduction = "umap", 
        group.by = "MMoCHi_annotation",
        raster = T,
        pt.size = 2,
        raster.dpi = c(1024,1024)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 

mochi_all_dimplot <- DimPlot(spleen, 
                             reduction = "umap", 
                             group.by = "MMoCHi_annotation",
                             raster = F,
                             alpha = 0.5,
                             raster.dpi = c(2048,2048),
                             pt.size = 0.01) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 

mochi_all_dimplot


mochi_all_dimplot <- ggrastr::rasterise(mochi_all_dimplot,
                                        layers='Point',
                                        dpi = 600)
mochi_all_dimplot

CD8_tcells <- c("cd8_temra", "cd8_naive", "cd8_tem", 
                "cd8_trm", "cd8_cm", "cd8_mait")
CD4_tcells <- c("cd4_cm", "cd4_naive", "cd4_tem", 
                "cd4_treg", "cd4_temra", "cd4_trm")
Gamma_Delta <- c("cd3_gd")

spleen@meta.data$cell_type_broad[spleen@meta.data$MMoCHi_annotation %in% CD4_tcells] <- "CD4+ T Cell"
spleen@meta.data$cell_type_broad[spleen@meta.data$MMoCHi_annotation %in% CD8_tcells] <- "CD8+ T Cell"
spleen@meta.data$cell_type_broad[spleen@meta.data$MMoCHi_annotation %in% Gamma_Delta] <- "Gamma-Delta T Cell"

unique(spleen@meta.data$cell_type_broad)

mochi_cd4_8_gd_dimplot <- DimPlot(spleen, 
                                  reduction = "umap", 
                                  group.by = "cell_type_broad",
                                  raster = F,
                                  alpha = 0.5,
                                  pt.size = 0.01,
                                  raster.dpi = c(2048,2048)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 

mochi_cd4_8_gd_dimplot

mochi_cd4_8_gd_dimplot <- ggrastr::rasterise(mochi_cd4_8_gd_dimplot,
                                             layers='Point',dpi = 600)
mochi_cd4_8_gd_dimplot


genes <- c("CD8A", "CD4")

annots <- AnnotationDbi::select(org.Hs.eg.db, keys=genes, 
                                columns="ENSEMBL", keytype="SYMBOL")

annots$ENSEMBL
intersect(row.names(spleen), annots$ENSEMBL)


cd8_featureplot <- FeaturePlot(spleen, ## Confirm CD4 and CD8 identity
                               features = c("ENSG00000153563"),
                               reduction = "umap",
                               pt.size = 0.01,
                               combine = T,
                               alpha = 0.5,
                               repel = T,
                               order = T,
                               raster = F,
                               raster.dpi = c(2048, 2048)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") &
  scale_color_viridis_b(option = "C", limits = c(0, 1), 
                        n.breaks = 8, direction = 1)


cd8_featureplot <- ggrastr::rasterise(cd8_featureplot,layers='Point',dpi = 600)
cd8_featureplot

cd4_featureplot <- FeaturePlot(spleen, ## Confirm CD4 and CD8 identity
                               features = c("ENSG00000010610"),
                               reduction = "umap",
                               pt.size = 0.01,
                               combine = T,
                               alpha = 0.6,
                               raster = F,
                               raster.dpi = c(2048, 2048),
                               order = T) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") &
  scale_color_viridis_b(option = "C", limits = c(0, 0.3), n.breaks = 8, direction = 1)

cd4_featureplot <- ggrastr::rasterise(cd4_featureplot,layers='Point',dpi = 600)
cd4_featureplot


wrap_plots(mochi_cd4_8_gd_dimplot, mochi_all_dimplot, 
           cd8_featureplot, cd4_featureplot,
           widths = c(1,1,1,1))


ggsave(
  "T-Cells in Spleen.pdf",
  path = "output_data\\umap_human_aging",
  width = 20,
  height = 6,
  units = "in"
)

## Spleen T-Cell Aging Analysis ------------------------------------------------

unique(spleen@meta.data$donor_age)

old <- c("56 years",
         "67 years",
         "70-74 years",
         "68 years",
         "64 years",
         "50-54 years",
         "40-44 years",
         "48 years",
         "55-60 years",
         "75 years",
         "73 years")

young <- c("20-24 years",
           "23 years",
           "25 years",
           "35-39 years",
           "33 years",
           "27 years",
           "25-30 years",
           "20 years")


spleen@meta.data$age_group[spleen@meta.data$donor_age %in% old] <- "older"
spleen@meta.data$age_group[spleen@meta.data$donor_age %in% young] <- "younger"
spleen@meta.data$age_group <- factor(spleen@meta.data$age_group, levels = c("younger", "older"))

## Now this will allow to work with <40y vs >40y cohort

unique(spleen@meta.data$age_group)
unique(spleen@meta.data$cell_type)

head(spleen@meta.data[["age_group"]], 20)
head(spleen@meta.data[["donor_age"]], 20)

DimPlot(spleen, 
        reduction = "umap", 
        group.by = "age_group",
        raster = T,
        pt.size = 4,
        raster.dpi = c(1024,1024)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 

DimPlot(spleen, 
        reduction = "umap", 
        group.by = "cell_type",
        raster = T,
        label = F,
        pt.size = 4,
        raster.dpi = c(1024,1024)) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8)) 


Idents(spleen) <- "MMoCHi_annotation"
Idents(spleen)

CD8_cells_mochi <- subset(spleen, idents = CD8_tcells)
CD8_cells_mochi

cols <- c(
  "cd3_gd" = "gray",
  "cd8_trm" = "orange3",
  "cd4_cm" = "deepskyblue2",
  "cd8_naive" = "deeppink",
  "cd4_tem" = "aquamarine",
  "cd4_trm" = "skyblue4",
  "cd8_tem" = "salmon",
  "cd4_treg" = "cornflowerblue",
  "cd8_temra" = "red",
  "cd4_naive" = "cyan2",
  "cd8_mait" = "orange",
  "cd4_temra" = "slateblue2",
  "cd8_cm" = "yellow3"
)

CD8_cells_mochi$MMoCHi_annotation <- stringr::str_wrap(
  CD8_cells_mochi$MMoCHi_annotation, width = 15)


cd8_dimplot <- DimPlot(CD8_cells_mochi, 
                       reduction = "umap", 
                       group.by = "MMoCHi_annotation",
                       raster = F,
                       alpha = 0.5,
                       label = T,
                       raster.dpi = c(2048,2048),
                       pt.size = 0.01) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "left") & xlim(-10, 20) & 
  ylim(-10, 20) & scale_color_manual(values = cols)

cd8_dimplot <- ggrastr::rasterise(cd8_dimplot, 
                                  layers='Point',dpi = 600)
cd8_dimplot

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
                 "CCNA2") ##get the Human IL-10 signature

annots <- AnnotationDbi::select(org.Hs.eg.db, keys=il_10_human, 
                                columns="ENSEMBL", keytype="SYMBOL")

il_10_ensembl <- annots$ENSEMBL
il_10_ensembl
il_10_ensembl <- intersect(il_10_ensembl, row.names(spleen))

signatures <- list(
  il_10_sig = il_10_ensembl)

CD8_cells_mochi <- AddModuleScore_UCell(CD8_cells_mochi, 
                                        features=signatures, name=NULL)

cd8_feature_plot <- FeaturePlot(CD8_cells_mochi, 
                                features = c("il_10_sig"),
                                reduction = "umap",
                                pt.size = 0.01,
                                alpha = 0.5,
                                cols = c("azure2", "red"),
                                split.by = "age_group",
                                combine = T,
                                repel = T,
                                keep.scale = "feature",
                                order = T,
                                raster = F,
                                raster.dpi = c(512, 512)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") & 
  scale_color_viridis_b(option = "B",
                        limits = c(0,0.4),
                        n.breaks = 10) & xlim(-10, 20) & ylim(-10, 20)

cd8_feature_plot_1 <- ggrastr::rasterise(cd8_feature_plot[[1]], 
                                         layers = "Point", 
                                         dpi = 600)
cd8_feature_plot_1


cd8_feature_plot_2 <- ggrastr::rasterise(cd8_feature_plot[[2]], 
                                         layers = "Point", 
                                         dpi = 600)
cd8_feature_plot_2

wrap_plots(cd8_dimplot, 
           cd8_feature_plot_1,
           cd8_feature_plot_2,
           widths = c(1, 1, 1))

ggsave(
  "CD8 T-Cells Spleen IL-10 Signature MMoCHi Annot.pdf",
  path = "output_data\\umap_human_aging",
  width = 15,
  height = 4,
  units = "in"
)

CD4_cells_mochi <- subset(spleen, idents = CD4_tcells)

CD4_cells_mochi$MMoCHi_annotation <- stringr::str_wrap(
  CD4_cells_mochi$MMoCHi_annotation, width = 15)


cd4_dimplot <- DimPlot(CD4_cells_mochi, 
                       reduction = "umap", 
                       group.by = "MMoCHi_annotation",
                       raster = F,
                       alpha = 0.5,
                       raster.dpi = c(2048,2048),
                       pt.size = 0.01) & 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "left") & xlim(-10, 20) & ylim(-10, 20) & 
  scale_color_manual(values = cols)

cd4_dimplot <- ggrastr::rasterise(cd4_dimplot, layers = "Point", dpi = 600)
cd4_dimplot

signatures <- list(
  il_10_sig = il_10_ensembl)

CD4_cells_mochi <- AddModuleScore_UCell(CD4_cells_mochi, features=signatures, name=NULL)

cd4_feature_plot <- FeaturePlot(CD4_cells_mochi, 
                                features = c("il_10_sig"),
                                reduction = "umap",
                                pt.size = 0.01,
                                alpha = 0.5,
                                cols = c("azure2", "red"),
                                split.by = "age_group",
                                combine = T,
                                repel = T,
                                keep.scale = "feature",
                                order = T,
                                raster = F,
                                raster.dpi = c(512, 512)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") & 
  scale_color_viridis_b(option = "B",
                        limits = c(0,0.3),
                        n.breaks = 10) & xlim(-10, 20) & ylim(-10, 20)

cd4_feature_plot_1 <- ggrastr::rasterise(cd4_feature_plot[[1]], 
                                         layers = "Point", dpi = 600)
cd4_feature_plot_1

cd4_feature_plot_2 <- ggrastr::rasterise(cd4_feature_plot[[2]], 
                                         layers = "Point", dpi = 600)
cd4_feature_plot_2


wrap_plots(cd4_dimplot, 
           cd4_feature_plot_1, 
           cd4_feature_plot_2,
           widths = c(1, 1, 1))

ggsave(
  "CD4 T-Cells Spleen IL-10 Signature MMoCHi Annot.pdf",
  path = "output_data\\umap_human_aging/",
  width = 15,
  height = 4,
  units = "in"
)

unique(spleen@meta.data$MMoCHi_annotation)


cols <- c(
  "cd3_gd" = "gray",
  "cd8_trm" = "orange3",
  "cd4_cm" = "deepskyblue2",
  "cd8_naive" = "deeppink",
  "cd4_tem" = "aquamarine",
  "cd4_trm" = "skyblue4",
  "cd8_tem" = "salmon",
  "cd4_treg" = "cornflowerblue",
  "cd8_temra" = "red",
  "cd4_naive" = "cyan2",
  "cd8_mait" = "orange",
  "cd4_temra" = "slateblue2",
  "cd8_cm" = "yellow3"
)

all_dimplot <- DimPlot(spleen, 
                       reduction = "umap", 
                       group.by = "MMoCHi_annotation",
                       raster = F,
                       alpha = 0.4,
                       label = T,
                       label.size = 3,
                       pt.size = 0.01,
                       raster.dpi = c(2048,2038)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8, lineheight = 1),
        legend.position = "left",
        legend.key.spacing.y = unit(6, "pt")) &
  xlim(c(-10,25)) & ylim(c(-10, 20)) & scale_color_manual(values = cols)

all_dimplot <- ggrastr::rasterise(all_dimplot, layers = "Point", dpi = 600)
all_dimplot

signatures <- list(
  il_10_sig = il_10_ensembl)

spleen <- AddModuleScore_UCell(spleen, features=signatures, name=NULL)

spleen_feature_plot <- FeaturePlot(spleen, 
                                   features = c("il_10_sig"),
                                   reduction = "umap",
                                   pt.size = 0.01,
                                   alpha = 0.4,
                                   split.by = "age_group",
                                   combine = T,
                                   repel = T,
                                   keep.scale = "feature",
                                   order = T,
                                   raster = F,
                                   raster.dpi = c(512, 512)) + 
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=8),
        legend.position = "right") & 
  scale_color_viridis_b(option = "B",
                        limits = c(0,0.35),
                        n.breaks = 10) & xlim(c(-10,25)) & ylim(c(-10, 20))

spleen_feature_plot

spleen_feature_plot_1 <- ggrastr::rasterise(spleen_feature_plot[[1]], 
                                            layers = "Point", dpi = 600)
spleen_feature_plot_1

spleen_feature_plot_2 <- ggrastr::rasterise(spleen_feature_plot[[2]], 
                                            layers = "Point", dpi = 600)
spleen_feature_plot_2


spleen_umap <- wrap_plots(all_dimplot, 
                          spleen_feature_plot_1,
                          spleen_feature_plot_2,
                          widths = c(1, 1, 1))

spleen_umap

saveRDS(spleen_umap, "spleen_umap.rds")

ggsave(
  "All T-Cells Spleen IL-10 Signature MMoCHi Annot.pdf",
  path = "output_data\\umap_human_aging/",
  width = 20,
  height = 7,
  units = "in"
)


## Plot UMAPs Together ---------------------------------------------------------


spleen_umap <- readRDS("spleen_umap.rds")

blood_umap <- readRDS("blood_umap.rds")

bm_umap <- readRDS("bm_umap.rds")


blood_umap / spleen_umap / bm_umap


ggsave(
  "umap_all_tissues1.pdf",
  path = "output_data\\umap_human_aging/",
  width = 10,
  height = 10,
  units = "in"
)


writeLines(capture.output(sessionInfo()), "sessionInfo_human_umaps.txt")

