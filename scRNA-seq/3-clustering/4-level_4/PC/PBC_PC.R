### Integration and lineage marker exploration, as well as singleR of PC+GC (18/06/21)

## Libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(Seurat)
library(harmony)
library(writexl)
library(ggpubr)

suppressMessages(require(SingleCellExperiment))
suppressMessages(require(CHETAH))
suppressMessages(require(SingleR))
suppressMessages(require(celldex))
suppressMessages(require(scRNAseq))
suppressMessages(require(cowplot))
suppressMessages(require(scater))
suppressMessages(require(scuttle))
suppressMessages(require(scran))
suppressMessages(require(bluster))


### loading files
  # >PC
setwd("/Users/spalominoe/Documents/Single-Cell/TONSIL/00.Plasma Cell/New Data 2/")
PC <- readRDS("PC_level4_90621.rds")

  # >BC LZ
setwd("/Users/spalominoe/Documents/Single-Cell/TONSIL/00.Plasma Cell/PC_from_GC/")
prePBC <- readRDS("prePBC_LZ-clust4.rds")
prePBC <- subset(prePBC, subset = assay == "5P",invert=TRUE)

DimPlot(prePBC)
FeaturePlot(prePBC, c("CD9", "BCL2A1"),cols = c("lightgray","blue"), reduction="umap", label = TRUE)

prePBC@meta.data$predicted.id <- as.character(prePBC@meta.data$seurat_clusters)
prePBC@meta.data$predicted.id[prePBC@meta.data$predicted.id=="9"] <- "CD9"
prePBC@meta.data$predicted.id[prePBC@meta.data$predicted.id=="4"] <- "BCL2A1"

####
####
#### In this script I'm going to join plasma cells + prePBC_LZ4
data <- merge(PC, y = prePBC, add.cell.ids = c("PC", "PBC"),  project = "PC_PBC")
head(colnames(data))

table(data@meta.data$predicted.id)

#### Normalize data, and adjust 3P and multiome
seurat_list <- SplitObject(data, split.by = "assay")
seurat_list <- seurat_list[c("3P", "multiome")]
seurat_list <- purrr::map(
  seurat_list,
  FindVariableFeatures,
  nfeatures = 5000
)

hvg <- purrr::map(seurat_list, VariableFeatures)
shared_hvg <- intersect(hvg$`3P`, hvg$multiome)

### variables, features, and reductions
data@meta.data$assay <- as.factor(data@meta.data$assay)

data <- data%>%
  ScaleData(features = shared_hvg) %>%
  RunPCA(features = shared_hvg) %>%
  RunHarmony(group.by.vars = "assay", reduction ="pca", dims = 1:15)

#
data <- RunUMAP(data, reduction = "harmony", dims = 1:15)
data <- FindNeighbors(data, reduction = "harmony", dims = 1:15)

#
DimPlot(data, group.by= "assay", reduction="umap")
DimPlot(data, group.by= "gem_id", reduction="umap")

#
data <- FindClusters(data, resolution = 0.6)
DimPlot(data, reduction = "umap")
data <- FindVariableFeatures(data)

DimPlot(data, group.by = "age_group")
DimPlot(data, reduction = "umap", group.by = "predicted.id", label = TRUE)
FeaturePlot(data,c("pct_mt","pct_ribosomal", "nCount_RNA","nFeature_RNA"),cols = c("lightgray","blue"), reduction="umap")

sort(table(data$predicted.id))

data@meta.data$PBC_PC <- as.character(data@meta.data$predicted.id)

Idents(data) <- "PBC_PC"
save.image(file = "PBC_PC2.RData")


DimPlot(data, cells.highlight = (WhichCells(data, idents = c("IgG plasmablast"))))


# Find markers
markers.PBC_PC <- FindAllMarkers(
  data, id = data@meta.data$predicted.id,
  logfc.threshold = 0.5,
  test.use = "wilcox",
  only.pos = TRUE,
  verbose = TRUE
)
markers.PBC_PC <- markers.PBC_PC %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), group_by = TRUE) %>%
  ungroup()



avg <- AverageExpression(data, features = NULL, add.ident = NULL, group.by = "PBC_PC" , return.seurat = TRUE, verbose = TRUE)

counts <- GetAssayData(avg, assay = "RNA")


## B-cell like markers
HLA <- markers.PBC_PC[grepl("^HLA", markers.PBC_PC$gene),]
HLA <- unique(HLA$gene)

Bcell <- c( "CD79B", "CD52", "CD19", "CD22", "CD83",  "CD72", "CD38", "REL", "PAX5", "BACH2", "CIITA", "FCMR", "AFF3", "SPIB", "MS4A1", "TRIM22", "LMO2", 
            "RGS13", "EML6", "CCDC88A", "SYNE2", 
            "MEF2B", "CD22", "EBF1", "CD22", "SERPINA9", "CD30")
Bcell <- c(Bcell, HLA)

Bcell_like <- markers.PBC_PC[markers.PBC_PC$gene %in% Bcell,]
Bcell_like <- Bcell_like %>% group_by(cluster) %>% top_n(50, avg_log2FC)

library(ComplexHeatmap)
Bcell_values <- counts[rownames(counts) %in% Bcell,]
Heatmap(t(Bcell_values))

### PC like markers
Pcell <- c("IRF4", "XBP1", "PRDM1", "JCHAIN", "MZB1", "NFKB")
Pcell_like <- markers.PBC_PC[markers.PBC_PC$gene %in% Pcell,]
Pcell_like <- Pcell_like %>% group_by(cluster) %>% top_n(100, avg_log2FC)

Pcell_values <- counts[rownames(counts) %in% Pcell,]
Heatmap(t(Pcell_values))

### Cycling markers
Cycling <- c("MKI67","NEIL1","FGD6","LMO2", "AICDA", "MME", "SUGCT", "PCNA")
Cycling_like <- markers.PBC_PC[markers.PBC_PC$gene %in% Cycling,]
Cycling_like <- Cycling_like %>% group_by(cluster) %>% top_n(100, avg_log2FC)

Cycling_values <- counts[rownames(counts) %in% Cycling,]
Heatmap(t(Cycling_values))

ALL <-  c(Bcell, Pcell)
ALL_values <- counts[rownames(counts) %in% ALL,]
Heatmap(t(ALL_values))


PC_shiny <- seurat2shiny(data, assay = data@active.assay, slot      = "data"             ,
                         reduction = "umap")


saveRDS(PC_shiny, file = "PC_shiny.rds")
saveRDS(data, file = "PC_GC.rds")

GC0_igG <- FindMarkers(data, ident.1 = "GC Derived precursor 0", ident.2 = "IgG plasmablast", min.pct = 0.5, logfc.threshold = 0.5)
GC1_igG <- FindMarkers(data, ident.1 = "GC Derived precursor 1", ident.2 = "IgG plasmablast", min.pct = 0.5, logfc.threshold = 0.5)
GC2_igG <- FindMarkers(data, ident.1 = "GC Derived precursor 2", ident.2 = "IgG plasmablast", min.pct = 0.5, logfc.threshold = 0.5)
GC3_igG <- FindMarkers(data, ident.1 = "GC Derived precursor 3", ident.2 = "IgG plasmablast", min.pct = 0.5, logfc.threshold = 0.5)

library(xlsx)
write.xlsx(GC2_igG, file="GC_precursors_vs_IgG_plasmablasts.xlsx", sheetName = "2_IgG", 
           col.names = TRUE, row.names = TRUE, append = TRUE)




####==========================================
# Subclustering cluster BCL2A1 (proliferative)
data_sub <- FindSubCluster(
  data,
  cluster = "BCL2A1",
  graph.name = "RNA_snn",
  resolution = 0.1,
  subcluster.name = "annotation_level_4"
)
Idents(data_sub) <- "annotation_level_4"

table(data_sub@meta.data$annotation_level_4)

DimPlot(data_sub, group.by="annotation_level_4", label = TRUE)


doubts <- subset(x = data, idents = c("CD9_0", "CD9_1", "Proliferating plasmablast"))
doubts2 <- subset(x = data, idents = c("CD9_1", "Proliferating plasmablast"))


# We get avarage expresion of each marker across clusters
data_sub@meta.data$merged <- as.character(data_sub@meta.data$annotation_level_4)

Idents(data_sub) < -"merged"

avg <- AverageExpression(data_sub, features = NULL, add.ident = NULL, group.by = "predicted.id" , return.seurat = TRUE, verbose = TRUE)

HLA <- markers.PBC_PC[grepl("^HLA", markers.PBC_PC$gene),]
HLA <- unique(HLA$gene)

Bcell <- c( "CD79B", "CD52", "CD19", "CD22", "CD83",  "CD72", "CD38", "REL", "PAX5", "BACH2", "CIITA")
Bcell <- c(Bcell, HLA)

Bcell_like <- markers.PBC_PC[markers.PBC_PC$gene %in% Bcell,]
Bcell_like <- Bcell_like %>% group_by(cluster) %>% top_n(50, avg_log2FC)


DoHeatmap(avg, features = Bcell_like$gene, group.by = "ident", size = 3, angle = 0, combine = TRUE, draw.lines = FALSE) + 
  theme(axis.text.x = element_text(size = 7), axis.line = element_line(colour = "#ffffff")) + 
  scale_fill_gradient2(low = '#1000ff', mid = "#ffffff", high = '#aa0101', space = "Lab", na.value = "#ffffff", midpoint = 0, guide = "colourbar", aesthetics = "fill") + 
  coord_flip(expand = TRUE, clip = "on")

### PC like markers
Pcell <- c("IRF4", "IRF8", "XBP1", "PRDM1")
Pcell_like <- markers.PBC_PC[markers.PBC_PC$gene %in% Pcell,]
Pcell_like <- Pcell_like %>% group_by(cluster) %>% top_n(100, avg_log2FC)

DoHeatmap(avg, features = Pcell_like$gene, group.by = "ident", size = 3, angle = 0, combine = TRUE, draw.lines = FALSE) + 
  theme(axis.text.x = element_text(size = 7), axis.line = element_line(colour = "#ffffff")) + 
  scale_fill_gradient2(low = '#1000ff', mid = "#ffffff", high = '#aa0101', space = "Lab", na.value = "#ffffff", midpoint = 0, guide = "colourbar", aesthetics = "fill") + 
  coord_flip(expand = TRUE, clip = "on")
#############
#############
#############

# single R . We are going to transfer the labels from the PC dataset to the PBC dataset
# > Prepare data
    # USED AS REFERENCE
PCsingleR <- as.SingleCellExperiment(PC)
table(PCsingleR$predicted.id)
sort(table(PCsingleR$predicted.id))

  # USED AS TEST
BCsingleR <- as.SingleCellExperiment(prePBC)
sort(table(BCsingleR$predicted.id))

#----
# > Perform predictions
predictions <- SingleR(test=BCsingleR, 
                       ref=PCsingleR, labels=PCsingleR$predicted.id)

plotScoreHeatmap(predictions)
table(predictions$labels)

plotDeltaDistribution(predictions)

#---
# > Cluster level annotation
dec <- modelGeneVarByPoisson(BCsingleR)
sce <- denoisePCA(BCsingleR, dec, subset.row=getTopHVGs(dec, n=5000))

colLabels(sce) <- clusterRows(reducedDim(sce), NNGraphParam())

set.seed(117)
sce <- runTSNE(sce, dimred="PCA")
plotTSNE(sce, colour_by="label")
plotTSNE(sce, colour_by=I(predictions$labels), text_colour="red")

SingleR(sce, PCsingleR, clusters=colLabels(sce), labels=PCsingleR$predicted.id)

save.image(file = "PBC_PC.RData")

##
library(scater)
collected <- list()
all.markers <- metadata(predictions)$de.genes

BCsingleR$labels <- predictions$labels

setwd("/Users/spalominoe/Documents/Single-Cell/TONSIL/00.Plasma Cell/PC_from_GC/")

pdf("Proliferative.pdf", width=20,height=60)
proliferating <- plotHeatmap(BCsingleR, order_columns_by="labels",
            features=(unique(unlist(all.markers$`Proliferating plasmablast`))), cex.lab=3, cex.axis=3, pointsize=3, main= "Proliferating plasmablasts")

dev.off()


pdf("GC Derived precursor 0.pdf", width=20,height=60)
GC_pre0 <- plotHeatmap(BCsingleR, order_columns_by="labels",
                             features=(unique(unlist(all.markers$`GC Derived precursor 0`))), cex.lab=2, cex.axis=2, pointsize=2, main= "GC Derived precursor 0")
dev.off()


pdf("GC Derived precursor 1.pdf", width=20,height=60)
GC_pre0 <- plotHeatmap(BCsingleR, order_columns_by="labels",
                       features=(unique(unlist(all.markers$`GC Derived precursor 1`))), cex.lab=2, cex.axis=2, pointsize=2, main= "GC Derived precursor 1")
dev.off()

pdf("IgG plasmablast.pdf", width=20,height=60)
IgG_plasmablast <- plotHeatmap(BCsingleR, order_columns_by="labels",
                       features=(unique(unlist(all.markers$`IgG plasmablast`))), cex.lab=2, cex.axis=2, pointsize=2, main= "IgG plasmablast")
dev.off()



pdf("IgD.pdf", width=20,height=60)
IgD <- plotHeatmap(BCsingleR, order_columns_by="labels",
                               features=(unique(unlist(all.markers$`IgD`))), cex.lab=2, cex.axis=2, pointsize=2, main= "IgD")
dev.off()
