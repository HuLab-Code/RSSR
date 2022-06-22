
library(dplyr)
library(Seurat)
library(patchwork)

# load the dataset
RGC_r1 <- CreateSeuratObject(counts = read1, project = "RGC_r1") # RGC_repeat1
RGC_r2 <- CreateSeuratObject(counts = read2, project = "RGC_r2") # RGC_repeat2
RGC <- merge(x = RGC_r1, y = RGC_r2)

# add the metadata
RGC <- AddMetaData(object = RGC, metadata = stage, col.name = 'stage')
RGC <- AddMetaData(object = RGC, metadata = round, col.name = 'round')

# round <- factor(c(rep("round1",ncol(read1)),
#                   rep("round2",ncol(read2))))
# 
# stage <- factor(c(rep("Reg",table(stage1)[[1]]), 
#                   rep("Sur",table(stage1)[[2]]),
#                   rep("Reg",table(stage2)[[1]]),
#                   rep("Sur",table(stage2)[[2]])))

RGC[["percent.mt"]] <- PercentageFeatureSet(RGC, pattern = "^mt-")

VlnPlot(RGC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "round", ncol = 3)
VlnPlot(RGC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "stage", ncol = 3)

plot1 <- FeatureScatter(RGC, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt",group.by = "stage")
plot2 <- FeatureScatter(RGC, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA",group.by = "stage")
plot1 + plot2
# ------------------------ 
meta <- data.frame(RGC@meta.data)
counts <- data.frame(RGC@assays$RNA@counts)

Pten <- colSums(counts[rownames(counts)=="Pten",])
Nefl <- colSums(counts[rownames(counts)=="Nefl",])
Nefm <- colSums(counts[rownames(counts)=="Nefm",])
Thy1 <- colSums(counts[rownames(counts)=="Thy1",])
Rbpms <- colSums(counts[rownames(counts)=="Rbpms",])
Slc17a6 <- colSums(counts[rownames(counts)=="Slc17a6",])
Sncg <- colSums(counts[rownames(counts)=="Sncg",])
Pou4f1 <- colSums(counts[rownames(counts)=="Pou4f1",])
Pou4f2 <- colSums(counts[rownames(counts)=="Pou4f2",])
Pou4f3 <- colSums(counts[rownames(counts)=="Pou4f3",])

RGC[["Pten_TPM"]] <- Pten
RGC[["Nefl_TPM"]] <- Nefl
RGC[["Nefm_TPM"]] <- Nefm
RGC[["Thy1_TPM"]] <- Thy1
RGC[["Rbpms_TPM"]] <- Rbpms
RGC[["Slc17a6_TPM"]] <- Slc17a6
RGC[["Sncg_TPM"]] <- Sncg
RGC[["Pou4f1_TPM"]] <- Pou4f1
RGC[["Pou4f2_TPM"]] <- Pou4f2
RGC[["Pou4f3_TPM"]] <- Pou4f3

marker <- counts[rownames(counts)=="Rbpms"|
                     rownames(counts)=="Slc17a6"|
                     rownames(counts)=="Sncg"|
                     rownames(counts)=="Thy1"|
                     rownames(counts)=="Nefl"|
                     rownames(counts)=="Nefm"|
                     rownames(counts)=="Pou4f1"|
                     rownames(counts)=="Pou4f2"|
                     rownames(counts)=="Pou4f3",]
RGC[["marker_TPM"]] <- colSums(marker)

# ------------------------ 
Pten_m <- meta[meta$Pten_TPM < 30,]
Rbpms_m <- meta[meta$Rbpms_TPM > 0,]
Slc17a6_m <- meta[meta$Slc17a6_TPM > 0,]
Sncg_m <- meta[meta$Sncg_TPM > 0,]
Thy1_m <- meta[meta$Thy1_TPM > 0,]
Nefl_m <- meta[meta$Nefl_TPM > 0,]
Nefm_m <- meta[meta$Nefm_TPM > 0,]
Pou4f1_m <- meta[meta$Pou4f1_TPM > 0,]
Pou4f2_m <- meta[meta$Pou4f2_TPM > 0,]
Pou4f3_m <- meta[meta$Pou4f3_TPM > 0,]

# ------------------------ Selecting the cells for analysis---------------------
RGC <- subset(RGC, subset = Pten_TPM < 30)
RGC <- subset(RGC, subset = marker_TPM > 0)
RGC <- subset(RGC, subset = nCount_RNA > 442813)
RGC <- subset(RGC, subset = percent.mt < 15)
RGC <- subset(RGC, subset = nFeature_RNA > 714)
meta <- data.frame(RGC@meta.data)

VlnPlot(RGC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "stage_r", ncol = 3)
plot1 <- FeatureScatter(RGC, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "stage_r")
plot2 <- FeatureScatter(RGC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "stage_r")
plot1 + plot2
# -------------------------------- Normalizing ---------------------------------
RGC <- NormalizeData(RGC)
RGC <- FindVariableFeatures(RGC, selection.method = "vst", nfeatures = 2000)
RGC <- ScaleData(RGC, features = rownames(RGC))
RGC <- RunPCA(RGC, features = VariableFeatures(object = RGC))
VizDimLoadings(RGC, dims = 1:2, reduction = "pca")
DimPlot(RGC, reduction = "pca")
ElbowPlot(RGC)
RGC <- FindNeighbors(RGC,features = VariableFeatures(object = RGC), dims = 1:20)
# ------------------ choose 5 clusters for further analysis  -------------------
RGC <- FindClusters(RGC, resolution = .2) 
head(RGC[[]])
RGC <- RunUMAP(object = RGC, dims = 1:15) # 15
RGC <- CellCycleScoring(RGC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(RGC, reduction = "umap", group.by = "seurat_clusters",pt.size = 1.8)
DimPlot(RGC, reduction = "umap", group.by = "stage", cols = c("#FB6A4A","#78C679"),pt.size = 1.8)
FeaturePlot(RGC,reduction = "umap", features = rev("nFeature_RNA"), pt.size = 1.8)
FeaturePlot(RGC,reduction = "umap", features = rev("percent.mt"), pt.size = 1.8)
# -------------------- The DEGs between different clusters --------------------- 
Idents(RGC) <- "seurat_clusters"
levels(RGC)
clusterRGC <- FindAllMarkers(object = RGC,logfc.threshold = 0.25)
clusterRGC_sig <- clusterRGC[clusterRGC$avg_log2FC > 0 &
                               clusterRGC$p_val_adj < 0.05,]

clusterRGC %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> top50
DoHeatmap(RGC, features = top10$gene) + NoLegend()
# -------------------- The DEGs between different RGCs --------------------- 
Idents(RGC) <- "stage"
levels(RGC)
RGC_diff <- FindAllMarkers(object = RGC, logfc.threshold = 0.)
RGC_diff_sig <- RGC_diff[RGC_diff$avg_log2FC > 0.25 & RGC_diff$p_val_adj < 0.05,]
