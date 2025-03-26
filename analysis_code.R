#load in the necessary libraries for analysis 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#loading the dataset 
filename <- file.choose()
seurat_data1 <- readRDS(filename)
#look into the data
print(seurat_data1)
summary(seurat_data1)

#preform quality control 
#idenfity and filter low quality cells looking at number of genes detected, total counts per cell, proportion of mitochondrial genes
#[[ can add columns to object metadata 
#PercentageFeatureSet() calculates the mitochondrial QC metrics 
#pattern adjusted based on system of interest

seurat_data1[["percent.mt"]] <- PercentageFeatureSet(seurat_data1, pattern = "^mt-")
#show QC metrics for the first 5 cells 
head(seurat_data1@meta.data, 5)
head(rownames(seurat_data1))

#visualize QC metrics as a violin plot 
VlnPlot(seurat_data1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#filter based on QC metrics 

plot1 <- FeatureScatter(seurat_data1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_data1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#calculating thresholds based on IQR 
#example thresholds from satija lab is filter cells that have unique feature counts over 2500 and less than 200 
#filter cells that have >5% mitochondrial counts 

lower_bound_genes <- quantile(seurat_data1$nFeature_RNA, 0.25) - 1.5 * IQR(seurat_data1$nFeature_RNA)
lower_bound_genes
upper_bound_genes <- quantile(seurat_data1$nFeature_RNA,0.75) + 1.5 * IQR(seurat_data1$nFeature_RNA)
upper_bound_genes

lower_bound_counts <- quantile(seurat_data1$nCount_RNA, 0.25) - 1.5 * IQR(seurat_data1$nCount_RNA)
lower_bound_counts
upper_bound_counts <- quantile(seurat_data1$nCount_RNA, 0.75) + 1.5 * IQR(seurat_data1$nCount_RNA)
upper_bound_counts

#now applying the thresholds 
seurat_data1 <- subset(seurat_data1, subset = nFeature_RNA > lower_bound_genes & 
                         nFeature_RNA < upper_bound_genes &
                         nCount_RNA > lower_bound_counts &
                         percent.mt < 5)

#visualize QC metrics and scatterplots post filtering to ensure the low quality cells have been removed 

VlnPlot(seurat_data1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#filter based on QC metrics 

plot1 <- FeatureScatter(seurat_data1, feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("Counts vs. Mitochondrial Percentage (Filtered)")
plot2 <- FeatureScatter(seurat_data1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +  ggtitle("Counts vs. Features (Filtered)")
plot1 + plot2

#normalize the data
seurat_data1 <- SCTransform(seurat_data1, vars.to.regress = "percent.mt", verbose = FALSE)

#dimensionality reduction 
#pca 
seurat_data1 <- RunPCA(seurat_data1, verbose = FALSE)
#looking into the pca results 
#print(seurat_data1[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(seurat_data1, dims = 1:2, reduction = "pca")

#png("pca.png", width = 800, height = 600)
DimPlot(seurat_data1, reduction= "pca") + NoLegend()
ElbowPlot(seurat_data1)
#dev.off()

#umap plot
seurat_data1 <- RunUMAP(seurat_data1, dims = 1:30, verbose = FALSE)
seurat_data1 <- FindNeighbors(seurat_data1, dims = 1:30, verbose = FALSE)
seurat_data1 <- FindClusters(seurat_data1, verbose = FALSE)
DimPlot(seurat_data1, label = TRUE)

seurat_data1.markers <- FindAllMarkers(seurat_data1, only.pos = TRUE)
seurat_data1.markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(avg_log2FC > 1)

#visualize the canonical marker genes as violin plots 
VlnPlot(seurat_data1, features = c("Cd3e", "Cd8a", "Cd4", "Cd14", "Lyz2", "Nphs2", "Pecam1", 
                                   "Kdr", "Cd74", "Batf3", "Ccr2", "Cd81", "Itgam",
                                   "C1qc"), 
        pt.size = 0.2, ncol = 4)

FeaturePlot(seurat_data1, features = c("Cd3e", "Cd8a", "Cd4", "Cd14", "Lyz2", "Nphs2", "Pecam1", "Kdr",
                                       "Cd74", "Batf3", "Ccr2", "Cd81", "Itgam",
                                       "C1qc"))

#generating a expression heatmap for given cells and features, so here we are plotting the top 20 markers for each cluster 
seurat_data1.markers %>% 
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>% 
  slice_head(n=10) %>%
  ungroup() -> top10
DoHeatmap(seurat_data1, features = top10$gene) + NoLegend()





