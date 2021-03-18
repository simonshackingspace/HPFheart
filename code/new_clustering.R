library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)
library(ggplot2)
library(CountClust)

counts <- Read10X(data.dir = "/home/zz2565/data/thesis_data/sc_heart_seurat_matrix/",
                  gene.column=1)   # Seurat function to read in 10x count data

# add meta data
meta_data = read.csv("/home/zz2565/data/thesis_data/diff_seqs/meta_data_singlet.csv", 
                     header=TRUE,
                     stringsAsFactors=T,
                     sep=",",
                     row.names=1)

# create seurat object and run pca and tsne
pbmc <- CreateSeuratObject(counts = t(counts), project="Human Heart scRNA")

# inspect data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT.")
VlnPlot(pbmc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

#QC and selecting cells for further analysis
pbmc <- subset(pbmc,
               subset = nFeature_RNA > 1000 &
                 nFeature_RNA < 5000 &
                 percent.mt < 10)
pbmc_out <- pbmc

# normalize data
pbmc <- NormalizeData(pbmc)

# identification of highly variable features
pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", 
                             nfeatures = 2000)

# scale the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# UMAP
pbmc <- RunUMAP(pbmc, 
                features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "umap", label=T, label.size=3)

# pca
pbmc <- RunPCA(pbmc,
               features = VariableFeatures(object = pbmc))
# determine number of cluster
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc)
JackStrawPlot(pbmc)
ElbowPlot(pbmc)
# pca cluster
pbmc <- FindNeighbors(pbmc, dims = 1:5)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# UMAP with clusters
DimPlot(pbmc, reduction = "umap", label=T, label.size=3)

# compare cluster 2 and 4
cluster2.markers = FindMarkers(pbmc, ident.1=2, ident.2=4,
                               min.pct=0.25, 
                               logfc.threshold = 0.25)
head(cluster2.markers, 10)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)
    top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
markers <- c()
for (i in 1:length(levels(top10$cluster))){
  clusterX <- top10[top10$cluster==i-1,]
  markers <- cbind(markers, clusterX$gene)
}
colnames(markers) <- c(0:13)
write.csv(markers, "/home/zz2565/Documents/thesis_code/thesis_output/markers.csv", row.names=F)

FeaturePlot(pbmc, features = c("ACTA2","OGN","IRX4", "HES4", "NR2F2"), min.cutoff = "q9")
FeaturePlot(pbmc, features = c("IF", "ACTN2", "TNNI3", "CASP3", "ACTA2", "CKAP4","OGN", "ALDH1A2", "ITLN1"), min.cutoff = "q9")
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
FeaturePlot(pbmc, features = c("MYH6", "FLT4"), min.cutoff = "q9")

celltypes <- scan("/home/zz2565/Documents/thesis_code/thesis_output/celltypes.txt", 
                  character(),
                  sep="\t")
pbmc <- RenameIdents(object = pbmc,
             '0' = celltypes[1],
             '1' = celltypes[2],
             '2' = celltypes[3],
             '3' = celltypes[4],
             '4' = celltypes[5],
             '5' = celltypes[6],
             '6' = celltypes[7],
             '7' = celltypes[8],
             '8' = celltypes[9],
             '9' = celltypes[10],
             '10' = celltypes[11])
Idents(pbmc_out) <- Idents(pbmc)
pbmc <- subset(pbmc, idents=c("immune cells", "erythrocyte"), invert=T)
DimPlot(pbmc, reduction = "umap", label=T, label.size=3)

pbmc_out <- subset(pbmc_out, idents=c("immune cells", "erythrocyte"), invert=T)
data_to_write_out <- as.data.frame(t(as.matrix(GetAssayData(pbmc_out))))
write.csv(data_to_write_out, "/home/zz2565/data/thesis_data/diff_seqs/sc_data_filtered.csv")
ct_file <- as.data.frame(Idents(pbmc_out))
write.table(ct_file, "/home/zz2565/data/thesis_data/diff_seqs/celltypes.csv", sep=",", col.names=F)
