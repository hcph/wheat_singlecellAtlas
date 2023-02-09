###Step 1: Run cellranger to produce expression matrix
#built genome
cellranger mkref --memgb=40 --genome=ak58_genome --fasta=161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta --genes=161010_Chinese_Spring_v1.0_pseudomolecules_parts.filter.gtf
#count
cellranger count --id=ajbk_15000_cell_range_AK58R1 \
                   --transcriptome=/public/home/xinwang/snRNA/db/wheat_genomeR2 \
                   --fastqs=/public/home/xinwang/snRNA/fq \
                   --sample=huanongxiaomai1 \
                   --expect-cells=15000 \
                   --localcores=8 \
                   --localmem=64 \
                   --include-introns 
###Step 2: Cell clustering using Seurat in R
library(Seurat)
library(patchwork)
library(dplyr)
wheat.data<-Read10X_h5("filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
pbmc<- CreateSeuratObject(counts=wheat.data, project = "pig_single", min.cells=3, min.features=200)
pbmc
pdf("AK58_metrics.pdf")
#QC metrics
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#visualized QC merics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#visualize feature-feature relationship
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 100000 & percent.mt < 5)
#normalizing
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#identification of highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#scaling the data and remove unwanted sources of variation, a standard ore-processing step prior to dimensional reduction techniques like PCA
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt") #2000
#pbmc <- ScaleData(pbmc, features = all.genes) 
#perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))  
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
#Determine the 'dimensionality' of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc)
dev.off()
save(pbmc,file="pbmc.all.Rdata")
#cluster the cells
load(file = 'pbmc.all.Rdata')
#adjust the dim and resolution
rm(list = ls())
library(ggplot2)
load(file = 'pbmc.all.Rdata')
sce <- test
library(clustree)
sce <- FindClusters(
  object = sce,
  resolution = c(seq(.1,1.6,.2))
)
png("1dim49-2dim49.clustree.png", width=1000,height=700)
clustree(sce@meta.data, prefix = "RNA_snn_res.")
dev.off()
#Choose dim=49, resolution=0.8 for further analysis
test <- FindNeighbors(pbmc, dims = 1:49)  #adjust dim
test <- FindClusters(test, resolution = 0.8)
#adjust res
# Look at cluster IDs of the first 5 cells
head(Idents(test), 5)
#table(pbmc$seurat_clusters) 
#Run non-linear dimensional reduction (UMAP/tSNE)
test <- RunUMAP(test, dims = 1:49) 
png("all_1dim49-2dim49-0.9.png", width=800, height=800)
DimPlot(test, reduction = 'umap', label=TRUE, pt.size = 2.0)
dev.off()
save(test,file = 'basic.sce.pbmc.Rdata')
save(pbmc,file='pbmc.final.Rdata')
#draw tsne map
test <- RunTSNE(test,dims=1:49,check_duplicates = FALSE)
pdf("TSNE.pdf",width=1000, heigh=800)
TSNEPlot(test, label=TRUE, cols = allcolour, pt.size = 0.6, label.size=5, label.box=F)
dev.off()
#using ggplot to draw umap figure
sam<-AddMetaData(test,test@reductions$umap@cell.embeddings,col.name = colnames(test@reductions$umap@cell.embeddings))
head(sam@meta.data)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
              "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
              "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
              "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
class_avg <- sam@meta.data %>%
    group_by(RNA_snn_res.0.8) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
 umap <-  ggplot(sam@meta.data ,aes(x=UMAP_1,y=UMAP_2))+
    geom_point(aes(color=RNA_snn_res.0.8))+
    scale_color_manual(values = allcolour)+
    geom_text(aes(label = RNA_snn_res.0.8), data = class_avg)+
    theme(text=element_text(family="Arial",size=18)) +
    theme(panel.background = element_rect(fill='white', colour='black'), 
                     panel.grid=element_blank(), axis.title = element_text(color='black',
                                                      family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
                     axis.ticks = element_line(color='black'), 
                     axis.ticks.margin = unit(0.6,"lines"),
                     axis.line = element_line(colour = "black"), 
                     axis.title.x=element_text(colour='black', size=18),
                     axis.title.y=element_text(colour='black', size=18),
                     axis.text=element_text(colour='black',size=18),
                     legend.title=element_blank(),
                     legend.text=element_text(family="Arial", size=18),
                     legend.key=element_blank())+
    theme(plot.title = element_text(size=22,colour = "black",face = "bold"))  + 
    guides(colour = guide_legend(override.aes = list(size=5)))
##Draw 3D UMAP figure
#install.packages("plotly")
library(plotly)
library(Seurat)
# Re-run UMAPs that you have accurate calculations for all UMAP(s)
pbmc_small <- RunUMAP(pbmc, dims = 1:49,n.components = 3L)
head(pbmc_small[["umap"]]@cell.embeddings)
umap_1 <- pbmc_small[["umap"]]@cell.embeddings[,1]
umap_2 <- pbmc_small[["umap"]]@cell.embeddings[,2]
umap_3 <- pbmc_small[["umap"]]@cell.embeddings[,3]
plot.data <- FetchData(object = pbmc_small, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "RNA_snn_res.0.8"))
plot.data$label <- paste(plot.data$RNA_snn_res.0.8)
plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~RNA_snn_res.0.8, 
        colors = c("0"="#DC143C","1"="#0000FF","2"="#20B2AA","3"="#FFA500","4"="#9370DB","5"="#98FB98",
                     "6"="#F08080","7"="#1E90FF","8"="#7CFC00","9"="#FFFF00",
                   "10"="#808000","11"="#FF00FF","12"="#FA8072","13"="#7B68EE",
                   "14"="#9400D3","15"="#800080","16"="#A0522D","17"="#D2B48C",
                   "18"="#D2691E","19"="#87CEEB","20"="#40E0D0","21"="#5F9EA0"),
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 3, width=2), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names
#Find markers
pbmc <- test
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
top2 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
VlnPlot(pbmc, features =top2$gene, 
#save markers
write.csv(pbmc.markers,"1dim49-2dim49-0.7-all.marker.csv")
top<- pbmc.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top,"1dim49-2dim49-0.7-top50_marker.csv")
#top10 genes heatmap
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(pbmc))
plot1 = DoHeatmap(pbmc, features = top10, group.by = "seurat_clusters", group.bar = T, size = 4)
pdf("1dim49-2dim49-0.7_top10marker.pdf", width=900,height=500)
plot1
dev.off()
top100 <- pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 100, wt = avg_log2FC)
write.csv(top100,"1dim49-2dim49-0.7-top100_marker.csv")
top200 <- pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 200, wt = avg_log2FC)
write.csv(top200,"1dim49-2dim49-0.7-top200_marker.csv")
pdf("1dim49-2dim49-0.7_top2marker.pdf", width=1000,height=1000)
DotPlot(pbmc, features = unique(top2$gene)) + RotatedAxis()
dev.off()
