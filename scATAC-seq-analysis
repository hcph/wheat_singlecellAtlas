#Prepare genome feature file
library(GenomicFeatures)
library(BSgenomewheatpart)
spompe <- makeTxDbFromGFF('/public/workspace/xuyongyue/ZLH/scATAC-seq/161010_Chinese_Spring_v1.0_pseudomolecules_parts.filter.gtf')
#saveDb(spompe, file="/public/workspace/xuyongyue/ZLH/scATAC-seq/wheat_txdb.sqlite")
wheat_txdb <- loadDb("/public/workspace/xuyongyue/ZLH/scATAC-seq/wheat_txdb.sqlite")
seqlevels(wheat_txdb) <- c("chr1A_part1","chr1A_part2","chr1B_part1","chr1B_part2","chr1D_part1",
                           "chr1D_part2","chr2A_part1","chr2A_part2","chr2B_part1","chr2B_part2","chr2D_part1","chr2D_part2",
                           "chr3A_part1","chr3A_part2","chr3B_part1","chr3B_part2","chr3D_part1","chr3D_part2","chr4A_part1",
                           "chr4A_part2","chr4B_part1","chr4B_part2","chr4D_part1","chr4D_part2","chr5A_part1","chr5A_part2",
                           "chr5B_part1","chr5B_part2","chr5D_part1","chr5D_part2","chr6A_part1","chr6A_part2","chr6B_part1",
                           "chr6B_part2","chr6D_part1","chr6D_part2","chr7A_part1","chr7A_part2","chr7B_part1","chr7B_part2",
                           "chr7D_part1","chr7D_part2","chrUn")
gene.ranges <- genes(wheat_txdb); length(gene.ranges$gene_id)
gir = GRanges(symbol=gene.ranges$gene_id, gene.ranges)  #将gene_id名字改为symbol
length(gene.ranges$gene_id)
tss.ranges <- resize(gene.ranges, 1, "start")
exon.ranges <- exons(wheat_txdb)
geneAnnotation = createGeneAnnotation(
  TSS = tss.ranges,
  exons = exon.ranges,
  genes = gir
)
genomeAnnotation = createGenomeAnnotation(genome = wheat, filter=F)
#bulid BSgenome object
#download faToTwoBit toolkit
#wget  -c "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit"
#chmod 755 faToTwoBit
#faToTwoBit 161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta wheat.2bit
#seed file from http://www.bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf
#forgeBSgenomeDataPkg("/public/workspace/xuyongyue/ZLH/scATAC-seq/wheat_seed")
#build、check、install
#R CMD build /public/workspace/xuyongyue/ZLH/scATAC-seq/BSgenome.wheat_part
#R CMD check /public/workspace/xuyongyue/ZLH/scATAC-seq/BSgenomewheatpart_1.0.0.tar.gz
#R CMD INSTALL -l /public/workspace/xuyongyue/R/x86_64-pc-linux-gnu-library/3.6.0 /public/workspace/xuyongyue/ZLH/scATAC-seq/BSgenomewheatpart_1.0.0.tar.gz
library(BSgenomewheatpart)
genomeAnnotation = createGenomeAnnotation(genome = wheat, filter=F)

##scATAC-seq data analysis
bytlib load R-3.6.0
devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.0", repos = BiocManager::repositories())
library(ArchR)
library(Seurat)
addArchRChrPrefix(chrPrefix = FALSE)
inputFiles="/public/workspace/xuyongyue/ZLH/scATAC-seq/R2/ajbk_fragments.tsv.gz"
addArchRThreads(threads = 1)
subThreading = FALSE
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames ="ajbk",
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  minTSS = 4, 
  minFrags = 3000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  excludeChr = "chrUn"
)

#Inferring scATAC-seq Doublets with ArchR
#ArrorFiles = c("huada.arrow", "ajbk.arrow")
ArrorFiles = c("ajbk.arrow")
doubScores <- addDoubletScores(
  input = "ajbk.arrow",
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

#Creating an ArchRProject
projHeme1 <- ArchRProject(
  ArrowFiles = "ajbk.arrow",
  outputDirectory = "AK58",
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
projHeme1
#Manipulating An ArchRProject
#We can access the cell names associated with each cell:
head(projHeme1$cellNames)
#We can access the sample names associated with each cell:
head(projHeme1$Sample)
#We can access the TSS Enrichment Scores for each cell:
quantile(projHeme1$TSSEnrichment)
# Subsetting an ArchRProject by cells
#We can subset the project numerically, for example taking the first 100 cells in the project:
projHeme1[1:100, ]
#We can subset the project based on certain cell names:
projHeme1[projHeme1$cellNames[1:100], ]
#We can subset the project to keep all cells corresponding to a specific sample:
idxSample <- BiocGenerics::which(projHeme1$Sample %in% "ajbk")
cellsSample <- projHeme1$cellNames[idxSample]
projHeme1[cellsSample, ]
#We can subset the project to only keep cells that meet a specific cutoff for the TSS enrichment score
idxPass <- which(projHeme1$TSSEnrichment >= 8)
cellsPass <- projHeme1$cellNames[idxPass]
projHeme1[cellsPass, ]
#Plotting QC metrics - log10(Unique Fragments) vs TSS enrichment score
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)
#Make a ridge plot for each sample for the TSS enrichment scores.
p1 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
#Example 2. Make a violin plot for each sample for the TSS enrichment scores.
p2 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
#Make a ridge plot for each sample for the log10(unique nuclear fragments).
p3 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)
#Example 4. Make a violin plot for each sample for the log10(unique nuclear fragments).
p4 <- plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 4, height = 4)
# Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles
p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
#Saving and Loading an ArchRProject
saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-ProjHeme1", load = FALSE)
#projHeme1<-readRDS("Save-ProjHeme1/Save-ArchR-Project.rds")

#Filtering Doublets from an ArchRProject
#projHeme1<-readRDS(file="Save-ProjHeme1/Save-ArchR-Project.rds")
projHeme2 <- filterDoublets(projHeme1)
#Filtering 616 cells from ArchRProject!
#ajbk : 616 of 7851 (7.8%)
projHeme2
#numberOfCells(1): 7235
#medianTSS(1): 8.868
#medianFrags(1): 11060

#Dimensionality Reduction with ArchR
projHeme2 <- addIterativeLSI(
  ArchRProj = projHeme2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI4", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 15000, 
  dimsToUse = 1:30
)
#Clustering using Seurat’s FindClusters() function
projHeme2 <- addClusters(
  input = projHeme2,
  reducedDims = "IterativeLSI4",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.4
)
#To better understand which samples reside in which clusters, we can create a cluster confusion matrix across each sample using the confusionMatrix() function
cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
cM
#Single-cell Embeddings
#Uniform Manifold Approximation and Projection (UMAP)
projHeme2 <- addUMAP(
  ArchRProj = projHeme2, 
  reducedDims = "IterativeLSI4", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.7, 
  metric = "cosine"
)
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters-0.7.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "ScranClusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-ScranClustersR-0.7.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

#t-Stocastic Neighbor Embedding (t-SNE)
projHeme2 <- addTSNE(
  ArchRProj = projHeme2, 
  reducedDims = "IterativeLSI4", 
  name = "TSNE", 
  perplexity = 30
)
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-TSNE-Sample-Clusters-0.7.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "ScranClusters", embedding = "TSNE")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-tSNE-Sample-ScranClusters-0.7.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
#Dimensionality Reduction After Harmony
projHeme2<-projHeme2[-which(projHeme2@cellColData@listData$Clusters %in% c("C1","C2","C3","C4","C5")),]
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "R1Plot-UMAP-Sample-Clusters-0.7.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

#Calculating Gene Scores in ArchR
#Identifying Marker Genes
markersGS <- getMarkerFeatures(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6
write.csv(markerList,"marger0.7--fdr0.01-lof2fc1.25.csv")
markerGenes  <- c(
  "TraesCS1B02G400500",  "TraesCS2D02G349900",  "TraesCS4B02G079900",  "TraesCS2D02G349500",  "TraesCS6B02G012800",  "TraesCS4D02G078800",  "TraesCS5A02G385400",  "TraesCS2B02G369900",  "TraesCS2A02G351700",  "TraesCS4D02G120700",  "TraesCS5D02G395400",  "TraesCS2A02G351400",  "TraesCS3B02G108200",  "TraesCS7A02G453100",
  "TraesCS7B02G468200",  "TraesCS2B02G236500",  "TraesCS2D02G217100",  "TraesCS7D02G531200",  "TraesCS7B02G353600",  "TraesCS7D02G442500",  "TraesCS7A02G544900",  "TraesCS4A02G425300",  "TraesCS7A02G063700",  "TraesCS7B02G468300",  "TraesCS7D02G234800",  "TraesCS4A02G486000",  "TraesCS6D02G177400",  "TraesCS7D02G531200",  "TraesCS6A02G254500",  "TraesCS7A02G544900",  "TraesCS7B02G468200",
  "TraesCS2A02G314800",  "TraesCS7A02G207100",  "TraesCS5D02G110600",  "TraesCS7B02G159600",  "TraesCS7A02G190600",  "TraesCS2B02G536100",  "TraesCS2B02G396800",  "TraesCS6D02G389500",  "TraesCS6A02G222100",  "TraesCS4A02G016400",  "TraesCS2D02G206300",  "TraesCS2A02G252700",  "TraesCS3D02G003600",  "TraesCS7B02G095500",  "TraesCS7A02G190600",  "TraesCS7B02G063400",  "TraesCS4A02G036800",
  "TraesCS4B02G269800",  "TraesCS4D02G269100",  "TraesCS6A02G238000",  "TraesCS7D02G019600",  "TraesCS1A02G129300",  "TraesCS6D02G157200",  "TraesCS3A02G186500",  "TraesCS3B02G216000",  "TraesCS6A02G171800",  "TraesCS2D02G583500",  "TraesCS3B02G277500",  "TraesCS3D02G248700",  "TraesCS4A02G388600",  "TraesCS4B02G162700",  "TraesCS4D02G166800",  "TraesCS5A02G405800",  "TraesCS5D02G157300",  "TraesCS5B02G486000",  "TraesCS4A02G107400",
  "TraesCS7D02G236300",  "TraesCS2D02G287300",  "TraesCS2B02G305900",  "TraesCS5A02G022200",  "TraesCS5A02G230500",  "TraesCS5D02G237300",  "TraesCS4D02G316800",  "TraesCS4B02G320300",  "TraesCS4B02G217400",  "TraesCS5B02G229000",  "TraesCS5D02G152400",  "TraesCS5A02G138700",  "TraesCS5B02G138300",  "TraesCS2B02G152000",  "TraesCS6A02G381900",  "TraesCS6D02G366600",
  "TraesCS7D02G181600",  "TraesCS3D02G395000",  "TraesCS3B02G039000",  "TraesCS3A02G400600",  "TraesCS3D02G036100",  "TraesCS3A02G044800",  "TraesCS6B02G183500",  "TraesCS7B02G385700",  "TraesCS4B02G330800",  "TraesCS7D02G263800",  "TraesCS1D02G176000",  "TraesCS3B02G209600",  "TraesCS3A02G179900",  "TraesCS4D02G190900",  "TraesCS4A02G114700",  "TraesCS3B02G209800",
  "TraesCS4B02G189500",  "TraesCS3A02G180100",  "TraesCS3D02G184800",  "TraesCS3B02G209700",  "TraesCS4D02G166800",  "TraesCS2B02G534500",  "TraesCS2D02G506900",  "TraesCS7A02G349000",  "TraesCS6A02G210900",  "TraesCS6D02G194000",  "TraesCS7A02G391200",  "TraesCS2B02G378100",  "TraesCS6B02G240000",  "TraesCS2D02G357600",  "TraesCS6D02G205300",  "TraesCS2D02G285300",  "TraesCS2A02G286700",
  "TraesCS2B02G303600",  "TraesCS3A02G218000",  "TraesCS1A02G122500",  "TraesCS1D02G218000",  "TraesCS1A02G215500",  "TraesCS1B02G229000",  "TraesCS6A02G260400",  "TraesCS2B02G118100",  "TraesCS1D02G082600",  "TraesCS7A02G357800",  "TraesCS4A02G440700",  "TraesCS4B02G210900",  "TraesCS4A02G093500",  "TraesCS4D02G211600",  "TraesCS7D02G084100")
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  #labelMarkers = markerGenes,
  returnMatrix = TRUE,
  transpose = TRUE,
  clusterCols = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap0.7", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)

#Visualizing Marker Genes on an Embedding
p <- plotEmbedding(
  ArchRProj = projHeme2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-WO-Imputation0.7.pdf", 
        ArchRProj = projHeme2, 
        addDOC = FALSE, width = 5, height = 5)

#Marker Genes Imputation with MAGIC
projHeme2 <- addImputeWeights(projHeme2, reducedDims = "IterativeLSI4")
p <- plotEmbedding(
  ArchRProj = projHeme2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme2)
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = projHeme2, 
        addDOC = FALSE, width = 5, height = 5)
#Track Plotting with ArchRBrowser
p <- plotBrowserTrack(
  ArchRProj = projHeme2, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$CD14)
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = projHeme2, 
        addDOC = FALSE, width = 5, height = 5)
#Launching the ArchRBrowser
#ArchRBrowser(projHeme2)

#save object
saveRDS(projHeme2,file="projHeme2.rds")
setwd("/public/workspace/xuyongyue/ZLH/scATAC-seq/newest")
pdf("test1.pdf")
projHeme2<-readRDS(file="projHeme2.rds")

#Defining Cluster Identity with scRNA-seq
#Cross-platform linkage of scATAC-seq cells with scRNA-seq cells
seRNA<-readRDS(file="/public/workspace/xuyongyue/ZLH/ajbk/pmbc.final.rds")
seRNA
DimPlot(seRNA, reduction = 'umap', label=TRUE, pt.size = 2.0)
colnames(colData(seRNA))
table(colData(seRNA)$BioClassification)
#Unconstrained Integration
projHeme2 <- addGeneIntegrationMatrix(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI4",
  seRNA = seRNA,
  addToArrow = TRUE,
  groupRNA = "seurat_clusters",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  force = TRUE
)
#Constrained Integration
cM <- as.matrix(confusionMatrix(projHeme2$Clusters, projHeme2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
#preClust      
[1,] "4"      "C10"
[2,] "1"      "C11"
[3,] "4"      "C6" 
[4,] "13"     "C9" 
[5,] "7"      "C14"
[6,] "0"      "C12"
[7,] "17"     "C7" 
[8,] "11"     "C8" 
[9,] "6"      "C15"
[10,] "15"     "C13"
unique(unique(projHeme2$predictedGroup_Un))
#From scRNA
cTNK <- "4"
cTNK
#cNonTNK <- paste0(c(0, 1, 4, 6, 10:11, 13, 15, 17, 20:21), collapse="|")
cNonTNK <- paste0(c(0, 1, 6, 7, 11, 13, 15, 17), collapse="|")
cNonTNK
#Assign scATAC to these categories
clustTNK <- rownames(cM)[grep(paste("^",cTNK,"$", sep=""), preClust)]
clustTNK
clustNonTNK <- rownames(cM)[grep(paste("^",cNonTNK,"$", sep=""), preClust)]
clustNonTNK
#RNA get cells in these categories
seRNA<- as.SingleCellExperiment(seRNA)
rnaTNK <- colnames(seRNA)[grep(cTNK, colData(seRNA)$seurat_clusters)]
head(rnaTNK)
rnaNonTNK <- colnames(seRNA)[grep(cNonTNK, colData(seRNA)$seurat_clusters)]
head(rnaNonTNK)
groupList <- SimpleList(
  TNK = SimpleList(
    ATAC = projHeme2$cellNames[projHeme2$Clusters %in% clustTNK],
    RNA = rnaTNK
  ),
  NonTNK = SimpleList(
    ATAC = projHeme2$cellNames[projHeme2$Clusters %in% clustNonTNK],
    RNA = rnaNonTNK
  )    
)
projHeme2 <- addGeneIntegrationMatrix(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI4",
  seRNA = seRNA,
  addToArrow = FALSE, 
  groupList = groupList,
  groupRNA = "seurat_clusters",
  nameCell = "predictedCell_Co",
  nameGroup = "predictedGroup_Co",
  nameScore = "predictedScore_Co"
)
#Comparing Unconstrained and Constrained Integrations
pal <- paletteDiscrete(values = colData(seRNA)$seurat_clusters)
pal
p1 <- plotEmbedding(
  projHeme2, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un", 
  pal = pal
)
p2 <- plotEmbedding(
  projHeme2, 
  colorBy = "cellColData", 
  name = "predictedGroup_Co", 
  pal = pal
)
plotPDF(p1,p2, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme2", load = FALSE)

#Adding Pseudo-scRNA-seq profiles for each scATAC-seq cell
projHeme3 <- addGeneIntegrationMatrix(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI4",
  seRNA = seRNA,
  addToArrow = TRUE,
  force= TRUE,
  groupList = groupList,
  groupRNA = "seurat_clusters",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)
getAvailableMatrices(projHeme3)
projHeme3 <- addImputeWeights(projHeme3,reducedDims = "IterativeLSI4")
markerGenes  <- c("TraesCS2D02G349500", "TraesCS7B02G468200", "TraesCS7B02G159600", "TraesCS4A02G016400", "TraesCS2B02G536100", "TraesCS6A02G222100", 
                  "TraesCS4D02G232000", "TraesCS4A02G074700", "TraesCS4A02G036800", "TraesCS6D02G157200", "TraesCS3A02G186500", "TraesCS5D02G425200", 
                  "TraesCS4A02G388600", "TraesCS4B02G162700", "TraesCS5A02G405800", "TraesCS5D02G157300", "TraesCS7D02G236300", "TraesCS4D02G316800", 
                  "TraesCS5A02G230500", "TraesCS2B02G152000", "TraesCS3D02G395000", "TraesCS6B02G183500", "TraesCS3B02G209600", 
                  "TraesCS4D02G166800", "TraesCS2B02G378100", "TraesCS6A02G210900", "TraesCS2A02G286700", "TraesCS1A02G215500", "TraesCS4A02G093500")

p1 <- plotEmbedding(
  ArchRProj = projHeme3, 
  colorBy = "GeneIntegrationMatrix", 
  name = markerGenes, 
  continuousSet = "horizonExtra",
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme3)
)
p2 <- plotEmbedding(
  ArchRProj = projHeme3, 
  colorBy = "GeneScoreMatrix", 
  continuousSet = "horizonExtra",
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme3)
)
p1c <- lapply(p1, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

p2c <- lapply(p2, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

do.call(cowplot::plot_grid, c(list(ncol = 3), p1c))

do.call(cowplot::plot_grid, c(list(ncol = 3), p2c))

plotPDF(plotList = p1, 
        name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = projHeme3, 
        addDOC = FALSE, width = 5, height = 5)

#Labeling scATAC-seq clusters with scRNA-seq information
cM <- confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup)
labelOld <- rownames(cM)
labelOld
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew
remapClust <- c(
  "0" = "epidermis/cortex I",
  "1" = "immature pericycle cells (IPC)",
  "2" = "meristem I",
  "3" = "proximal meristem",
  "4" = "xylem pole pericycle (XPP)",
  "5" = "provascular cell",
  "6" = "root hair",
  "7" = "meristem ",
  "8" = "epidermis/root hair",
  "9" = "phloem pole pericycle (PPP)",
  "10" = "endodermis I (casparian strip)",
  "11" = "metaxylem",
  "12" = "epidermis/cortex II",
  "13" = "proxylem",
  "14" = "immature sieve elements",
  "15" = "root cap",
  "16" = "prophloem",
  "17" = "endodermis II (casparian strip)",
  "18" = "stem cell niche (SCN)",
  "19" = "root border cell",
  "20" = "columella",
  "21" = "companion cell")
remapClust <- remapClust[names(remapClust) %in% labelNew]
labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)
labelNew2
projHeme3$Clusters2 <- mapLabels(projHeme3$Clusters, newLabels = labelNew2, oldLabels = labelOld)
p1 <- plotEmbedding(projHeme3, colorBy = "cellColData", name = "Clusters2")
p1
plotPDF(p1, name = "Plot-UMAP-Remap-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = projHeme3, outputDirectory = "Save-ProjHeme3", load = FALSE)

#Pseudo-bulk Replicates in ArchR
#Making Pseudo-bulk Replicates
saveRDS(projHeme3,"projHeme3.rds")
projHeme3 <- readRDS("projHeme3.rds")
#projHeme3 <- readRDS(file = "/public/workspace/xuyongyue/ZLH/scATAC-seq/newest/Save-ProjHeme3/Save-ArchR-Project.rds"
projHeme4 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "Clusters2")

#Calling Peaks with ArchR
# Calling Peaks w/ Macs2
pathToMacs2 <- findMacs2()
projHeme4 <- addReproduciblePeakSet(
  ArchRProj = projHeme4, 
  groupBy = "Clusters2", 
  pathToMacs2 = pathToMacs2,
  genomeSize = 17e9
)
getPeakSet(projHeme4)
#Add Peak Matrix
saveArchRProject(ArchRProj = projHeme4, outputDirectory = "Save-ProjHeme4", load = FALSE)
#projHeme4 <- readRDS(file="/public/workspace/xuyongyue/ZLH/scATAC-seq/new/Save-ProjHeme4/Save-ArchR-Project.rds")
projHeme5 <- addPeakMatrix(projHeme4)
getAvailableMatrices(projHeme5)
saveRDS(projHeme5,"projHeme5.rds")

#Identifying Marker Peaks with ArchR
#Identifying Marker Peaks with ArchR
#Our scRNA labels
table(projHeme5$Clusters2)
markersPeaks <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerList
#access this from the list via the $ accessor.
#markerList$companion_cell
colnames(markersPeaks) <- c("cortex","endodermisII","epidermis","ground_meristem",
                            "immature_sieve_elements","metaxylem","pericycle","proxylem","root_hair",
                            "vascular_parenchyma_initial_cell")
#Plotting Marker Peaks in ArchR
# Marker Peak Heatmaps
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  #returnMatrix = TRUE,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",  #默认参数是FDR<=0.1
  transpose = TRUE
)

#Marker Peak MA and Volcano Plots
pma <- plotMarkers(seMarker = markersPeaks, name = "cortex", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pma
pv <- plotMarkers(seMarker = markersPeaks, name = "cortex", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
pv
plotPDF(pma, pv, name = "cortex-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
#Marker Peaks in Browser Tracks
p <- plotBrowserTrack(
  ArchRProj = projHeme5, 
  groupBy = "Clusters2", 
  geneSymbol = c("TraesCS2D02G349500"),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["proximal_meristem"],
  upstream = 50000,
  downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$TraesCS7B02G095500)
plotPDF(p, name = "TraesCS7B02G095500_Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
#Pairwise Testing Between Groups
#the peaks that are higher in the group passed to useGroups will have positive fold change values 
#while the peaks that are higher in the group passed to bgdGroups will have negative fold change values
markerTest <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "cortex Ⅰ (parenchyma",
  bgdGroups = "epidermis/LRC"
)
pma <- plotMarkers(seMarker = markerTest, name = "cortex Ⅰ (parenchyma", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma
pv <- plotMarkers(seMarker = markerTest, name = "epidermis/LRC", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
pv
plotPDF(pma, pv, name = "cortex-vs-epidermis-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

#Motif and Feature Enrichment with ArchR
#Motif Enrichment in Differential Peaks
library("JASPAR2018")
projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "jaspar2018", name = "Motif")
motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 1.5"
)
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
ggUp
motifsDo <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)
df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
plotPDF(ggUp, ggDo, name = "cortex-vs-epidermis-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE, max.overlaps = 1000000)
#Motif Enrichment in Marker Peaks
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
#ustom Enrichment
EncodePeaks <- c(
  Encode_K562_GATA1 = "https://www.encodeproject.org/files/ENCFF632NQI/@@download/ENCFF632NQI.bed.gz",
  Encode_GM12878_CEBPB = "https://www.encodeproject.org/files/ENCFF761MGJ/@@download/ENCFF761MGJ.bed.gz",
  Encode_K562_Ebf1 = "https://www.encodeproject.org/files/ENCFF868VSY/@@download/ENCFF868VSY.bed.gz",
  Encode_K562_Pax5 = "https://www.encodeproject.org/files/ENCFF339KUO/@@download/ENCFF339KUO.bed.gz"
)
addPeakSet(
  ArchRProj = NULL,
  peakSet = NULL,
  genomeAnnotation = getGenomeAnnotation(ArchRProj),
  force = FALSE
)

# Motif Deviations
# lets make sure we have added motif annotations to our ArchRProject
if("Motif" %ni% names(projHeme5@peakAnnotation)){
  projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")
}
projHeme5 <- addBgdPeaks(projHeme5)
projHeme5 <- addDeviationsMatrix(
  ArchRProj = projHeme5, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs
p <- plotGroups(ArchRProj = projHeme5, 
                groupBy = "Clusters2", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(projHeme5)
)
p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
p <- plotEmbedding(
  ArchRProj = projHeme5, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme5)
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
markerRNA <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA <- markerRNA[markerRNA %ni% c("SREBF1","CEBPA-DT")]
markerRNA
p <- plotEmbedding(
  ArchRProj = projHeme5, 
  colorBy = "GeneScoreMatrix", 
  name = sort(markerRNA), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme5)
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
markerRNA <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "GeneIntegrationMatrix")
markerRNA <- markerRNA[markerRNA %ni% c("SREBF1","CEBPA-DT")]
markerRNA
p <- plotEmbedding(
  ArchRProj = projHeme5, 
  colorBy = "GeneIntegrationMatrix", 
  name = sort(markerRNA), 
  embedding = "UMAP",
  continuousSet = "blueYellow",
  imputeWeights = getImputeWeights(projHeme5)
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

