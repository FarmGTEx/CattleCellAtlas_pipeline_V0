##Reference website: https://satijalab.org/seurat/articles/get_started.html

library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(ddqcr)
library(doubletFinder)
library(SingleR)
library(celldex)
library(ggplot2)
library(reshape2)
ref = HumanPrimaryCellAtlasData()
setwd("{your_workspace}")

###read the count matrix
sc.data<-Read10X(data.dir="./sample_id/outs/filtered_feature_bc_matrix")
sc <- CreateSeuratObject(counts = sc.data, project = "sc", min.cells = 3)

###remove low quality cell
##remove cells with low nFeature, nCount and high mito_genes based on the ddqc method
mt_gene<-read.table(file="./bta_genename_final.txt",header=F,sep="\t",quote="",row.names=1)  #change with different species
mt_select_gene<-intersect(rownames(mt_gene),rownames(sc))
sc[["percent.mt"]] <-PercentageFeatureSet(sc, features = mt_select_gene)
sc <- initialQC(sc)
sc[["percent.mt"]] <-PercentageFeatureSet(sc, features = mt_select_gene)
df.qc <- ddqc.metrics(sc)
sc <- filterData(sc, df.qc)
#sc<- subset(sc, subset = percent.mt < 10) #depends on your data quality
##remove doublets
#pre-process
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(sc), 10)
all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)
sc <- RunPCA(sc, features = VariableFeatures(object = sc))
sc <- FindNeighbors(sc, dims = 1:10)
sc <- FindClusters(sc)
sc <- RunUMAP(sc, dims = 1:10)
#calculate the percentage of heterologous doublets
sweep.data <- paramSweep_v3(sc, PCs=1:10)
sweep.stats <- summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)
#calculate the percentage of homotypic doublets
homotypic.prop=modelHomotypic(sc@meta.data$seurat_clusters)
nExp_poi=round((ncol(sc)*8*1e-6)*length(sc$orig.ident))
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
#remove doublets
sc=doubletFinder_v3(sc, PCs = 1:10, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])),nExp = nExp_poi.adj, reuse.pANN = FALSE)
sc@meta.data$DF_hi.lo<- sc@meta.data[,9] 
Doublet<-table(sc@meta.data$DF_hi.lo=="Doublet")
Idents(sc) <- "DF_hi.lo"
sc <-subset(x = sc, idents="Singlet")
saveRDS(sc,file="./sample.rds")

###pre-process of cell clustering
##read different rds files in a given tissue (Suppose there are multiple samples in each tissue)
sc1 <- readRDS("./sample1.rds")
sc2 <- readRDS("./sample2.rds")
##merge different samples
sc1$orig.ident<-"sc1"
sc2$orig.ident<-"sc2"
sc <- merge(sc1, y = sc2, add.cell.ids = c("sc1", "sc2"), project = "sc")

###Normalizing the merged data
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)

###Feature selection
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
##Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sc), 10)
##plot variable features with and without labels
plot1 <- VariableFeaturePlot(sc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

###Scailing data
all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)

###linear dimensional reduction
sc <- RunPCA(sc, features = VariableFeatures(object = sc))
DimHeatmap(sc, dims = 1:15, cells = 500, balanced = TRUE)

###remove batch effect
sc <- RunHarmony(object = sc, group.by.vars=c("orig.ident") , plot_convergence = FALSE) #depend on your batch effects

###Determine the dimensionality
sc <- JackStraw(sc, num.replicate = 100, dims = 30)
sc <- ScoreJackStraw(sc, dim = 1:30)
Plot1 = ElbowPlot(sc, ndims = 30)
Plot2 = JackStrawPlot(sc, dim = 1:30)
Plot1 + Plot2
sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30)  #select suitable nPCs(dims) from Plot1 + Plot2

###cell cluster
sc <- FindClusters(sc, resolution = 0.5)
sc <- RunUMAP(sc, reduction = "harmony", dims = 1:30)  #follow your selected nPCs
DimPlot(sc, reduction = "umap", group.by = "orig.ident" , pt.size=1, label = FALSE , label.size = 3)
DimPlot(sc, reduction = "umap", label = TRUE)
DimPlot(sc, label = TRUE, split.by = "orig.ident")  + NoLegend()
saveRDS(sc, file= "./sample.rds"))

###Marker genes extraction
sc.markers <- FindAllMarkers(sc, min.pct = 0.25, logfc.threshold = 0.25)
sc.markers %>%
group_by(cluster) %>%
slice_max(n = 2, order_by = avg_log2FC)
write.csv(sc.markers, file= "./sample_marker.csv")

###Auto-annotation by singleR
norm_count = GetAssayData(sc, slot = "data")
pred<- SingleR(test = norm_count, ref = ref, labels = ref$label.main)
sc@meta.data$labels <-pred$labels
DimPlot(sc, group.by = c("seurat_clusters", "labels"),reduction = "umap")
aa=table(pred$labels,sc$seurat_clusters)
aa=apply(aa,2,function(x) x/sum(x))
df=as.data.frame(melt(aa))
df$Var2=as.factor(df$Var2)
ggplot(df, aes(Var2, Var1)) + geom_point(aes(size = value), colour = "green") + theme_bw()
saveRDS(sc, file= "./sample_auto_annotation.rds"))
