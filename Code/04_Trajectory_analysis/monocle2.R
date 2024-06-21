

#### ===============================================================
#### monocle2
#### ===============================================================
#### Load Packages
library(monocle)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggsci)
library(harmony)
library(ggplot2)
library(ggsci)
options(mc.cores = 30)
cores=30

#### 1. Load Data
cell='Spinous cells'
setwd("~/dat/Epithelial/07remun/Spinous")
data <- readRDS('celltype.rds')


#### 2. Build CDS objects
expr_matrix <- as(as.matrix(data@assays$RNA@counts), 'sparseMatrix')
p_data <- data@meta.data 
f_data <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                    phenoData = pd,
                    featureData = fd,
                    lowerDetectionLimit = 0.5,
                    expressionFamily = negbinomial.size())


#### 3. Quality Control
#### Estimate Size Factor and Dispersion
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#### Filter Low-quality Cells
cds <- detectGenes(cds, min_expr = 0.1) 
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10)) 


#### 4. Identify DEGs
diff <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr="~seurat_clusters",cores=cores) 
deg <- subset(diff, qval < 0.01) 
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)
write.table(deg,file=paste0(cell, ".monocle.DEG.xls"),col.names=T,row.names=F,sep="\t",quote=F)


#### 5. Trajectory Construction Gene Visualization
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()

#### 6. Reduce Dimension
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

#### 7. Constructing a Quasi-timeline Trajectory and Arranging Cells within a Quasi-timeline
cds <- orderCells(cds)

Time_diff <- differentialGeneTest(cds[ordergene, ], cores = cores, fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)]
write.csv(Time_diff, "Time_diff_all.csv", row.names = F)

#### 8. Save
saveRDS(cds, paste0(cell, "_endmtcells.rds"))

#### 9. Plot
p1 <- plot_cell_trajectory(cds, color_by = "seurat_clusters")  + scale_color_npg() 
p2 <- plot_cell_trajectory(cds, color_by = "cell_type")  + scale_color_npg() 
p3 <- plot_cell_trajectory(cds, color_by = "State")  + scale_color_nejm()
p4 <- plot_cell_trajectory(cds, color_by = "Pseudotime") 
plotc <- p1|p2|p3|p4
ggsave("TI_pseudotimeplot.pdf", plot = plotc, width = 16, height = 8)

p <- plot_cell_trajectory(cds, color_by = "seurat_clusters")  + scale_color_npg() +facet_wrap("~seurat_clusters", nrow=1)
ggsave("TI_facet.pdf", p, width = 16, height = 8)

#### Changes in Gene Expression with Cellular Status
keygenes <- head(ordergene,10)
cds_subset <- cds[keygenes,]
p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters")
p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "cell_type")
p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
p4 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
plotc <- p1|p2|p3|p4
ggsave("Genes_pseudotimeplot.pdf", plot = plotc, width = 16, height = 8)#ָ������

p1 <- plot_genes_jitter(cds[keygenes,], grouping = "State", color_by = "State")   # seurat_clusters or cell_type or State or Pseudotime ...
p2 <- plot_genes_violin(cds[keygenes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(cds[keygenes,], color_by = "State")
plotc <- p1|p2|p3
ggsave("Genes_Jitterplot.pdf", plot = plotc, width = 16, height = 20)

#### Heatmap
Time_genes <- Time_diff[,1]
p <- plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)
hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(Time_diff_sig, "Time_diff_sig.csv", row.names = F)

