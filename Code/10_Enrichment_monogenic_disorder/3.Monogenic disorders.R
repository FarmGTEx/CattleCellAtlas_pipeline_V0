library(dplyr)
library(Seurat)
library(patchwork)

##extract all celltypes in each tissue
library(readxl)
celltype <- read_excel("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/20230912 Cattle-cell_type_annotation-BO.xlsx",sheet=3)
celltype <- data.frame(celltype)
str <- "/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/"
setwd("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds")
list0 <- list.files(pattern = ".rds")
list0 <- sapply(list0, function(x) unlist(strsplit(x, "\\."))[1])
list0 <- data.frame(list0)
tissue <- list0$list0
list1<-NULL
for(i in tissue){
  list1[[i]]<-paste0("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/", i, ".rds")
}
list2<-NULL
for(i in tissue){
  list2[[i]]<-paste0("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/", i, "_anno.rds")
}

for(i in 1:length(tissue)) {
   sc <- readRDS(list1[[i]])
   tissue_cell <- celltype[celltype$Tissue == tissue[[i]],]
   cell_type <- tissue_cell$Cell.type
   Idents(sc) <- "seurat_clusters"
   table(Idents(sc))
   names(cell_type) <- levels(sc)
   sc <- RenameIdents(sc, cell_type)
   sc@meta.data$CellType <- Idents(sc)
   saveRDS(sc, file = list2[[i]])
}


#count_CPM
library(edgeR)
str <- "/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/"
setwd("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds")
list0 <- list.files(pattern = ".rds")
list0 <- sapply(list0, function(x) unlist(strsplit(x, "\\."))[1])
list0 <- data.frame(list0)
list0 <- list0$list0
#list0 <- c("Sublingual gland","Ovary","Oviduct","Placenta","Lung","PBMC","Spleen","Esophage","Duodenum","Ileum","Colon","Trachea")
list1<-NULL
for(i in list0){
  list1[[i]]<-paste0("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/", i, ".rds")
}
for (i in 1:length(list0)) { 
    sc <- readRDS(list1[[i]])
    Idents(sc) <- "CellType"
    celltype <- unique(sc$CellType)
    for(j in 1:length(celltype)) {
        sc1 <- subset(x=sc, idents = celltype[[j]]) 
        dat <- data.frame(sc1@assays$RNA@counts)
        zero_count <- apply(dat, 1, function(row) sum(row == 0))
        zero_percentage <- zero_count / ncol(dat)
        dat <- dat[zero_percentage <= 0.5, ]
        cpm <- edgeR::cpm(dat)
        mean_cpm <- rowMeans(cpm)
        write.csv(mean_cpm, file = paste0(str, list0[[i]], "_", celltype[[j]], ".csv"))
    }
}

#combine the same celltype from different tissues
library(readxl)
all_celltype <- read_excel("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CellType_group.xlsx", sheet = 1)
all_celltype <- data.frame(all_celltype)
list0 <- all_celltype$Celltype
list1 <- paste0("anno_",list0)
for (i in 1:length(list0)) { 
    file_list <- list.files(path = "/faststorage/project/cattle_gtexs/Global_atlas/Global_cell_atlas/Cell_cycle", pattern = list1[[i]], full.names = TRUE)   
    merged_data <- data.frame()
    for (file in file_list) {
        data <- read.csv(file)
        if (nrow(merged_data) == 0) {
           merged_data <- data
           } else {
           merged_data <- merge(merged_data, data, by = "X", all = TRUE)
        }
    }
    merged_data[is.na(merged_data)]=0
    write.csv(merged_data, file = paste0("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/gene_filter/test/20/combine/", list0[[i]], ".csv"), row.names = FALSE)
    #}
}

#calculate the mean tpm and combine different celltype together   
str <- "/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/gene_filter/test/20/combine/"  
setwd("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/gene_filter/test/20/combine")
file <- list.files(pattern = ".csv")
list0 <- sapply(file, function(x) unlist(strsplit(x, "\\."))[1])
list0 <- data.frame(list0)
list0 <- list0$list0
list1 <- list.files(pattern = ".csv")

merged_data <- read.csv(list1[[1]])
colnames(merged_data)[2] <- list0[[1]]
colnames(merged_data)[1] <- "gene"
for (i in 2:length(list0)) { 
    dat <- read.csv(list1[[i]])
    if (ncol(dat) >2) {
       rownames(dat) <- dat$X
       dat <- dat[,-1]
       dat <- rowMeans(dat) 
       dat <- data.frame(dat)
       colnames(dat) <- list0[[i]]
       dat$gene <- rownames(dat)
       } else {
       colnames(dat)[2] <- list0[[i]]
       colnames(dat)[1] <- "gene"
    }
    merged_data <- merge(merged_data, dat, by = "gene", all = TRUE)
}
merged_data[is.na(merged_data)]=0
rownames(merged_data) <- merged_data$gene
merged_data <- merged_data[,-1]
write.csv(merged_data, file = "/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/gene_filter/test/20/combine_filter20_genes_celltype_cpm_mean.csv")

##tau test
library(readxl)
cpm <- read.csv("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/all_genes/combine_all_genes_celltype_cpm_mean.csv")
rownames(cpm) <- cpm$gene
cpm <- cpm[,-1:-2]
group <- read_excel("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CellType_group.xlsx",sheet = 1)
group <- data.frame(group)
rownames <- group$Celltype
group <- group[,-1]
group <- data.frame(group)
rownames(group) <- rownames
tau <- read.csv("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/gene_filter/test/20/output/tau_specificity_filter20.csv")
tau_select <- tau[tau$X0>0.95,]
rownames(tau_select) <- tau_select$gene
cpm_select <- cpm[rownames(cpm) %in% rownames(tau_select),]    
dat=log2(cpm_select+1)    
dat_scale=t(scale(t(dat)))
tiff(filename = '/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/gene_filter/test/20/output/tau_filter20.tiff',width =3000,height=3000,res=300)
pheatmap(dat_scale, cluster_rows = TRUE, cluster_cols = TRUE, annotation_col = group, color = colorRampPalette(colors = c("blue", "white","red"))(100), show_rownames = FALSE, show_colnames = FALSE)
dev.off()

##calculate the correlation between cell type groups and disorder groups
#extract high average zscore genes from each cell type group
library(readxl)
cpm <- read.csv("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global_atlas/All_rds/annotation_rds/CellType/CPM/all_genes/combine_all_genes_celltype_cpm_mean.csv")
rownames(cpm) <- cpm$gene
cpm <- cpm[,-1:-2]
zscore <- read.csv("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global_atlas/All_rds/annotation_rds/CellType/CPM/all_genes/output/zscore_specificity.csv")
rownames(zscore) <- zscore$gene
zscore <- zscore[,-1]
disorder <- read.csv("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global_atlas/All_rds/annotation_rds/CellType/omia_genes.csv")
disease <- unique(disorder$Disease_Group) 
group <- read_excel("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global_atlas/All_rds/annotation_rds/CellType/CellType_group.xlsx", sheet = 1)
group <- data.frame(group)
list0 <- unique(group$Group)
for (i in 1:length(list0)){
    celltype <- group[group$Group == list0[[i]],]
    celltype <- celltype$Celltype
    merged_zscore <- data.frame()
    for(j in 1:length(celltype)){  
        a = which(names(zscore)==celltype[[j]])
        select_zscore <- zscore[zscore[,a] >0.75, ] 
        if (nrow(merged_zscore) == 0) {
           merged_zscore <- select_zscore
           } else {
           merged_zscore <- rbind(merged_zscore, select_zscore)
        } 
    } 
    write.csv(merged_zscore, file = paste0("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global_atlas/All_rds/annotation_rds/CellType/CPM/all_genes/output/select_zscore/", list0[[i]], ".csv"))
} 

#fisher test between cell type group and disorder groups
str1 <- "/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global_atlas/All_rds/annotation_rds/CellType/CPM/all_genes/output/select_zscore/"
setwd("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global_atlas/All_rds/annotation_rds/CellType/CPM/all_genes/output/select_zscore")
file <- list.files(pattern = ".csv")
list0 <- sapply(file, function(x) unlist(strsplit(x, "\\."))[1])
list0 <- data.frame(list0)
list0 <- list0$list0    
list1<-NULL
for(i in list0){
  list1[[i]]<-paste0(str1, i, ".csv")
}    
p_value <- array(data = NA,dim = c(7,10))
for(i in 1:length(list1)) {
    zscore <- read.csv(list1[[i]])
    zscore <- na.omit(zscore)
    rownames(zscore) <- zscore$X
    zscore <- zscore[,-1]
    zscore_genes <- rownames(zscore)
    for(j in 1:length(disease)) {
        selet_disorder <- disorder[disorder$Disease_Group == disease[j],]
        disorder_genes <- selet_disorder$Gene
        same_genes <- intersect(zscore_genes,disorder_genes)
        a = length(same_genes)
        b = length(disorder_genes)
        c = length(zscore_genes)
        d = 23945
        if (((b*c)/d) >= 5) {
            chisq_result <- chisq.test(matrix(c(a,b-a,c-a,d-b-c+a),nrow=2)) 
            p <- chisq_result$p.value
            } else if (((b*c)/d) >= 1) {
                chisq_result <- chisq.test(matrix(c(a,b-a,c-a,d-b-c+a),nrow=2), correct=TRUE)
                p <- chisq_result$p.value
                } else {
                  fisher_result <- chisq.test(matrix(c(a,b-a,c-a,d-b-c+a),nrow=2))
                  p <- fisher_result$p.value
                  }
            p_value[i,j] <- p
        }
}
colnames(p_value) <- disease 
rownames(p_value) <- list0

#OR
OR_value <- array(data = NA,dim = c(7,10))
for(i in 1:length(list1)) {
    zscore <- read.csv(list1[[i]])
    zscore <- na.omit(zscore)
    rownames(zscore) <- zscore$X
    zscore <- zscore[,-1]
    zscore_genes <- rownames(zscore)
    for(j in 1:length(disease)) {
        selet_disorder <- disorder[disorder$Disease_Group == disease[j],]
        disorder_genes <- selet_disorder$Gene
        same_genes <- intersect(zscore_genes,disorder_genes)
        a = length(same_genes)
        b = length(disorder_genes)
        c = length(zscore_genes)
        d = 23945 
        OR <- (a/(b-a))/((c-a)/(d-b-c+a))
        OR_value[i,j] <- OR
        }
}
colnames(OR_value) <- disease 
rownames(OR_value) <- list0

##FDR correction
list0 <- colnames(p_value) 
merged_FDR <- data.frame()   
for(i in 1:10) {  
    FDR <- p.adjust(p_value[,i], method ="BH")  
    FDR <- data.frame(FDR)
    colnames(FDR) <- list0[[i]]
    if (ncol(merged_FDR) == 0) {
           merged_FDR <- FDR
           } else {
           merged_FDR <- cbind(merged_FDR, FDR)
        }
}

library(reshape2)
OR_value <- data.frame(OR_value)
OR_value <- cbind(RowNames = rownames(OR_value), OR_value)
plot <- melt(OR_value, id.vars="RowNames")
names(plot) <- c("CellType", "Disorder", "OR")

write.csv(p_value, file = "/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global_atlas/All_rds/annotation_rds/CellType/CPM/all_genes/output/p_value.csv")
write.csv(merged_FDR, file = "/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global_atlas/All_rds/annotation_rds/CellType/CPM/all_genes/output/FDR.csv")
write.csv(OR_value, file = "/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global_atlas/All_rds/annotation_rds/CellType/CPM/all_genes/output/OR.csv")
write.csv(plot, file = "/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global_atlas/All_rds/annotation_rds/CellType/CPM/all_genes/output/plot.csv")

#plot
library(ggplot2)
dat <- read.csv("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/all_genes/output/plot.csv")
dat$Significance <- factor(dat$Significance)
dat$log2OR <- log2(dat$OR + 1)
tiff(filename = '/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/all_genes/output/cellgroup_disorder.tiff',width =3000,height=2000,res=300)
ggplot() + 
geom_point(data = dat,aes(x = Disorder, y = CellType, color=log2OR, size=Significance)) +
theme_classic() +
ylab("Cell type groups")+
xlab("Disorder groups")+
theme(axis.title=element_text(size = 18,face="bold",colour = "black"),
      axis.text.x = element_text(angle=45, hjust=1, vjust=1),
      axis.text=element_text(size = 15,face="bold",colour = "black"),  
      legend.title = element_text(face="bold",size=14),
      legend.text = element_text(face="bold",size=12))+
      scale_size_discrete(name = "Significance", breaks = c("0", "1", "2"), labels = c("NS/NA","Nominal significance","Bonferroni corrected")) +
      #scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33')) +
      scale_color_gradientn(values = seq(0,1,0.2),colours = c('gray','red'))
      #theme(text = element_text("serif"))
dev.off()

#select celltype and recalculate zscore
library(readxl)
cpm <- read.csv("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/all_genes/combine_all_genes_celltype_cpm_mean.csv")
rownames(cpm) <- cpm$gene
cpm <- cpm[,-1:-2]
group <- read_excel("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CellType_group.xlsx", sheet = 1)
group <- data.frame(group)
celltype_group <- unique(group$Group)
select_celltype <- group[group$Group == celltype_group[[4]],]
select_celltype <- select_celltype$Celltype
select_cpm <- cpm[, colnames(cpm) %in% select_celltype]
write.csv(select_cpm, "/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/all_genes/output/select_celltype/Immune_zscore.csv")
#Then re-calculate the zcsore

##extract interesting cell type and disease
library(readxl)
library(reshape2)
cpm <- read.csv("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/all_genes/combine_all_genes_celltype_cpm_mean.csv")
rownames(cpm) <- cpm$gene
cpm <- cpm[,-1:-2]
group <- read_excel("/faststorage/project/cattle_gtexs/Downstream_analysis/Monogenic_disease/CellType_group.xlsx", sheet = 1)
group <- data.frame(group)
celltype_group <- unique(group$Group)
disorder <- read.csv("/faststorage/project/cattle_gtexs/Downstream_analysis/Monogenic_disease/omia_genes.csv")
disease <- unique(disorder$Disease_Group) 
select_celltype <- group[group$Group == celltype_group[[4]],]
select_celltype <- select_celltype$Celltype
selet_disorder <- disorder[disorder$Disease_Group == disease[3],]
disorder_genes <- selet_disorder$Gene
select_cpm <- cpm[rownames(cpm) %in% disorder_genes, colnames(cpm) %in% select_celltype]
zscore <- read.csv("/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global atlas/All_rds/annotation_rds/CellType/CPM/all_genes/output/select_celltype/output/Epithelial_zscore_specificity.csv")
rownames(zscore) <- zscore$X
zscore <- zscore[,-1]
select_zscore <- zscore[rownames(zscore) %in% disorder_genes, colnames(zscore) %in% select_celltype]
log_cpm <- log10(select_cpm+1)
log_cpm <- cbind(RowNames = rownames(log_cpm), log_cpm)
log_cpm <- melt(log_cpm, id.vars="RowNames")
names(log_cpm) <- c("genes", "CellType", "LogCPM")
select_zscore <- cbind(RowNames = rownames(select_zscore), select_zscore)
select_zscore <- melt(select_zscore, id.vars="RowNames")
names(select_zscore) <- c("genes", "CellType", "zscore")
log_cpm$zscore <- select_zscore$zscore
write.csv(log_cpm, file = "/faststorage/project/cattle_gtexs/Downstream_analysis/Monogenic_disease/select_celltype/Immune_cpm_zscore.csv")
library(ggplot2)
tiff(filename = '/usr/home/ironbank/researchdata/gtex/cattle_scdata/Global_atlas/All_rds/annotation_rds/CellType/CPM/all_genes/output/select_celltype/output/Epithelial.tiff',width =4000,height=2500,res=300)
ggplot() + 
geom_point(data = log_cpm,aes(x = CellType, y = genes, color=zscore, size=LogCPM)) +
theme_classic() +
ylab("Skin disorder genes")+
xlab("Epithelial cell")+
theme(axis.title=element_text(size = 18,face="bold",colour = "black"),
      axis.text.x = element_text(angle=45, hjust=1, vjust=1),
      axis.text=element_text(size = 15,face="bold",colour = "black"),  
      legend.title = element_text(face="bold",size=14),
      legend.text = element_text(face="bold",size=12))+
      scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
#theme(text = element_text("serif"))
dev.off()