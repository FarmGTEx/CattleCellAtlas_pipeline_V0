library(scPagwas)
library(rtracklayer)
library(dplyr)
library(Seurat)
library(patchwork)

##Enrichment between traits and cell types within tissues
#prepare sc data
str <- "~/dat/cattle_scdata/Global atlas/All_rds/annotation_rds/"
setwd("~/dat/cattle_scdata/Global atlas/All_rds/annotation_rds")
list0 <- list.files(pattern = ".rds")
list0 <- sapply(list0, function(x) unlist(strsplit(x, "\\."))[1])
list0 <- data.frame(list0)
tissue <- list0$list0
list1<-NULL
for(i in tissue){
  list1[[i]]<-paste0(str, i, ".rds")
}

#prepare genome annotation data
gtf_df<- rtracklayer::import("~/dat/reference/Bos_taurus.ARS-UCD1.2.110.gtf")
gtf_df <- as.data.frame(gtf_df)
gtf_df <- gtf_df[,c("seqnames","start","end","type","gene_name")]
gtf_df <- gtf_df[gtf_df$type=="gene",]
block_annotation<-gtf_df[,c(1,2,3,5)]
colnames(block_annotation)<-c("chrom", "start","end","label")
block_annotation$chrom <- paste0("chr", block_annotation$chrom)
                
#prepare ld file
ld <- read.table("produc.ld.ld", sep = "", header = T)
lapply(unique(ld$CHR_A), function(i){
  a<-data.table(ld[ld$CHR_A == i,])
  file_name <- paste0("~/dat/complextraits/",i,".Rds")
  saveRDS(a, file = file_name)
})
ld<-lapply(as.character(1:29),function(chrom){
  chrom_ld_file_path <- paste("~/dat/complextraits/ld_chrom", '/', chrom, '.Rds', sep = '')
  ld_data <- readRDS(chrom_ld_file_path)[, .(SNP_A, SNP_B, R)]
  return(ld_data)
})
ld <- setNames(ld, paste0("chr", c(1:29)))

#prepare pathway data
merge_dat <- data.frame()
for (i in 1:length(Genes_by_pathway_kegg)) {
  list1 <- data.frame(Genes_by_pathway_kegg[i])
  list1$group <- names(Genes_by_pathway_kegg[i])
  names(list1)[1] <- "gene"
  merge_dat <- rbind(merge_dat, list1)
}
                
#transfer the Homologous genes by using Biomart
pathway <- read.csv("~/dat/complextraits/scpagwas/cattle_human_pathway.csv")
pathway <- pathway[,-1]
kegg <- split(pathway$gene, pathway$group)

for (k in 1:length(tissue)) {
  sc <- readRDS(list1[[k]])
  sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
  sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(sc)
  sc <- ScaleData(sc, features = all.genes)
  #prepare gwasdata
  gwas <- read.table("~/dat/complextraits/06production/fct/fct_china.mlma" ,sep = "\t")
  gwas <- gwas[,c(1,2,3,8,7,6)]
  names(gwas) <- c("chrom","pos","rsid","se","beta","maf")
  ##separate the pipeline
  Pagwas <- list()
  Pagwas <- Single_data_input(
    Pagwas = Pagwas,
    assay = "RNA",
    Single_data = sc,
    #Pathway_list = Genes_by_pathway_kegg
    Pathway_list = kegg
  )
  sc <- sc[, colnames(Pagwas$data_mat)]
  
  Pagwas <- Pathway_pcascore_run(
    Pagwas = Pagwas,
    #Pathway_list = Genes_by_pathway_kegg
    Pathway_list = kegg
  )

  Pagwas <- GWAS_summary_input(
    Pagwas = Pagwas,
    gwas_data = gwas,
    maf_filter = 0.1
  )

  Pagwas$snp_gene_df <- SnpToGene(
    gwas_data = Pagwas$gwas_data,
    block_annotation = block_annotation,
    marg = 10000
  )

  Pagwas <- Pathway_annotation_input(
    Pagwas = Pagwas,
    block_annotation = block_annotation
  )
        
  Pagwas <- Link_pathway_blocks_gwas(
    Pagwas = Pagwas,
    chrom_ld = ld,
    singlecell = F,
    celltype = T,
    backingpath="~/dat/complextraits/scpagwas/scPagwastest_output")
    Pagwas$lm_results <- Pagwas_perform_regression(Pathway_ld_gwas_data = Pagwas$Pathway_ld_gwas_data)
  Pagwas <- Boot_evaluate(Pagwas, bootstrap_iters = 200, part = 0.5)
  Pagwas$Pathway_ld_gwas_data <- NULL
  results <- data.frame(Pagwas$bootstrap_results)
  write.csv(results, paste0("~/dat/complextraits/scpagwas/out/", "fct_", tissue[[k]], ".csv"))
}

##Enrichment between all celltypes and traits
#prepare sc data
sc <- readRDS("~/dat/cattle_scdata/Global atlas/All_rds/annotation_rds/all_anno.rds")
Idents(sc) <- "CellType"
celltype <- unique(Idents(sc))
sc[["CellName"]] <- colnames(sc)
percent <- 0.1
celltype_cells <- colnames(sc)[sc@meta.data$CellType == celltype[[1]]]
if (length(celltype_cells) > 500) {
num_cells <- round(percent * length(celltype_cells))
} else {
num_cells <- length(celltype_cells)
}
selected_cells = NULL
selected_cells <- c(sample(celltype_cells, num_cells))
for(j in 2:length(celltype)) {
celltype_cells1 <- colnames(sc)[sc@meta.data$CellType == celltype[[j]]]
if (length(celltype_cells1) > 500) {
num_cells1 <- round(percent * length(celltype_cells1))
} else {
num_cells1 <- length(celltype_cells1)
}
selected_cells <- c(sample(selected_cells), sample(celltype_cells1, num_cells1))
}
sc1 <- subset(sc, subset = CellName %in% selected_cells)

merge_sc <- NormalizeData(merge_sc, normalization.method = "LogNormalize", scale.factor = 10000)
merge_sc <- FindVariableFeatures(merge_sc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merge_sc)
merge_sc <- ScaleData(merge_sc, features = all.genes)

##separate the pipeline
Pagwas <- list()
Pagwas <- Single_data_input(
  Pagwas = Pagwas,
  assay = "RNA",
  Single_data = merge_sc,
  #Pathway_list = Genes_by_pathway_kegg
  Pathway_list = kegg
)
merge_sc <- merge_sc[, colnames(Pagwas$data_mat)]
  
Pagwas <- Pathway_pcascore_run(
  Pagwas = Pagwas,
  #Pathway_list = Genes_by_pathway_kegg
  Pathway_list = kegg
)

#prepare gwasdata
str1 <- "~/dat/complextraits/"
group <- list.dirs(str1, full.names = FALSE, recursive = FALSE)
group <- group[c(1:2)]
#merge_dat <- data.frame()
for(i in 1:length(group)) {
str2 <- paste0(str1, group[[i]], "/")
list2 <- list.dirs(str2, full.names = FALSE, recursive = FALSE)
list3<-NULL
for(j in list2){
list3[[j]]<-paste0(str2, j, "/", j, "_china.mlma")
}
for(j in 1:length(list2)) {
gwas <- read.table(list3[[j]] ,sep = "\t")
gwas <- gwas[,c(1,2,3,8,7,6)]
names(gwas) <- c("chrom","pos","rsid","se","beta","maf")

Pagwas <- GWAS_summary_input(
  Pagwas = Pagwas,
  gwas_data = gwas,
  maf_filter = 0.1
)

Pagwas$snp_gene_df <- SnpToGene(
  gwas_data = Pagwas$gwas_data,
  block_annotation = block_annotation,
  marg = 10000
)

Pagwas <- Pathway_annotation_input(
  Pagwas = Pagwas,
  block_annotation = block_annotation
)
        
Pagwas <- Link_pathway_blocks_gwas(
  Pagwas = Pagwas,
  chrom_ld = ld,
  singlecell = F,
  celltype = T,
backingpath="~/dat/complextraits/scpagwas/scPagwastest_output")
Pagwas$lm_results <- Pagwas_perform_regression(Pathway_ld_gwas_data = Pagwas$Pathway_ld_gwas_data)
Pagwas <- Boot_evaluate(Pagwas, bootstrap_iters = 200, part = 0.5)
Pagwas$Pathway_ld_gwas_data <- NULL
results <- data.frame(Pagwas$bootstrap_results)
write.csv(results, paste0("~/dat/complextraits/scpagwas/out/traits_celltype/", list2[[j]], "_CellLineage.csv"))
