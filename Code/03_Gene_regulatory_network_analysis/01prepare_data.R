#! /work/home/sdxgroup01/00envs/anaconda3/envs/R4.2.0/bin/R

library(optparse)
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(openxlsx)
library(tools)
library(grid)

op_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, action = "store", help = "The input of Seurat RDS", metavar = "rds"),
  # make_option(c("-a", "--anno"), type = "character", default = NULL, action = "store", help = "The cell type of annotation",metavar="annotation"),
  make_option(c("-o", "--outpath"), type = "character", default = "./", action = "store", help = "The outpath of output file", metavar = "opath")
)
parser <- OptionParser(option_list = op_list)
opt <- parse_args(parser)

tissue <- strsplit(basename(opt$input), "\\.")[[1]][1]
obj <- readRDS(opt$input)

homologenes <- read.table("~/01scrna/data/Gene_homologenes.txt", header = T)
homo <- distinct(homologenes, Gene_name, Human_gene_name)
tmp_matrix <- t(as.matrix(obj@assays$RNA@counts))
tmp_matrix1 <- tmp_matrix[,colnames(tmp_matrix) %in% homo$Gene_name]
tmp <- homo$Human_gene_name[match(colnames(tmp_matrix1), homo$Gene_name)]
colnames(tmp_matrix1) <- tmp

write.csv(tmp_matrix1, file = paste0(opt$outpath, "/", tissue, "_exp", ".csv"), quote = F)
write.csv(obj@meta.data, paste0(opt$outpath, "/", tissue, "_metadata.csv"), quote = F)
