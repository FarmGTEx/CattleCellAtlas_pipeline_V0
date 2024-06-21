
#### ================================================
#### AddModuleScore
#### ================================================

#### Load Packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggsci)
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(ggpubr)


#### 1. Load Data
setwd("~/dat/immune/06macro_merge")
sce <- readRDS("celltype.rds")
#### M1/M2
M1 <- c("SOCS1", "NOS2", "TNF", "CXCL9", "CXCL10", "CXCL11", "CD86", "IL1A", "IL1B", "IL6", "CCL5", "IRF5", "IRF1", "CCR7")
M2 <- c("IL4R", "CCL4", "CCL18", "CCL22", "MARCO", "VEGFA", "CTSA", "CTSB", "TGFB1", "MMP9", "CLEC7A", "MSR1", "IRF4", "CD163", "TGM2", "MRC1")
gene <- list(M1, M2)
#### APS
# MHC1 <- c('AZGP1','BOLA','JSP.1,''BOLA-NC1','CLEC4A','IKBKB','LNPEP','MFSD6','TAP1','ABCB9',
#     'ACE','CALR','ERAP1','ERAP2','HFE','IDE','MR1','PDIA3','SAR1B','TAPBP','TAPBPL')
# MHC2 <- c('CD74','CTSD','CTSF','CTSL','CTSS','CTSV','DNM2','FCGR2B','BOLA-DMA','BOLA-DMB','BOLA-DOA','BOLA-DOB','BOLA-DQA5',
#     'BOLA-DQB','BLA-DQB','BOLA-DRA','DSB','LGMN','PIKFYVE','TRAF6','ARL8B','MARCHF1','MARCHF8','PYCARD','THBS1','TREM2')


#### 2. Score the gene sets
df <- AddModuleScore(
    object = sce,
    features = gene,
    ctrl = 100, 
    name = c("M1", "M2")
)


#### 3. Save
dd <- FetchData(df, vars = c("tissue", "M11", "M22"))
d1 <- dd[, c(1, 2)]
d2 <- dd[, c(1, 3)]
d1$type <- "M11"
d2$type <- "M22"
names(d1) <- c("tissue", "MHC", "type")
names(d2) <- c("tissue", "MHC", "type")
df <- rbind(d1, d2)
write.csv(df, "M1_vs_M2_score.csv")


#### 4. Plot
p <- ggplot(df, aes(x = tissue, y = MHC)) +
    geom_boxplot(aes(fill = type), outlier.alpha = 0) +
    theme_bw() +
    stat_compare_means(
        aes(group = type),
        label = "p.signif",
        method = "t.test",
        hide.ns = T,
        size = 5,
        label.y = 1.7
    )
ggsave("M1_vs_M2_tissue.pdf", p, width = 15, height = 4)
