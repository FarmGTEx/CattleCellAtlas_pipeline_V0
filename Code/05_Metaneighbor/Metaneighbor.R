

#### =========================================================
#### metaNeighbor
#### =========================================================
#### Load Packages
library(MetaNeighbor)
library(SummarizedExperiment)
library(stringr)
library(ComplexHeatmap)

#### 1. Load Data
setwd("~/dat/Immune/03heatmap")
sce <- readRDS("celltype.rds")

#### 2. Create Dataset
count_matrix <- as.matrix(sce@assays$RNA@counts)
meta_data <- sce@meta.data
mn_data <- SummarizedExperiment(
        assays = list(exprs = count_matrix),
        colData = meta_data,
        rowData = NULL
)

#### 3. Search Variable Genes
var_genes <- variableGenes(dat = mn_data, exp_labels = mn_data$tissue)

#### 4. MetaNeighbor analysis
result <- MetaNeighborUS(
        var_genes = var_genes,
        dat = mn_data,
        study_id = as.factor(mn_data$cell_type),
        cell_type = as.factor(mn_data$tissue),
        fast_version = TRUE
)

#### 5. Save
write.csv(result, "Metaneighbor_tissue.csv", quote = F)

row.names(result) <- sapply(strsplit(row.names(result), "\\|"), function(x) x[2])
colnames(result) <- sapply(strsplit(colnames(result), "\\|"), function(x) x[2])

#### 6. Plot
cairo_pdf("Metaneighbor_tissue.pdf")
ComplexHeatmap::pheatmap(result,
        color = colorRampPalette(c("#ffffffc0", "#fefebd", "#fdb674", "#d73027"))(100),
        show_colnames = F,
        show_rownames = F,
        heatmap_legend_param = list(title = "AUROC"),
        fontsize = 10
)
dev.off()
