
#### ====================================================================
#### hdWGCNA
#### ====================================================================
#### Load Packages
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(qlcMatrix)
library(harmony)
library(enrichR)
library(GeneOverlap)

theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 10)


#### ==================================================================
#### 1. Load Data
#### ==================================================================
setwd("/work/home/sdxgroup01/01user/zwj/scRNA/Epithelial/07remun/Spinous")
data <- readRDS("celltype.rds")
Idents(data) <- "seurat_clusters"
cell <- "spinous"


#### ==================================================================
#### 2. Obtain information of target cell_type
#### ==================================================================
#### 2.1 Setting Seurat Objects for WGCNA
seurat_obj <- SetupForWGCNA(
    data,
    gene_select = "fraction", # or variable, default 2000
    fraction = 0.05, # Select genes expressed in at least 5% of cells in the dataset
    wgcna_name = paste0(cell, "_hdWGCNA") # the name of the hdWGCNA experiment
)
table(seurat_obj@meta.data$tissue)
table(seurat_obj@meta.data$seurat_clusters)

#### 2.2 Building Metacells
seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("cell_type", "seurat_clusters"), # specify the columns in seurat_obj@meta.data to group by
    reduction = "harmony", # select the dimensionality reduction to perform KNN on
    k = 25, # nearest-neighbors parameter
    max_shared = 10, # maximum number of shared cells between two metacells
    ident.group = "seurat_clusters", # set the Idents of the metacell seurat object
    min_cells = 50 # minimum cell count threshold in group
)


#### ==================================================================
#### 3. Co-expression Network Analysis
#### ==================================================================
#### 3.1 Establish Expression Matrix
seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = unique(seurat_obj@meta.data$seurat_clusters), # the name of the group of interest in the group.by column
    group.by = "seurat_clusters", # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = "RNA", # using RNA assay
    slot = "data" # using normalized data
    # use_metacells = TRUE # use the metacells (TRUE) or the full expression matrix (FALSE)
)

#### 3.2 Test Different Soft Powers
seurat_obj <- TestSoftPowers(
    seurat_obj,
    networkType = "signed" # can also use "unsigned" or "signed hybrid"
)
plot_list <- PlotSoftPowers(seurat_obj, point_size = 5, text_size = 3)
wrap_plots(plot_list, ncol = 2)
ggsave("./01_softPowers.pdf", dpi = 300)

#### 3.3 Get Threshold
power_table <- GetPowerTable(seurat_obj)
a <- power_table$SFT.R.sq
i <- 1
for (b in a) {
    if (b < 0.8) {
        i <- i + 1
        print(i)
    } else if (b > 0.8) {
        break
    }
}
select_soft_power <- power_table$Power[i]

#### 3.4 Construct Co-expression Network
seurat_obj <- ConstructNetwork(
    seurat_obj,
    soft_power = select_soft_power,
    setDatExpr = FALSE,
    overwrite_tom = TRUE,
    tom_name = cell # name of the topoligical overlap matrix written to disk
)
cairo_pdf("./02_Dendrogram.pdf")
PlotDendrogram(seurat_obj, main = paste0(cell, " hdWGCNA Dendrogram"))
dev.off()


#### ==================================================================
#### 4. Co expression network module analysis (module feature genes and connectivity)
#### ==================================================================
#### 4.1 Compute Coordination Module Feature Genes
#### need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
#### compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars = "tissue"
)
#### Harmonized Module Eigengenes
hMEs <- GetMEs(seurat_obj)
#### Module eEigengenes
MEs <- GetMEs(seurat_obj, harmonized = FALSE)

#### 4.2 Calculate correlation between modules and phenotypes
cur_traits <- c("nCount_RNA", "nFeature_RNA")
seurat_obj <- ModuleTraitCorrelation(
    seurat_obj,
    traits = cur_traits,
    group.by = "tissue"
)
mt_cor <- GetModuleTraitCorrelation(seurat_obj)

cairo_pdf("./03_module_trait_correlation.pdf")
PlotModuleTraitCorrelation(
    seurat_obj,
    label = "fdr",
    label_symbol = "stars",
    text_size = 2,
    text_digits = 2,
    text_color = "black",
    high_color = "red",
    mid_color = "white",
    low_color = "blue",
    plot_max = 0.2,
    combine = T
)
dev.off()


#### 4.3 Identify Hub Genes
#### 4.3.1 Compute Eigengene-based Connectivity (kME)
seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = "cell_type", group_name = cell
)

#### 4.3.2 Rename hdWGCNA Module
seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name = paste0(cell, "-M")
)

#### 4.3.3 Get the Module Assignment Table:
modules <- GetModules(seurat_obj)
write.csv(modules, "./module_gene.csv", quote = FALSE)
#### Plot
cairo_pdf("./04_PlotKMEs.pdf", width = 25, height = 15)
PlotKMEs(seurat_obj, n_hubs = 10, ncol = 4, text_size = 6, plot_widths = c(3, 2))
dev.off()
#### Get Hub Genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
write.csv(hub_df, "./hub_gene.csv", quote = FALSE)

#### Save
saveRDS(seurat_obj, file = paste0(cell, "_hdWGCNA.rds"))

#### 4.3.4 Calculate Hub Gene Score
seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = 25,
    method = "Seurat"
)

#### 4.4 Make a Featureplot of hMEs for Each Module
plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features = "hMEs", # plot the hMEs
    order = TRUE # order so the points with highest hMEs are on top
)
wrap_plots(plot_list)
ggsave("./05_module_features.pdf", height = 10, width = 15, dpi = 300)

#### Make a Featureplot of hMEs Score for Each Module
plot_list <- ModuleFeaturePlot(
    seurat_obj,
    reduction = "umap",
    features = "scores", # plot the hub gene scores
    order = "shuffle", # order so cells are shuffled
    ucell = TRUE # depending on Seurat vs UCell for gene scoring
)
wrap_plots(plot_list)
ggsave("./06_module_scores.pdf", height = 10, width = 15, dpi = 300)

#### 4.5 Module Association Diagram
cairo_pdf("./07_module_relation.pdf")
ModuleCorrelogram(seurat_obj, exclude_grey = TRUE, features = "hMEs")
dev.off()

#### 4.6 Module Feature Genes
#### get hMEs from Seurat Object
MEs <- GetMEs(seurat_obj, harmonized = TRUE)
MEs <- MEs[, order(names(MEs))]
mods <- colnames(MEs)
mods <- mods[mods != "grey"]
#### Add hMEs to Seurat Metadata
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
#### Plot
p <- DotPlot(seurat_obj, features = mods, group.by = "tissue")+
    coord_flip() +
    RotatedAxis() +
    scale_color_gradient2(high = "red", mid = "grey95", low = "blue")
ggsave("./08_gene_expression.pdf", dpi = 300)

#### Plot hMEs using Seurat VlnPlot function
plot_list <- lapply(mods, function(x) {
    print(x)
    p <- VlnPlot(
        seurat_obj,
        features = x,
        group.by = "tissue",
        pt.size = 0 # don't show actual data points
    )
    p <- p + geom_boxplot(width = .25, fill = "white")
    p <- p + xlab("") + ylab("hME") + NoLegend() +
        theme(
            axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            plot.title = element_text(size = 20)
        )
    p
})
wrap_plots(plot_list)
ggsave("./09_gene_expression_vinplot.pdf", dpi = 300, width = 20, height = 15)

#### 4.7 Differential Module Feature Genes (DME) Analysis
DMEs_all <- FindAllDMEs(seurat_obj, group.by = "tissue")
p <- PlotDMEsVolcano(
    seurat_obj,
    DMEs_all,
    plot_labels = FALSE,
    show_cutoff = FALSE
)
p + facet_wrap(~group)
ggsave("./08_DMEs.pdf", height = 20, width = 10, dpi = 300)

#### 4.8 Module Network Diagram
#### Single
ModuleNetworkPlot(
    seurat_obj,
    mods = "all",
    outdir = "ModuleNetworks",
    plot_size = c(6, 6),
    label_center = FALSE,
    edge.alpha = 0.25,
    vertex.label.cex = 1,
    vertex.size = 6
)

#### Combine
cairo_pdf("./09_module_net.pdf")
HubGeneNetworkPlot(
    seurat_obj,
    n_hubs = 3, n_other = 5,
    edge_prop = 0.75,
    mods = "all",
    vertex.label.cex = 0.5,
    hub.vertex.size = 6,
    other.vertex.size = 2,
)
dev.off()

#### 4.9 UMAP Network Diagram
seurat_obj <- RunModuleUMAP(
    seurat_obj,
    n_hubs = 10, # number of hub genes to include for the UMAP embedding
    n_neighbors = 15, # neighbors parameter for UMAP
    min_dist = 0.1 # min distance between points in UMAP space
)
#### Get the Hub Genes UMAP 
umap_df <- GetModuleUMAP(seurat_obj)
p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(
        color = umap_df$color, # color each point by WGCNA module
        size = umap_df$kME * 2 # size of each point based on intramodular connectivity
    ) +
    umap_theme()
ggsave("./10umap_module.pdf", p, dpi = 300)

cairo_pdf("./11_umap_net.pdf")
ModuleUMAPPlot(
    seurat_obj,
    edge.alpha = 0.25,
    sample_edges = TRUE,
    edge_prop = 0.1, # proportion of edges to sample (20% here)
    label_hubs = 1, # how many hub genes to plot per module?
    keep_grey_edges = FALSE
)
dev.off()


#### ==================================================================
#### 6. Module Cell (tissue) Marker Gene Overlap Analysis
#### ==================================================================
#### 6.1 Calculate Tissue Marker Genes
Idents(seurat_obj) <- seurat_obj$tissue
markers <- Seurat::FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    logfc.threshold = 1
)

#### 6.2 Calculate Overlap Genes Between Tissues and Modules
overlap_df <- OverlapModulesDEGs(
    seurat_obj,
    deg_df = markers,
    fc_cutoff = 1 # log fold change cutoff for overlap analysis
)
plot_list <- OverlapBarPlot(overlap_df)
wrap_plots(plot_list)
ggsave("./12_overlap_gene.pdf", width = 15, height = 30, dpi = 300)

#### Plot Odds Ratio of the Overlap 
p <- OverlapDotPlot(overlap_df, plot_var = "odds_ratio") +
    ggtitle(paste0("Overlap of modules & ", cell, " markers"))
ggsave("./13_overlap_gene_all.pdf", p, width = 15, height = 30, dpi = 300)
