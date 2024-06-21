#! /work/home/sdxgroup01/00envs/anaconda3/envs/sc_cellchat/bin/Rscript
library(optparse)
library(Seurat)
library(CellChat)
library(cowplot)
options(stringsAsFactors = FALSE) 

op_list <- list(
make_option(c("-i", "--input"), type = "character", default = NULL, action = "store", help = "The input of Seurat RDS",metavar="rds"),
make_option(c("-o", "--outpath"), type = "character", default = "./", action = "store", help = "The outpath of output file",metavar="opath")
)
parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)

CellCommunication.CellChat <- function(infile, outdir, species = NULL){
  outpdf <- file.path(outdir, paste0(tissue,"communication.cellchat_global.pdf"))
  pdf(outpdf,width = 8, height =10)
  obj <- readRDS(infile)
  data.input <- Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
  ## add C in seurat cluster
  meta <- data.frame(labels = obj$"cell_type", row.names = names(obj$label)) 
  cellchat <- CellChat::createCellChat(object = data.input, meta = meta, group.by = "labels")

  # CellChat database
  ## current human and mouse are available
  if (species == "human"){
    CellChatDB <- CellChat::CellChatDB.human  
  }else if (species == "mouse"){
    CellChatDB <- CellChat::CellChatDB.mouse
  } else {
    print("Only support human or mouse!!!")
  }
  plot_tmp <- CellChat::showDatabaseCategory(CellChatDB)
  print(plot_tmp)
  # Show the structure of the database
  dplyr::glimpse(CellChatDB$interaction)
  # use a subset of CellChatDB for cell-cell communication analysis
  #CellChatDB.use <- CellChat::subsetDB(CellChatDB) 
  cellchat@DB <- CellChatDB

  # preprocessing the expression data
  cellchat <- CellChat::subsetData(cellchat)
  #future::plan("multiprocess", workers = 4) # do parallel
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  
  # Part II: Inference of cell-cell communication network
  ## Compute the communication probability and infer cellular communication network
  cellchat <- CellChat::computeCommunProb(cellchat)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  df.net <- CellChat::subsetCommunication(cellchat)
  write.table(df.net, file.path(outdir, paste0(tissue,"_CellChat.inferred_network_ligands.tsv")))
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  df.netp <- CellChat::subsetCommunication(cellchat, slot.name = "netP")
  write.table(df.netp, file.path(outdir, paste0(tissue,"_CellChat.inferred_network_pathway.tsv")))
  ## Calculate the aggregated cell-cell communication network
  cellchat <- CellChat::aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  p1 <- CellChat::netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  p2 <- CellChat::netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
 ##
  outpdf <- paste0("single_cell_chat_number_weight.pdf")
  pdf(outpdf, width = 12, height = 15)
  mat1 <- cellchat@net$count
  mat2 <- cellchat@net$weight
  numRows <- ceiling(length(groupSize) / 1)
  par(mfrow = c(1, 2), xpd = TRUE)
  for (i in 1:nrow(mat1)) {
    mat3 <- matrix(0, nrow = nrow(mat1), ncol = ncol(mat1), dimnames = dimnames(mat1))
    mat3[i, ] <- mat1[i, ]
    netVisual_circle(mat3, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2, vertex.label.cex = 0.7, arrow.size = 0.1, edge.weight.max = max(mat1), title.name = paste0("Number of interactions in ",rownames(mat1)[i]))

    mat4 <- matrix(0, nrow = nrow(mat2), ncol = ncol(mat2), dimnames = dimnames(mat2))
    mat4[i, ] <- mat2[i, ]
    netVisual_circle(mat4, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2, vertex.label.cex = 0.7, arrow.size = 0.1, edge.weight.max = max(mat2), title.name = paste0("Interaction weights/strength ",rownames(mat2))[i])
    }
  dev.off()

  # Part III: Visualization of cell-cell communication network
  # plot a specified pathway
  pathways.show <- head(cellchat@netP$pathways,n=20)
  outpdf <- file.path(outdir, paste0(tissue,"communication.cellchat_output_contributation.pdf"))
  pdf(outpdf,width = 8, height =10)
  # vertex.receiver = as.numeric(which(rownames(mat1) %in% c("Luminal cells1", "Luminal cells2", "Luminal cells3", "Luminal cells4", "Myoepithelial")))
  vertex.receiver = c(1,2,3)  #define a numeric vector giving the index of the celltype as targes
  for (pathway in pathways.show){
    par(mfrow = c(1,1), xpd=TRUE)
    p1 <- CellChat::netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy") 
    #p2 <- CellChat::netVisual_aggregate(cellchat, signaling = pathway, layout = "chord")
    #p2 <- CellChat::netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
    #plot base on single object, could not use par
    p3 <- CellChat::netVisual_heatmap(cellchat, signaling = pathway, color.heatmap = "Reds")
    print(p3)
    ### Compute the contribution  
    p_contribute <- CellChat::netAnalysis_contribution(cellchat, signaling = pathway)
    print(p_contribute)

    pairLR <- CellChat::extractEnrichedLR(cellchat, signaling = pathway, geneLR.return = FALSE)
    #### show one ligand-receptor pair
    LR.show <- pairLR[1,] 
    #vertex.receiver = seq(1,4) #define a numeric vector giving the index of the celltype as targes 
    #p4 <- CellChat::netVisual_individual(cellchat, signaling = pathway, pairLR.use = LR.show, layout = "circle", vertex.receiver = vertex.receiver)
    p4 <- CellChat::netVisual_individual(cellchat, signaling = pathway, pairLR.use = LR.show, vertex.receiver = vertex.receiver)
    print(p4)
    #p5 <- CellChat::netVisual_individual(cellchat, signaling = pathway, pairLR.use = LR.show, layout = "chord")
    #print(p5)
    p6 <-CellChat::plotGeneExpression(cellchat, signaling = pathway)
    print(p6)
    p7 <- CellChat::netVisual_bubble(cellchat, sources.use = vertex.receiver, targets.use = c(1:nrow(mat1)), remove.isolate = FALSE) 
    print(p7)
  }
  dev.off()
  saveRDS(cellchat, file = file.path(outdir, paste0(tissue,"communication.cellchat_final.rds")))
  p8 <- CellChat::netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:3), signaling = cellchat@netP$pathways, remove.isolate = FALSE)
  print(p8)
  pairLR.use <- CellChat::extractEnrichedLR(cellchat, signaling = cellchat@netP$pathways)
  p9 <- CellChat::netVisual_bubble(cellchat, sources.use = c(1,2), targets.use = c(2:3), pairLR.use = pairLR.use, remove.isolate = TRUE)
  print(p9)
  p10 <- CellChat::netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(3:5), lab.cex = 0.5,legend.pos.y = 30)
  print(p10)
  p11 <- CellChat::netVisual_chord_gene(cellchat, sources.use = c(1,2), targets.use = c(3:4), slot.name = "netP", legend.pos.x = 10)
  print(p11)
  
  
  # Part IV: Systems analysis of cell-cell communication network
  outpdf <- file.path(outdir, paste0(tissue,"communication.cellchat_output_network.pdf"))
  pdf(outpdf,width = 12, height =12)
  cellchat <- CellChat::netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  CellChat::netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 6, font.size = 10)
  ## Visualize the dominant senders (sources) and receivers (targets) in a 2D space
  gg1 <- CellChat::netAnalysis_signalingRole_scatter(cellchat)
  gg2 <- CellChat::netAnalysis_signalingRole_scatter(cellchat,signaling = cellchat@netP$pathways)
  # gg2 <- CellChat::netAnalysis_signalingRole_scatter(cellchat,signaling = "SPP1")
  print(gg1 + gg2)
  ## Identify signals contributing most to outgoing or incoming signaling
  ht1 <- CellChat::netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width = 10, height =25)
  ht2 <- CellChat::netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width = 10, height =25)
  print(ht1+ht2)
  dev.off()

  ## Identify and visualize outgoing communication pattern of secreting cells
  library(NMF)
  library(ggalluvial)
  outpdf <- file.path(outdir, paste0(tissue,"communication.cellchat_output_selectK.pdf"))
  pdf(outpdf,width = 8, height =10)
  p12 <- CellChat::selectK(cellchat, pattern = "outgoing") #Please run the function `selectK` for selecting a suitable k
  print(p12)
  dev.off()
  #### use the count of drop suddenly in last step
  outpdf <- file.path(outdir, paste0(tissue,"communication.cellchat_output_Patterns.pdf"))
  pdf(outpdf,width = 8, height =10)
  nPatterns = 3
  cellchat <- CellChat::identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  p13 <- CellChat::netAnalysis_river(cellchat, pattern = "outgoing")
  print(p13)
  p14 <- CellChat::netAnalysis_dot(cellchat, pattern = "outgoing")
  print(p14)
  ### Identify and visualize incoming communication pattern of target cells
  p15 <- CellChat::selectK(cellchat, pattern = "incoming")
  print(p15)
  cellchat <- CellChat::identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
  p16 <- CellChat::netAnalysis_river(cellchat, pattern = "incoming")
  print(p16)
  p17 <-CellChat::netAnalysis_dot(cellchat, pattern = "incoming")
  print(p17)
  dev.off()
  ## Manifold and classification learning analysis of signaling network
  ### Identify signaling groups based on their functional similarity
  outpdf <- file.path(outdir, paste0(tissue,"communication.cellchat_output_netsimilarity.pdf"))
  pdf(outpdf,width = 8, height =10)
  par(mfrow = c(1, 2), xpd = TRUE)
  cellchat <- CellChat::computeNetSimilarity(cellchat, type = "functional")
  cellchat <- CellChat::netEmbedding(cellchat, type = "functional")
  cellchat <- CellChat::netClustering(cellchat, type = "functional")
  p18 <- CellChat::netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
  print(p18)
  ### Identify signaling groups based on structure similarity
  cellchat <- CellChat::computeNetSimilarity(cellchat, type = "structural")
  cellchat <- CellChat::netEmbedding(cellchat, type = "structural")
  cellchat <- CellChat::netClustering(cellchat, type = "structural")
  p19 <- CellChat::netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
  print(p19)
  dev.off()
  # Part V: Save the CellChat object
  
  
}
inrds <- opt$input
outdir <- opt$outpath
species <- "human"
tissue <- strsplit(basename(inrds),"\\.")[[1]][1]
result <- CellCommunication.CellChat(inrds, outdir, species)