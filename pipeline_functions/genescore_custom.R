# identifying the T cell susbet with Ucell
# signature gene from the paper suggested by Hirohisa Imbalance of Regulatory and Cytotoxic SARS-CoV-2-Reactive CD4+ T Cells in COVID-19
# https://www.sciencedirect.com/science/article/pii/S0092867420313076#bib23
# Single-cell transcriptomic analysis of allergen-specific T cells in allergy and asthma Th1 and Th2
# Supplementary Data at this location /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/Resource

genescore_custom_abhinav <- function(obj, assay="MAGIC_RNA", savedir, cores=12, objname, markername, markers_list, 
                                     reductions = "umap",x="UMAP_1",y="UMAP_2" , wide = 6, hite=5.5, dot_size=0.5){
  library(ArchR)
  library(ggplot2)
  source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/express_cell_front.R")
  DefaultAssay(obj) <- assay
  rm(markers)
  markers <- list()
  files <- markername
  Tcellsubset <- markers_list
  Tcellsubset <- rownames(obj)[match(Tcellsubset,rownames(obj), nomatch = 0)]
  message(paste(" ",Tcellsubset,sep=" "))
  markers[[files]] <- Tcellsubset

  obj <- AddModuleScore_UCell(obj, features = markers, ncores = cores)
  signature.names <- paste0(names(markers), "_UCell")
  p1 <- VlnPlot(obj, features = signature.names, group.by = "seurat_clusters")
  p3 <- VlnPlot(obj, features = signature.names, group.by = "seurat_clusters", pt.size = 0)
  p2 <- FeaturePlot(obj, features = signature.names, reduction = reductions)

  dir.create(paste(savedir,"VlnPlot",sep = ""),showWarnings = FALSE)
  pdf(paste(savedir,"VlnPlot/",objname,"_",assay,"_",markername,"_gene_score.pdf",sep = ""), width = 14, height = 12)
  print(p1)
  dev.off()

  dir.create(paste(savedir,"VlnPlot",sep = ""),showWarnings = FALSE)
  pdf(paste(savedir,"VlnPlot/",objname,"_",assay,"_",markername,"_gene_score_2.pdf",sep = ""), width = 14, height = 12)
  print(p3)
  dev.off()

  dir.create(paste(savedir,"featureplot",sep = ""),showWarnings = FALSE)
  rm(plot_list)
  plot_list <- list()
  p <- featureplot_front(obj, signature.names, reduction = reductions, x=x, y=y,size=dot_size) + scale_color_gradientn(colors = ArchRPalettes$solarExtra)

  pdf(paste(savedir,"featureplot/",objname,"_",assay,"_",markername,"_gene_score.pdf",sep = ""), width = wide, height = hite)
  print(p)
  dev.off()
  
  message("Plot saved at this location ",savedir)
  
  return(obj)
}


