library(Seurat)
PBC_Tcell <- readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/analysis/saveRDS/Tcell_subset_obj_4.RDS")

## 0,2,3,5,6,11,14
CD8_clus <- c(0, 2, 3, 5, 6, 11, 14)
CD8_cells <- rownames(PBC_Tcell@meta.data)[grep(paste0("^", CD8_clus, "$", collapse = "|"), PBC_Tcell@meta.data$seurat_clusters)]
CD8_Tcells <- subset(PBC_Tcell, cells = CD8_cells)
CD8_Tcells <- NormalizeData(CD8_Tcells)
CD8_Tcells <- FindVariableFeatures(CD8_Tcells)
CD8_Tcells <- ScaleData(CD8_Tcells)
CD8_Tcells <- RunPCA(CD8_Tcells)
CD8_Tcells <- RunHarmony(CD8_Tcells, group.by.vars = "orig.ident")

CD8_Tcells <- RunUMAP(
    CD8_Tcells,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/CD8/"
dir.create(paste0(savedir, "analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/CD8_integrated_umap.pdf"), width = 6, height = 5)
DimPlot(CD8_Tcells, reduction = "harmonyumap", label = TRUE, label.size = 4, group.by = "orig.ident")
dev.off()

CD8_Tcells <- FindNeighbors(CD8_Tcells, reduction = "harmony", dims = 1:30)
# CD8_Tcells <- FindClusters(CD8_Tcells, resolution = 0.8)
CD8_Tcells <- FindClusters(CD8_Tcells, resolution = 0.5)

dir.create(paste0(savedir, "analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/CD8_integrated_umap_res_0.5.pdf"), width = 6, height = 5)
DimPlot(CD8_Tcells, reduction = "harmonyumap", group.by = "seurat_clusters", label = TRUE, label.size = 4)
dev.off()


genes <- c("CD3E", "CD3G", "CD3D", "CD2", "CD4", "CD8A", "CD8B", "IL7R", "TCF7", "LEF1", "TBX21", "EOMES", "BCL6", "CD28", "FAS", "TNFSF8", "IL21", "IFNG", "GZMK", "TNF", "FOXP3", "MKI67", "CXCR6", "TYROBP", "LYZ")

dir.create(paste0(savedir, "analysis/featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/featureplots/Tcell_required.pdf"), width = 15, height = 15)
FeaturePlot(CD8_Tcells, genes, reduction = "harmonyumap", label = TRUE)
dev.off()

dir.create(paste0(savedir, "analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Tcell_nopoint_required.pdf"), width = 12, height = 12)
VlnPlot(CD8_Tcells, genes, pt.size = 0, group.by = "seurat_clusters")
dev.off()

Tcellmarkers <- FindAllMarkers(CD8_Tcells)
dir.create(paste0(savedir, "analysis/Table"), showWarnings = FALSE)
write.table(Tcellmarkers, paste0(savedir, "analysis/Table/PBC_CD8_Tcellmarkers.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

dir.create(paste0(savedir, "analysis/saveRDS"), showWarnings = FALSE)
saveRDS(CD8_Tcells, paste0(savedir, "analysis/saveRDS/PBC_CD8_Tcells.RDS"))
