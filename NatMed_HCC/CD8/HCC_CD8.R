#### Performing for CD8
library(Seurat)
setwd("/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/")
load("./rawdata/GSE206325_data_SingleCellExperiment_object.Rda")
metadata <- read.csv("./rawdata/GSE206325_full_HCC_cluster_annotation.csv")
seurat_obj <- as.Seurat(sce, counts = "counts", data = NULL)
metadata_req <- metadata[grep("CD8", metadata$type), ]

patterns <- metadata$cluster_ID
replacements <- metadata$type

library(stringr)
seurat_obj@meta.data$celltype <- seurat_obj@meta.data$cell_to_cluster
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    seurat_obj@meta.data$celltype <- str_replace_all(seurat_obj@meta.data$celltype, pattern, replacements[i])
}

patterns <- metadata$cluster_ID
replacements <- metadata$subgroup

library(stringr)
seurat_obj@meta.data$sub_celltype <- seurat_obj@meta.data$cell_to_cluster
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    seurat_obj@meta.data$sub_celltype <- str_replace_all(seurat_obj@meta.data$sub_celltype, pattern, replacements[i])
}

cellnames <- rownames(seurat_obj@meta.data[grep("^CD8$", seurat_obj@meta.data$celltype), ])
CD8 <- subset(seurat_obj, cells = cellnames)

CD8 <- NormalizeData(CD8, normalization.method = "LogNormalize", scale.factor = 10000)
CD8 <- FindVariableFeatures(CD8, selection.method = "vst", nfeatures = 3000)
CD8 <- ScaleData(CD8)
CD8 <- RunPCA(CD8)

### Integration using Harmony
library(harmony)
CD8 <- RunHarmony(CD8, group.by.vars = "cell_to_sample_ID")

CD8 <- RunUMAP(
    CD8,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/CD8/"
dir.create(paste0(savedir, "umap"), showWarning = FALSE)
pdf(paste0(savedir, "umap/celltype_integrated.pdf"))
DimPlot(CD8, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

CD8 <- FindNeighbors(CD8, reduction = "harmony", dims = 1:30)
CD8 <- FindClusters(CD8, resolution = 0.5)

pdf(paste0(savedir, "umap/celltype_integrated_clusters.pdf"))
DimPlot(CD8, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

pdf(paste0(savedir, "umap/sample_integrated_clusters.pdf"))
DimPlot(CD8, reduction = "harmonyumap", group.by = "cell_to_cluster", label = TRUE, label.size = 4)
dev.off()

markers <- FindAllMarkers(CD8)
write.table(markers, paste0(savedir, "Table/markers.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

patterns <- metadata$cluster_ID
replacements <- metadata$subgroup

library(stringr)
CD8@meta.data$sub_celltype2 <- CD8@meta.data$cell_to_cluster
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    CD8@meta.data$sub_celltype2 <- str_replace_all(CD8@meta.data$sub_celltype2, pattern, replacements[i])
    CD8@meta.data$sub_celltype2 <- gsub("-", "", CD8@meta.data$sub_celltype2)
}

pdf(paste0(savedir, "umap/sub_celltype_integrated_clusters.pdf"))
DimPlot(CD8, group.by = "sub_celltype2", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

genes <- c("CD3E", "CD3G", "CD3D", "CD2", "CD4", "CD8A", "CD8B", "IL7R", "TCF7", "LEF1", "TBX21", "EOMES", "BCL6", "CD28", "FAS", "TNFSF8", "IL21", "IFNG", "GZMK", "TNF", "FOXP3", "MKI67", "CXCR6", "TYROBP", "LYZ")

dir.create(paste0(savedir, "analysis/featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/featureplots/Tcell_required.pdf"), width = 15, height = 15)
FeaturePlot(CD8, genes, reduction = "harmonyumap", label = TRUE)
dev.off()

dir.create(paste0(savedir, "analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Tcell_nopoint_required.pdf"), width = 12, height = 12)
VlnPlot(CD8, genes, pt.size = 0, group.by = "seurat_clusters")
dev.off()

### removing these cluster
CD8_clus <- c(4, 10, 11, 5, 9, 12, 14, 15, 16)

CD8_cells <- rownames(CD8@meta.data)[grep(paste0("^", CD8_clus, "$", collapse = "|"), CD8@meta.data$seurat_clusters, invert = TRUE)]
CD8_2 <- subset(CD8, cells = CD8_cells)

CD8_2 <- NormalizeData(CD8_2, normalization.method = "LogNormalize", scale.factor = 10000)
CD8_2 <- FindVariableFeatures(CD8_2, selection.method = "vst", nfeatures = 3000)
CD8_2 <- ScaleData(CD8_2)
CD8_2 <- RunPCA(CD8_2)

### Integration using Harmony
library(harmony)
CD8_2 <- RunHarmony(CD8_2, group.by.vars = "cell_to_sample_ID")

CD8_2 <- RunUMAP(
    CD8_2,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/CD8_2/"
dir.create(paste0(savedir, "umap"), showWarning = FALSE, recursive = TRUE)
pdf(paste0(savedir, "umap/celltype_integrated.pdf"))
DimPlot(CD8_2, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

CD8_2 <- FindNeighbors(CD8_2, reduction = "harmony", dims = 1:30)
CD8_2 <- FindClusters(CD8_2, resolution = 0.5)

pdf(paste0(savedir, "umap/celltype_integrated_clusters.pdf"))
DimPlot(CD8_2, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

pdf(paste0(savedir, "umap/sample_integrated_clusters.pdf"))
DimPlot(CD8_2, reduction = "harmonyumap", group.by = "cell_to_cluster", label = TRUE, label.size = 4)
dev.off()

markers <- FindAllMarkers(CD8_2)
write.table(markers, paste0(savedir, "Table/markers.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

patterns <- metadata$cluster_ID
replacements <- metadata$subgroup

library(stringr)
CD8_2@meta.data$sub_celltype2 <- CD8_2@meta.data$cell_to_cluster
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    CD8_2@meta.data$sub_celltype2 <- str_replace_all(CD8_2@meta.data$sub_celltype2, pattern, replacements[i])
    CD8_2@meta.data$sub_celltype2 <- gsub("-", "", CD8_2@meta.data$sub_celltype2)
}

pdf(paste0(savedir, "umap/sub_celltype_integrated_clusters.pdf"))
DimPlot(CD8_2, group.by = "sub_celltype2", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

genes <- c("CD3E", "CD3G", "CD3D", "CD2", "CD4", "CD8A", "CD8B", "IL7R", "TCF7", "LEF1", "TBX21", "EOMES", "BCL6", "CD28", "FAS", "TNFSF8", "IL21", "IFNG", "GZMK", "TNF", "FOXP3", "MKI67", "CXCR6", "TYROBP", "LYZ")

dir.create(paste0(savedir, "analysis/featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/featureplots/Tcell_required.pdf"), width = 15, height = 15)
FeaturePlot(CD8_2, genes, reduction = "harmonyumap", label = TRUE)
dev.off()

dir.create(paste0(savedir, "analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Tcell_nopoint_required.pdf"), width = 12, height = 12)
VlnPlot(CD8_2, genes, pt.size = 0, group.by = "seurat_clusters")
dev.off()


# region removing more
CD8_clus <- c(9, 10, 12)

CD8_cells <- rownames(CD8_2@meta.data)[grep(paste0("^", CD8_clus, "$", collapse = "|"), CD8_2@meta.data$seurat_clusters, invert = TRUE)]
CD8_3 <- subset(CD8, cells = CD8_cells)

CD8_3 <- NormalizeData(CD8_3, normalization.method = "LogNormalize", scale.factor = 10000)
CD8_3 <- FindVariableFeatures(CD8_3, selection.method = "vst", nfeatures = 3000)
CD8_3 <- ScaleData(CD8_3)
CD8_3 <- RunPCA(CD8_3)

### Integration using Harmony
library(harmony)
CD8_3 <- RunHarmony(CD8_3, group.by.vars = "cell_to_sample_ID")

CD8_3 <- RunUMAP(
    CD8_3,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/CD8_3/"
dir.create(paste0(savedir, "umap"), showWarning = FALSE, recursive = TRUE)
pdf(paste0(savedir, "umap/celltype_integrated.pdf"))
DimPlot(CD8_3, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

CD8_3 <- FindNeighbors(CD8_3, reduction = "harmony", dims = 1:30)
CD8_3 <- FindClusters(CD8_3, resolution = 0.5)

pdf(paste0(savedir, "umap/celltype_integrated_clusters.pdf"))
DimPlot(CD8_3, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

pdf(paste0(savedir, "umap/sample_integrated_clusters.pdf"))
DimPlot(CD8_3, reduction = "harmonyumap", group.by = "cell_to_cluster", label = TRUE, label.size = 4)
dev.off()

markers <- FindAllMarkers(CD8_3)
dir.create(paste0(savedir, "Table"), showWarnings = FALSE)
write.table(markers, paste0(savedir, "Table/markers.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

patterns <- metadata$cluster_ID
replacements <- metadata$subgroup

library(stringr)
CD8_3@meta.data$sub_celltype2 <- CD8_3@meta.data$cell_to_cluster
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    CD8_3@meta.data$sub_celltype2 <- str_replace_all(CD8_3@meta.data$sub_celltype2, pattern, replacements[i])
    CD8_3@meta.data$sub_celltype2 <- gsub("-", "", CD8_3@meta.data$sub_celltype2)
}

pdf(paste0(savedir, "umap/sub_celltype_integrated_clusters.pdf"))
DimPlot(CD8_3, group.by = "sub_celltype2", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

genes <- c("CD3E", "CD3G", "CD3D", "CD2", "CD4", "CD8A", "CD8B", "IL7R", "TCF7", "LEF1", "TBX21", "EOMES", "BCL6", "CD28", "FAS", "TNFSF8", "IL21", "IFNG", "GZMK", "TNF", "FOXP3", "MKI67", "CXCR6", "TYROBP", "LYZ")

dir.create(paste0(savedir, "analysis/featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/featureplots/Tcell_required.pdf"), width = 15, height = 15)
FeaturePlot(CD8_3, genes, reduction = "harmonyumap", label = TRUE)
dev.off()

dir.create(paste0(savedir, "analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Tcell_nopoint_required.pdf"), width = 12, height = 12)
VlnPlot(CD8_3, genes, pt.size = 0, group.by = "seurat_clusters")
dev.off()

dir.create(paste0(savedir, "saveRDS"), showWarnings = FALSE)
saveRDS(CD8_3, paste0(savedir, "saveRDS/HCC_CD8.RDS"))
