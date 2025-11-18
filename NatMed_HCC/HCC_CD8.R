library(Seurat)
setwd("/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/")
laod("./rawdata/GSE206325_data_SingleCellExperiment_object.Rda")
metadata <- read.csv("./rawdata/GSE206325_full_HCC_cluster_annotation.csv")
seurat_obj <- as.Seurat(sce, counts = "counts", data = NULL)
metadata_req <- metadata[grep("CD4|Naive|Treg|Proliferating", metadata$type), ]

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

cellnames <- rownames(seurat_obj@meta.data[grep("^CD4$|^Naive$|^Treg$|^Proliferating$", seurat_obj@meta.data$celltype), ])
CD4 <- subset(seurat_obj, cells = cellnames)

CD4 <- NormalizeData(CD4, normalization.method = "LogNormalize", scale.factor = 10000)
CD4 <- FindVariableFeatures(CD4, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(CD4)
CD4 <- ScaleData(CD4, features = all.genes)
CD4 <- RunPCA(CD4, features = VariableFeatures(object = CD4))
CD4 <- FindNeighbors(CD4, dims = 1:30)
CD4 <- RunUMAP(CD4, dims = 1:30)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
dir.create(paste0(savedir, "umap"))
pdf(paste0(savedir, "umap/celltype.pdf"))
DimPlot(CD4, group.by = "celltype")
dev.off()

### Integration using Harmony
CD4 <- RunHarmony(CD4, group.by.vars = "cell_to_sample_ID")

CD4 <- RunUMAP(
    CD4,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)
savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "umap/celltype_integrated.pdf"))
DimPlot(CD4, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

CD4 <- FindNeighbors(CD4, reduction = "harmony", dims = 1:30)
CD4 <- FindClusters(CD4, resolution = 0.8)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "umap/celltype_integrated_clusters.pdf"))
DimPlot(CD4, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "umap/sample_integrated_clusters.pdf"))
DimPlot(CD4, reduction = "harmonyumap", group.by = "cell_to_cluster", label = TRUE, label.size = 4)
dev.off()

markers <- FindAllMarkers(CD4)

write.table(markers, paste0(savedir, "Table/markers.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

patterns <- metadata$cluster_ID
replacements <- metadata$subgroup

library(stringr)
CD4@meta.data$sub_celltype2 <- CD4@meta.data$cell_to_cluster
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    CD4@meta.data$sub_celltype2 <- str_replace_all(CD4@meta.data$sub_celltype2, pattern, replacements[i])
    CD4@meta.data$sub_celltype2 <- gsub("-", "", CD4@meta.data$sub_celltype2)
}

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "umap/sub_celltype_integrated_clusters.pdf"))
DimPlot(CD4, group.by = "sub_celltype2", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()
