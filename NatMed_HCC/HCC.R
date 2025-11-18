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

# region subsetting
cellnames <- rownames(CD4@meta.data[grep("^8$|^17$|^18$|^19$|^20$|^10$|^14$|^16$|^15$", CD4@meta.data$seurat_clusters, invert = TRUE), ])
CD4_2 <- subset(CD4, cells = cellnames)

CD4_2 <- NormalizeData(CD4_2, normalization.method = "LogNormalize", scale.factor = 10000)
CD4_2 <- FindVariableFeatures(CD4_2, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(CD4_2)
CD4_2 <- ScaleData(CD4_2, features = all.genes)
CD4_2 <- RunPCA(CD4_2, features = VariableFeatures(object = CD4_2))
CD4_2 <- FindNeighbors(CD4_2, dims = 1:30)
CD4_2 <- RunUMAP(CD4_2, dims = 1:30)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
dir.create(paste0(savedir, "umap"))
pdf(paste0(savedir, "umap/celltype_refined.pdf"))
DimPlot(CD4_2, group.by = "celltype")
dev.off()

### Integration using Harmony
CD4_2 <- RunHarmony(CD4_2, group.by.vars = "cell_to_sample_ID")

CD4_2 <- RunUMAP(
    CD4_2,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)
savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "umap/celltype_integrated_refined.pdf"))
DimPlot(CD4_2, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "umap/celltype_integrated_clusters_previous.pdf"))
DimPlot(CD4_2, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

CD4_2@meta.data$seurat_clusters2 <- (CD4_2@meta.data$seurat_clusters)

CD4_2 <- FindNeighbors(CD4_2, reduction = "harmony", dims = 1:30)
CD4_2 <- FindClusters(CD4_2, resolution = 0.4)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "umap/celltype_integrated_clusters_refined_0.4.pdf"))
DimPlot(CD4_2, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "umap/sample_integrated_clusters_refined.pdf"))
DimPlot(CD4_2,
    reduction = "harmonyumap", group.by = "cell_to_cluster",
    label = TRUE, label.size = 4
)
dev.off()

markers <- FindAllMarkers(CD4_2)

write.table(markers, paste0(savedir, "Table/markers_refined_0.4.txt"), quote = F, row.names = T, col.names = T, sep = "\t")
saveRDS(CD4_2, "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/saveRDS/CD4_2.RDS")

# endregion

# region subsetting 9 and 10
cellnames <- rownames(CD4_2@meta.data[grep("^9$|^10$", CD4_2@meta.data$seurat_clusters, invert = TRUE), ])
CD4_3 <- subset(CD4_2, cells = cellnames)

CD4_3 <- NormalizeData(CD4_3, normalization.method = "LogNormalize", scale.factor = 10000)
CD4_3 <- FindVariableFeatures(CD4_3, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(CD4_3)
CD4_3 <- ScaleData(CD4_3, features = all.genes)
CD4_3 <- RunPCA(CD4_3, features = VariableFeatures(object = CD4_3))
CD4_3 <- FindNeighbors(CD4_3, dims = 1:30)
CD4_3 <- RunUMAP(CD4_3, dims = 1:30)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
dir.create(paste0(savedir, "umap"))
pdf(paste0(savedir, "umap/celltype_refined_3.pdf"))
DimPlot(CD4_3, group.by = "celltype")
dev.off()

### Integration using Harmony
CD4_3 <- RunHarmony(CD4_3, group.by.vars = "cell_to_sample_ID")

CD4_3 <- RunUMAP(
    CD4_3,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "umap/celltype_integrated_refined_3.pdf"))
DimPlot(CD4_3, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

pdf(paste0(savedir, "umap/subtype_integrated_refined_3.pdf"))
DimPlot(CD4_3, group.by = "sub_celltype", reduction = "harmonyumap", label = TRUE, label.size = 6) + NoLegend()
dev.off()

# savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
# pdf(paste0(savedir, "umap/celltype_integrated_clusters_previous.pdf"))
# DimPlot(CD4_3, reduction = "harmonyumap", label = TRUE, label.size = 4)
# dev.off()

CD4_3@meta.data$seurat_clusters3 <- CD4_3@meta.data$seurat_clusters
CD4_3 <- FindNeighbors(CD4_3, reduction = "harmony", dims = 1:30)


genes <- c(
    "TCF7", "LEF1", "IL7R", "GPR183", "LTB", "CD28", "ICOS", "MAF", "TNFSF8", "CXCL13", "IL21",
    "CSF2", "IFNG", "TNF", "FOXP3", "IL2RA", "ENTPD1", "CTLA4", "PDCD1", "TIGIT", "HAVCR2", "GZMK",
    "CXCR5", "CXCR6", "CXCR3", "CCR7"
)

savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
res = c("0.8")
for i in 1:length(res){
    CD4_3 <- FindClusters(CD4_3, resolution = as.integer(res[i]), cluster.name = paste0("res",res[i]))
    
    savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
    pdf(paste0(savedir, "umap/clusters_refined_",res[i],"_integrated.pdf"))
    DimPlot(CD4_3, reduction = "harmonyumap", group.by = paste0("res",res[i]), label = TRUE, label.size = 6)
    dev.off()

    markers <- FindAllMarkers(CD4_3, group.by = paste0("res",res[i]))
    write.table(markers, paste0(savedir, "Table/markers_refined_",res[i],".txt"), 
    quote = F, row.names = T, col.names = T, sep = "\t")

    dir.create(paste0(savedir, "vlnplots"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "vlnplots/Tcell_refined_3_",res[i],".pdf"), width = 15, height = 15)
    VlnPlot(CD4_3, genes, pt.size = 0, group.by = paste0("res",res[i]))
    dev.off()
}

dir.create(paste0(savedir, "featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/Tcell_refined_3.pdf"), width = 15, height = 15)
FeaturePlot(CD4_3, genes, reduction = "harmonyumap", label = TRUE)
dev.off()

saveRDS(CD4_3, "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/saveRDS/CD4_3.RDS")

# endregion

# region subsetting 9 and 10
cellnames <- rownames(CD4_3@meta.data[grep("^17$", CD4_3@meta.data$res0.8, invert = TRUE), ])
CD4_4 <- subset(CD4_3, cells = cellnames)

CD4_4 <- NormalizeData(CD4_4, normalization.method = "LogNormalize", scale.factor = 10000)
CD4_4 <- FindVariableFeatures(CD4_4, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(CD4_4)
CD4_4 <- ScaleData(CD4_4, features = all.genes)
CD4_4 <- RunPCA(CD4_4, features = VariableFeatures(object = CD4_4))
CD4_4 <- FindNeighbors(CD4_4, dims = 1:30)
CD4_4 <- RunUMAP(CD4_4, dims = 1:30)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
dir.create(paste0(savedir, "umap"))
pdf(paste0(savedir, "umap/celltype_refined_3.pdf"))
DimPlot(CD4_4, group.by = "celltype")
dev.off()

### Integration using Harmony
CD4_4 <- RunHarmony(CD4_4, group.by.vars = "cell_to_sample_ID")

CD4_4 <- RunUMAP(
    CD4_4,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "umap/celltype_integrated_refined_4.pdf"))
DimPlot(CD4_4, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

pdf(paste0(savedir, "umap/subtype_integrated_refined_4.pdf"))
DimPlot(CD4_4, group.by = "sub_celltype", reduction = "harmonyumap", label = TRUE, label.size = 6) + NoLegend()
dev.off()

# savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
# pdf(paste0(savedir, "umap/celltype_integrated_clusters_previous.pdf"))
# DimPlot(CD4_4, reduction = "harmonyumap", label = TRUE, label.size = 4)
# dev.off()

# CD4_4@meta.data$seurat_clusters3 <- CD4_4@meta.data$seurat_clusters
CD4_4 <- FindNeighbors(CD4_4, reduction = "harmony", dims = 1:30)

genes <- c(
    "TCF7", "LEF1", "IL7R", "GPR183", "LTB", "CD28", "ICOS", "MAF", "TNFSF8", "CXCL13", "IL21",
    "CSF2", "IFNG", "TNF", "FOXP3", "IL2RA", "ENTPD1", "CTLA4", "PDCD1", "TIGIT", "HAVCR2", "GZMK",
    "CXCR5", "CXCR6", "CXCR3", "CCR7"
)

savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
res = c("0.4","0.8")
i=2
for i in 1:length(res){
    CD4_4 <- FindClusters(CD4_4, resolution = 0.8, cluster.name = paste0("res",res[i]))
    
    savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
    pdf(paste0(savedir, "umap/clusters_refined_",res[i],"_integrated_2.pdf"))
    print(DimPlot(CD4_4, reduction = "harmonyumap", group.by = paste0("res",res[i]), label = TRUE, label.size = 6))
    dev.off()

    markers <- FindAllMarkers(CD4_4, group.by = paste0("res",res[i]))
    write.table(markers, paste0(savedir, "Table/markers_refined_",res[i],"_2.txt"), 
    quote = F, row.names = T, col.names = T, sep = "\t")

    dir.create(paste0(savedir, "vlnplots"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "vlnplots/Tcell_refined_3_",res[i],"_2.pdf"), width = 15, height = 15)
    print(VlnPlot(CD4_4, genes, pt.size = 0, group.by = paste0("res",res[i])))
    dev.off()
}

dir.create(paste0(savedir, "featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/Tcell_refined_4_res0.8.pdf"), width = 15, height = 15)
FeaturePlot(CD4_4, genes, reduction = "harmonyumap", label = TRUE)
dev.off()

saveRDS(CD4_4, "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/saveRDS/CD4_4.RDS")

cellnames <- rownames(CD4_4@meta.data[grep("^17$", CD4_4@meta.data$res0.8, invert = TRUE), ])
CD4_5 <- subset(CD4_4, cells = cellnames)

dir.create(paste0(savedir, "featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/Tcell_refined_5_res0.8.pdf"), width = 15, height = 20)
FeaturePlot(CD4_5, genes, reduction = "harmonyumap", label = FALSE)
dev.off()

res = c("0.4","0.8")
dir.create(paste0(savedir, "vlnplots"), showWarnings = FALSE, recursive = TRUE)
for(i in 1:length(res)){
    pdf(paste0(savedir, "vlnplots/Tcell_refined_5_",res[i],"_2.pdf"), width = 15, height = 15)
    print(VlnPlot(CD4_5, genes, pt.size = 0, group.by = paste0("res",res[i])))
    dev.off()
}
    
savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "umap/clusters_refined_0.8_integrated_noC17.pdf"))
DimPlot(CD4_5, reduction = "harmonyumap", group.by = "res0.8", label = TRUE, label.size = 6)
dev.off()

pdf(paste0(savedir, "umap/clusters_refined_0.4_integrated_noC17.pdf"))
DimPlot(CD4_5, reduction = "harmonyumap", group.by = "res0.4", label = TRUE, label.size = 6)
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "umap/celltype_integrated_refined_noC17.pdf"))
DimPlot(CD4_5, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 6)
dev.off()

pdf(paste0(savedir, "umap/subtype_integrated_refined_noC17.pdf"))
DimPlot(CD4_5, group.by = "sub_celltype", reduction = "harmonyumap", label = TRUE, label.size = 6) + NoLegend()
dev.off()

cluster = c(0:8)
celltypes=c(
    "Post-stem CD4 T cells 1","Treg 1","Stem-like CD4 T cell","Tph",
    "Dual producer","Effector/memory","Treg 2","GZMK Th1","Post-stem CD4 T cells 2"
    )

patterns <- cluster
replacements <- celltypes

library(stringr)
CD4_5@meta.data$celltypes <- CD4_5@meta.data$res0.4
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    CD4_5@meta.data$celltypes <- str_replace_all(CD4_5@meta.data$celltypes, pattern, replacements[i])
}

savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "umap/celltype_CD4_HCC_Yuki_annot.pdf"))
DimPlot(CD4_5, group.by = "celltypes", reduction = "harmonyumap", label = TRUE, label.size = 6) + NoLegend()
dev.off()

saveRDS(CD4_5, "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/saveRDS/CD4_5.RDS")

#region CD4 Trajectory
library(Seurat)
CD4_obj=readRDS("/mnt/data/projects/NatMed_HCC/analysis/saveRDS/CD4_5.RDS")

#### Violin plots
Set1 <- c("TCF7", "LEF1", "SELL", "CCR7", "FOXP3", "TIGIT")
Set2 <- c("CTLA4", "IFNG", "TNF", "CCL5", "CD40LG", "GZMK")
Set3 <- c("HSPA1A", "HSPA1B", "HSPD1", "NR4A1", "HSPA6", "BAG3")
Set4 <- c("CXCL13", "PDCD1", "TNFSF8", "IL21", "TOX2", "MAF", "GPR183")
Set5 <- c("ZEB2", "ID2", "XCL1", "XCL2", "CCL4", "GZMA")

sets <- c("Set1","Set2","Set3","Set4","Set5")
savedir = "/mnt/data/projects/NatMed_HCC/analysis/"

for(i in 1:length(sets)){
    # pdf(paste0(savedir,"vlnplots/",sets[i],"_PBC_HCC_genes.pdf"))
    # print(VlnPlot(PBC_HCC,get(sets[i])))
    # dev.off()
    pdf(paste0(savedir,"vlnplots/",sets[i],"_HCC_genes_no_points.pdf"))
    print(VlnPlot(CD4_obj,get(sets[i]), pt.size = 0))
    dev.off()
}


### SlingShot
# Performing the trajectory it is better to take into consideration the global structure rather than the local structure. 
# So we use reduced dimension UMAP rather than tSNE
# SlingShot: multiple branching lineages.
# Slingshot breaks the inference problem into two steps, we are able to make use of appropriate methods for each task and 
# avoid the common
# trade-off between stability and the flexibility to detect complex structures.
# Using a cluster-based MST for lineage inference allows Slingshot to identify potentially complex global patterns in the data without
# being overly sensitive to individual data points.
# And our novel simultaneous principal curves method for pseudotime inference extends the stability and robustness properties of principal
# curves to the case of multiple branching lineages.
#
# 2) the inference of pseudotime variables for cells along each lineage (simultaneous principal curves) to smooth the curves.
# multiple lineage inference into two stages:
# 1. Identification of lineages, i.e., ordered sets of cell clusters, where all lineages share a starting cluster and each leads to a
# unique terminal cluster. This is achieved by constructing an MST on clusters of cells.
# 2. For each lineage, identification of pseudotimes, i.e., a one-dimensional variable representing each cell’s transcriptional progression
# toward the terminal state.
#
# Robustness to noise.
# The Monocle procedure, which constructs an MST on individual cells and orders them according to a PQ tree along the longest path of the
# MST, was the least stable of the methods we compared. The path drawn by Monocle was highly variable and sensitive to even small amounts
# of noise; this instability has been.
#
# the vertices of the piecewise linear path drawn by the cluster-based MST, multiple cells will often be assigned identical pseudotimes,
# corresponding to the value at the vertex. principal curve approach was the most stable method, but on more complex datasets, it has the
# obvious limitation of only characterizing a single lineage. It is for this reason that we chose to extend principal curves to accommodate
# multiple branching lineages.
#
# Slingshot allows the use of a shape-sensitive distance measure inspired by the Mahalanobis distance which scales the distance between
# cluster centers based on the covariance structure of the two clusters.

# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html
suppressPackageStartupMessages({
  library(slingshot)
  library(SingleCellExperiment)
  library(RColorBrewer)
  library(scales)
  library(viridis)
  library(UpSetR)
  library(pheatmap)
  library(msigdbr)
  library(fgsea)
  library(knitr)
  library(ggplot2)
  library(gridExtra)
  library(tradeSeq)
  library(Seurat)
  library(bioc2020trajectories)
})

# Save the objects as separate matrices for input in slingshot
## convert back to singleCellExperiment
DefaultAssay(CD4_obj) <- "originalexp"
sce <- as.SingleCellExperiment(CD4_obj, assay = "originalexp")

# The question is: should we fit a separate trajectory for each condition? We might expect the trajectory itself to be changed by the
# treatment if the treatment effect is systematically large. Otherwise, the treatment may impact the expression profile of some genes but
# the overall trajectory will be preserved.

# We can calculate the imbalance score Regions with a high score indicate that the local cell distribution according to treatment label is unbalanced
# compared the overall distribution.Here, we see that, while there are some small regions of imbalance, the global path along the development axis is well-balanced.
# This means that we can fit a global trajectory to the full dataset, so this approach ensures that our trajectory accounts for all cell types present in the overall data.

# Since in this there is no condition so we will not use this funciton
# scores <- bioc2020trajectories::imbalance_score(
#   rd = reducedDims(sce)$UMAP,
#   cl = colData(sce)$pheno$treatment_id,
#   k = 20, smooth = 40)

# The goal of slingshot is to use clusters of cells to uncover global structure and convert this structure into smooth lineages represented
# by one-dimensional variables, called “pseudotime.” We provide tools for learning cluster relationships in an unsupervised or semi-supervised
# manner and constructing smooth curves representing each lineage,

# The minimal input to slingshot is a matrix representing the cells in a reduced-dimensional space and a vector of cluster labels.
# With these two inputs, we then:
# 1. Identify the global lineage structure by constructing an minimum spanning tree (MST) on the clusters, with the getLineages function.
# 2. Construct smooth lineages and infer pseudotime variables by fitting simultaneous principal curves with the getCurves function.
# 3. Assess the output of each step with built-in visualization tools.

# However, we recommend clustering the cells even in datasets where only a single lineage is expected, as it allows for the potential discovery of novel branching events.
# The clusters identified in this step will be used to determine the global structure of the underlying lineages (that is, their number,
# when they branch off from one another, and the approximate locations of those branching events). This is different than the typical goal
# of clustering single-cell data, which is to identify all biologically relevant cell types present in the dataset. For example, when
# determining global lineage structure, there is no need to distinguish between immature and mature neurons since both cell types will,
# presumably, fall along the same segment of a lineage.

# These two steps can be run separately with the getLineages and getCurves functions, or together with the wrapper function, slingshot (recommended).

# sce_auto <- slingshot(sce,
#                  clusterLabels = 'seurat_clusters',
#                  reducedDim = 'UMAP', approx_points = 100)
#
# library(grDevices)
# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
#
# pdf(paste(savedir,"Trajectory/slingshot_trajectory.pdf",sep = ""))
# plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
# lines(SlingshotDataSet(sce), lwd=2, col='black')
# dev.off()
# the outgroup is an artificial cluster that is equidistant from all real clusters at some threshold value. If the original MST sans the
# outgroup contains an edge that is longer than twice the threshold, the addition of the outgroup will cause the MST to instead be routed
# through the outgroup.

## Starting Cluster
savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/trajectory/"
dir.create(savedir, showWarnings = FALSE)
dir.create(paste(savedir, "UMAP/", sep = ""), showWarnings = FALSE)
p <- DimPlot(CD4_obj, reduction = "harmonyumap", label = TRUE, group.by = "res0.4")
pdf(paste(savedir, "UMAP/CD4_res0.4.pdf", sep = ""))
p
dev.off()

# CD4_sce <- slingshot(sce,
#   clusterLabels = "seurat_clusters2",
#   reducedDim = "UMAP", approx_points = 100,
#   omega = TRUE,
#   omega_scale = 1.5
# )

CD4_sce_naive <- slingshot(sce,
  clusterLabels = colData(sce)$res0.4,
  reducedDim = "HARMONYUMAP",
  start.clus = c("2"),
  # end.clus = c("7"),
  approx_points = 150
  # omega = TRUE,omega_scale=1.5
)

library(ArchR)
library(scater)
embedded_orig <- embedCurves(CD4_sce, "UMAP")
rm(plot_list)
plot_list <- list()

for (i in 1:8) {
  embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
  embedded <- data.frame(embedded$s[embedded$ord, ])
  g <- plotUMAP(CD4_sce, colour_by = paste("slingPseudotime_", i, sep = ""), dimred = "HARMONYUMAP",)
  stopifnot(all(rownames(combine@reductions$umap@cell.embeddings) == rownames(g$data)))
  data <- merge(combine@reductions$umap@cell.embeddings, g$data, by = "row.names")
  colnames(data) <- c("cellname", "UMAP_1", "UMAP_2", "X", "Y", paste("Lineage", i, sep = ""))
  p <- ggplot(data, aes_string("UMAP_1", "UMAP_2", color = paste("Lineage", i, sep = ""))) +
    geom_point(size = 0.01) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
    geom_path(data = embedded, aes(x = UMAP_1, y = UMAP_2), color = "black", size = 1.2) +
    theme_bw()
  plot_list[[i]] <- p
}

require(gridExtra)
pdf(paste(savedir, "UMAP/sce_UMAP_splitted_unbias.pdf", sep = ""), width = 12, height = 9)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]], nrow = 3, ncol = 3)
dev.off()

library(ArchR)
library(scater)
embedded_orig <- embedCurves(CD4_sce_naive, "HARMONYUMAP")

rm(plot_list)
plot_list <- list()
for (i in 1:2) {
  embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
  embedded <- data.frame(embedded$s[embedded$ord, ])
  g <- plotUMAP(CD4_sce_naive, colour_by = paste("slingPseudotime_", i, sep = ""), dimred = "HARMONYUMAP")
  stopifnot(all(rownames(CD4_obj@reductions$umap@cell.embeddings) == rownames(g$data)))
  # data <- merge(embedded_orig@reductions$umap@cell.embeddings, g$data, by = "row.names")
  # colnames(data) <- c("cellname", "HARMONYUMAP_1", "HARMONYUMAP_2", "X", "Y", paste("Lineage", i, sep = ""))
  colnames(g$data) <- c("HARMONYUMAP_1","HARMONYUMAP_2",paste("Lineage", i, sep = ""))
  p <- ggplot(g$data, aes_string("HARMONYUMAP_1", "HARMONYUMAP_2", color = paste("Lineage", i, sep = ""))) +
    geom_point(size = 0.01) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
    geom_path(data = embedded, aes(x = harmonyUMAP_1, y = harmonyUMAP_2), color = "black", size = 1.2) +
    theme_bw()
  plot_list[[i]] <- p
}

require(gridExtra)
pdf(paste(savedir, "UMAP/HCC_CD4_clus_2_lineages.pdf", sep = ""), width = 8, height = 3.5)
grid.arrange(plot_list[[1]], plot_list[[2]], nrow = 1, ncol = 2)
dev.off()

pdf(paste(savedir, "UMAP/HCC_CD4_UMAP_lineage_2.pdf", sep = ""), width = 4.5, height = 4)
plot_list
dev.off()

#region cluster2 from 0.4
cellnames <- rownames(CD4_5@meta.data[grep("^2$", CD4_5@meta.data$res0.4, invert = FALSE), ])
CD4_clus2 <- subset(CD4_5, cells = cellnames)

CD4_clus2 <- NormalizeData(CD4_clus2, normalization.method = "LogNormalize", scale.factor = 10000)
CD4_clus2 <- FindVariableFeatures(CD4_clus2, selection.method = "vst", nfeatures = 3000)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/clus2/"
dir.create(paste0(savedir, "umap"))
pdf(paste0(savedir, "umap/celltype_refined_clus2.pdf"))
DimPlot(CD4_clus2, group.by = "celltype")
dev.off()

### Integration using Harmony
CD4_clus2 <- RunHarmony(CD4_clus2, group.by.vars = "cell_to_sample_ID")

CD4_clus2 <- RunUMAP(
    CD4_clus2,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

pdf(paste0(savedir, "umap/celltype_integrated_clus2.pdf"))
DimPlot(CD4_clus2, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 6)
dev.off()

pdf(paste0(savedir, "umap/subtypes_integrated_clus2.pdf"))
DimPlot(CD4_clus2, group.by = "sub_celltype", reduction = "harmonyumap", label = TRUE, label.size = 5) + NoLegend()
dev.off()

# savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
# pdf(paste0(savedir, "umap/celltype_integrated_clusters_previous.pdf"))
# DimPlot(CD4_clus2, reduction = "harmonyumap", label = TRUE, label.size = 4)
# dev.off()

# CD4_clus2@meta.data$seurat_clusters3 <- CD4_clus2@meta.data$seurat_clusters
CD4_clus2 <- FindNeighbors(CD4_clus2, reduction = "harmony", dims = 1:30)

genes <- c(
    "TCF7", "LEF1", "IL7R", "GPR183", "LTB", "CD28", "ICOS", "MAF", "TNFSF8", "CXCL13", "IL21",
    "CSF2", "IFNG", "TNF", "FOXP3", "IL2RA", "ENTPD1", "CTLA4", "PDCD1", "TIGIT", "HAVCR2", "GZMK",
    "CXCR5", "CXCR6", "CXCR3", "CCR7"
)

res = c("0.4","0.8")
# i=2
dir.create(paste0(savedir,"Table"),showWarnings = FALSE)
dir.create(paste0(savedir,"vlnplots"),showWarnings = FALSE)

pdf(paste0(savedir, "umap/clusters_refined_",res[i],"_integrated_2.pdf"))
print(DimPlot(CD4_clus2, reduction = "harmonyumap", group.by = paste0("res",res[i]), label = TRUE, label.size = 6))
dev.off()

for (i in 1:length(res)){
    CD4_clus2 <- FindClusters(CD4_clus2, resolution = as.numeric(res[i]), cluster.name = paste0("res",res[i]))
    
    pdf(paste0(savedir, "umap/clusters_refined_",res[i],"_integrated_2.pdf"))
    print(DimPlot(CD4_clus2, reduction = "harmonyumap", group.by = paste0("res",res[i]), label = TRUE, label.size = 6))
    dev.off()

    markers <- FindAllMarkers(CD4_clus2, group.by = paste0("res",res[i]))
    write.table(markers, paste0(savedir, "Table/markers_refined_",res[i],"_2.txt"), 
    quote = F, row.names = T, col.names = T, sep = "\t")

    dir.create(paste0(savedir, "vlnplots"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "vlnplots/Tcell_refined_3_",res[i],"_2.pdf"), width = 15, height = 15)
    print(VlnPlot(CD4_clus2, genes, pt.size = 0, group.by = paste0("res",res[i])))
    dev.off()
}

savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/clus2/"
for (i in 1:length(res)){
    dir.create(paste0(savedir, "vlnplots"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "vlnplots/Tcell_refined_3_",res[i],"_2.pdf"), width = 15, height = 15)
    print(VlnPlot(CD4_clus2, genes, pt.size = 0, group.by = paste0("res",res[i])))
    dev.off()
}

pdf(paste0(savedir, "umap/sample_clus2_",res[i],"_integrated_2.pdf"))
print(DimPlot(CD4_clus2, reduction = "harmonyumap", group.by = "cell_to_sample_ID", label = FALSE, label.size = 6)) + NoLegend()
dev.off()

dir.create(paste0(savedir, "featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/Tcell_refined_4_res0.8.pdf"), width = 15, height = 20)
FeaturePlot(CD4_clus2, genes, reduction = "harmonyumap", label = TRUE)
dev.off()

saveRDS(CD4_clus2, "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/saveRDS/CD4_clus2.RDS")

#endregion

#region cluster3 from 0.4
cellnames <- rownames(CD4_5@meta.data[grep("^3$", CD4_5@meta.data$res0.4, invert = FALSE), ])
CD4_clus3 <- subset(CD4_5, cells = cellnames)

CD4_clus3 <- NormalizeData(CD4_clus3, normalization.method = "LogNormalize", scale.factor = 10000)
CD4_clus3 <- FindVariableFeatures(CD4_clus3, selection.method = "vst", nfeatures = 3000)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/clus3/"
dir.create(paste0(savedir, "umap"))
pdf(paste0(savedir, "umap/celltype_refined_clus3.pdf"))
DimPlot(CD4_clus3, group.by = "celltype")
dev.off()

### Integration using Harmony
CD4_clus3 <- RunHarmony(CD4_clus3, group.by.vars = "cell_to_sample_ID")

CD4_clus3 <- RunUMAP(
    CD4_clus3,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

pdf(paste0(savedir, "umap/celltype_integrated_clus3.pdf"))
DimPlot(CD4_clus3, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 6)
dev.off()

pdf(paste0(savedir, "umap/subtypes_integrated_clus3.pdf"))
DimPlot(CD4_clus3, group.by = "sub_celltype", reduction = "harmonyumap", label = TRUE, label.size = 5) + NoLegend()
dev.off()

# savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
# pdf(paste0(savedir, "umap/celltype_integrated_clusters_previous.pdf"))
# DimPlot(CD4_clus3, reduction = "harmonyumap", label = TRUE, label.size = 4)
# dev.off()

# CD4_clus3@meta.data$seurat_clusters3 <- CD4_clus3@meta.data$seurat_clusters
CD4_clus3 <- FindNeighbors(CD4_clus3, reduction = "harmony", dims = 1:30)

genes <- c(
    "TCF7", "LEF1", "IL7R", "GPR183", "LTB", "CD28", "ICOS", "MAF", "TNFSF8", "CXCL13", "IL21",
    "CSF2", "IFNG", "TNF", "FOXP3", "IL2RA", "ENTPD1", "CTLA4", "PDCD1", "TIGIT", "HAVCR2", "GZMK",
    "CXCR5", "CXCR6", "CXCR3", "CCR7"
)

res = c("0.4","0.8")
# i=2
dir.create(paste0(savedir,"Table"),showWarnings = FALSE)
dir.create(paste0(savedir,"vlnplots"),showWarnings = FALSE)

pdf(paste0(savedir, "umap/clusters_refined_",res[i],"_integrated_2.pdf"))
print(DimPlot(CD4_clus3, reduction = "harmonyumap", group.by = paste0("res",res[i]), label = TRUE, label.size = 6))
dev.off()

for (i in 1:length(res)){
    CD4_clus3 <- FindClusters(CD4_clus3, resolution = as.numeric(res[i]), cluster.name = paste0("res",res[i]))
    
    pdf(paste0(savedir, "umap/clusters_refined_",res[i],"_integrated_2.pdf"))
    print(DimPlot(CD4_clus3, reduction = "harmonyumap", group.by = paste0("res",res[i]), label = TRUE, label.size = 6))
    dev.off()

    markers <- FindAllMarkers(CD4_clus3, group.by = paste0("res",res[i]))
    write.table(markers, paste0(savedir, "Table/markers_refined_",res[i],"_2.txt"), 
    quote = F, row.names = T, col.names = T, sep = "\t")

    dir.create(paste0(savedir, "vlnplots"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "vlnplots/Tcell_refined_3_",res[i],"_2.pdf"), width = 15, height = 15)
    print(VlnPlot(CD4_clus3, genes, pt.size = 0, group.by = paste0("res",res[i])))
    dev.off()
}

for (i in 1:length(res)){
    dir.create(paste0(savedir, "vlnplots"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(savedir, "vlnplots/Tcell_refined_3_",res[i],"_2.pdf"), width = 15, height = 15)
    print(VlnPlot(CD4_clus3, genes, pt.size = 0, group.by = paste0("res",res[i])))
    dev.off()
}

pdf(paste0(savedir, "umap/sample_clus3_",res[i],"_integrated_2.pdf"))
print(DimPlot(CD4_clus3, reduction = "harmonyumap", group.by = "cell_to_sample_ID", label = FALSE, label.size = 6)) + NoLegend()
dev.off()

dir.create(paste0(savedir, "featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/Tcell_refined_4_res0.8.pdf"), width = 15, height = 20)
FeaturePlot(CD4_clus3, genes, reduction = "harmonyumap", label = TRUE)
dev.off()

saveRDS(CD4_clus3, "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/saveRDS/CD4_clus3.RDS")
#endregion

#region Gene module score
# Loop through each column
library(Seurat)
HCC = readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/saveRDS/CD4_5.RDS")

DefaultAssay(HCC) <- "RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/.immune/Kyoto/resource/CD4_genes/", pattern = "txt", full.names = TRUE)
filename <- gsub(".txt", "", basename(files))
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[, 1]
    Tcellsubset <- rownames(HCC)[match(Tcellsubset, rownames(HCC), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

HCC <- AddModuleScore(HCC, features = markers, slot = "data")
colnames(HCC@meta.data)[17:33] <- paste0("CD4_", filename)

Idents(HCC) <- HCC@meta.data$res0.4

savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir, "vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/CD4_pancancer_HCC_module_score_genes.pdf"), width = 12, height = 12)
VlnPlot(HCC, paste0("CD4_", filename), pt.size = 0)
dev.off()

dir.create(paste0(savedir, "featureplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/CD4_pancancer_module_score_genes.pdf"), width = 12, height = 12)
FeaturePlot(HCC, paste0("CD4_", filename), reduction = "umap", label = TRUE)
dev.off()

### CD8
DefaultAssay(HCC) <- "RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/.immune/Kyoto/resource/CD8_genes/", pattern = "txt", full.names = TRUE)
filename <- gsub(".txt", "", basename(files))
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[, 1]
    Tcellsubset <- rownames(HCC)[match(Tcellsubset, rownames(HCC), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

HCC <- AddModuleScore(HCC, features = markers, slot = "data")
colnames(HCC@meta.data)[34:52] <- paste0("CD8_", filename)

# savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/HCC/"
dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir, "vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/CD8_pancancer_module_score_genes_violin.pdf"), width = 12, height = 12)
VlnPlot(HCC, paste0("CD8_", filename), pt.size = 0)
dev.off()

dir.create(paste0(savedir, "featureplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/CD8_pancancer_module_score_genes_featureplot.pdf"), width = 12, height = 12)
FeaturePlot(HCC, paste0("CD8_", filename), reduction = "umap", label = TRUE)
dev.off()

DefaultAssay(HCC) <- "RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/.immune/Kyoto/liver/resources/genelist/", pattern = "CD", full.names = TRUE)
filename <- paste0(basename(files), "_2")
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[, 1]
    Tcellsubset <- rownames(HCC)[match(Tcellsubset, rownames(HCC), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

HCC <- AddModuleScore(HCC, features = markers, slot = "data")

colnames(HCC@meta.data)[82:92] <- filename

savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"

dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir, "vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/Cell_2017_violinplot_HCC.pdf"), width = 12, height = 9)
VlnPlot(HCC, filename, pt.size = 0)
dev.off()


#endregion

#region Tregs
library(Seurat)
HCC = readRDS("/mnt/data/projects/NatMed_HCC/analysis/saveRDS/CD4_5.RDS")

savedir = "/mnt/data/projects/NatMed_HCC/analysis/Tregs/"
dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE, recursive = T)
pdf(paste0(savedir,"UMAP/CD4_Tcell.pdf"))
DimPlot(HCC, label = T, label.size = 8, reduction = "harmonyumap")
dev.off()

# cellnames <- rownames(seurat_obj@meta.data[grep("^CD4$|^Naive$|^Treg$|^Proliferating$", seurat_obj@meta.data$celltype), ])
CD4 <- subset(HCC, idents = c(1,6))

pdf(paste0(savedir,"UMAP/Treg_previous_UMAP.pdf"))
DimPlot(CD4, label = T, label.size = 8, reduction = "harmonyumap")
dev.off()

CD4 <- NormalizeData(CD4, normalization.method = "LogNormalize", scale.factor = 10000)
CD4 <- FindVariableFeatures(CD4, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(CD4)
CD4 <- ScaleData(CD4, features = all.genes)
CD4 <- RunPCA(CD4, features = VariableFeatures(object = CD4))
CD4 <- FindNeighbors(CD4, dims = 1:30)
CD4 <- RunUMAP(CD4, dims = 1:30)

pdf(paste0(savedir, "UMAP/Treg_unintegrated.pdf"))
DimPlot(CD4)
dev.off()

### Integration using Harmony
library(harmony)
CD4 <- RunHarmony(CD4, group.by.vars = "cell_to_sample_ID")

CD4 <- RunUMAP(
    CD4,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

pdf(paste0(savedir, "UMAP/Treg_integrated.pdf"))
DimPlot(CD4, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

pdf(paste0(savedir, "UMAP/Treg_integrated_samples.pdf"))
DimPlot(CD4, group.by = "cell_to_sample_ID", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

CD4 <- FindNeighbors(CD4, reduction = "harmony", dims = 1:30)
CD4 <- FindClusters(CD4, resolution = 0.6)

pdf(paste0(savedir, "UMAP/celltype_integrated_clusters.pdf"))
DimPlot(CD4, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/"
pdf(paste0(savedir, "UMAP/sample_integrated_clusters.pdf"))
DimPlot(CD4, reduction = "harmonyumap", group.by = "cell_to_cluster", label = TRUE, label.size = 4)
dev.off()

markers <- FindAllMarkers(CD4)

dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(markers, paste0(savedir, "Table/Tregs_markers.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

### removing the cluster 12 and 14
Tregs_cells=rownames(CD4@meta.data[grep("^12$|^14$",CD4@meta.data$seurat_clusters,invert=T),])

savedir = "/mnt/data/projects/NatMed_HCC/analysis/Tregs/removed_clus12_14/"
dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE, recursive = T)
pdf(paste0(savedir,"UMAP/CD4_Tcell.pdf"))
DimPlot(CD4, label = T, label.size = 8, reduction = "harmonyumap")
dev.off()

# cellnames <- rownames(seurat_obj@meta.data[grep("^CD4$|^Naive$|^Treg$|^Proliferating$", seurat_obj@meta.data$celltype), ])
Tregs <- subset(CD4, cells = Tregs_cells)

pdf(paste0(savedir,"UMAP/Treg_previous_UMAP.pdf"))
DimPlot(Tregs, label = T, label.size = 8, reduction = "harmonyumap")
dev.off()

Tregs <- NormalizeData(Tregs, normalization.method = "LogNormalize", scale.factor = 10000)
Tregs <- FindVariableFeatures(Tregs, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Tregs)
Tregs <- ScaleData(Tregs, features = all.genes)
Tregs <- RunPCA(Tregs, features = VariableFeatures(object = Tregs))
Tregs <- FindNeighbors(Tregs, dims = 1:30)
Tregs <- RunUMAP(Tregs, dims = 1:30)

pdf(paste0(savedir, "UMAP/Treg_unintegrated.pdf"))
DimPlot(Tregs)
dev.off()

### Integration using Harmony
library(harmony)
Tregs <- RunHarmony(Tregs, group.by.vars = "cell_to_sample_ID")

Tregs <- RunUMAP(
    Tregs,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

pdf(paste0(savedir, "UMAP/Treg_integrated.pdf"))
DimPlot(Tregs, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

pdf(paste0(savedir, "UMAP/Treg_integrated_samples.pdf"))
DimPlot(Tregs, group.by = "cell_to_sample_ID", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

Tregs <- FindNeighbors(Tregs, reduction = "harmony", dims = 1:30)
Tregs <- FindClusters(Tregs, resolution = 0.6, )

pdf(paste0(savedir, "UMAP/celltype_integrated_clusters.pdf"))
DimPlot(Tregs, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

markers <- FindAllMarkers(Tregs)

dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(markers, paste0(savedir, "Table/Tregs_markers_Res0.15.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

genes = c("FOXP3", "IKZF2", "TIGIT", "CTLA4", "ICOS", "HLA-DRB1", "TNFRSF9",
"TCF7", "LEF1", "CXCR6", "CXCR3", "IL10", "TBX21", "RORC",
"BCL6", "PDCD1", "CXCR5", "IL21R", "CCR7", "IL2RA", "FCRL3" )

dir.create(paste0(savedir,"vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/gene_vlnplots_res0.15.pdf"), width = 12, height = 15)
VlnPlot(Tregs, genes, pt.size = 0)
dev.off()

genes = c("FOXP3", "IKZF2", "TIGIT", "CTLA4", "ICOS", "HLA-DRB1", "TNFRSF9",
"TCF7", "LEF1", "CXCR6", "CXCR3", "IL10", "TBX21", "RORC",
"BCL6", "PDCD1", "CXCR5", "IL21R", "CCR7", "IL2RA", "FCRL3")

dir.create(paste0(savedir,"vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/gene_vlnplots_res0.6.pdf"), width = 12, height = 15)
VlnPlot(Tregs, genes, pt.size = 0, group.by = "originalexp_snn_res.0.6")
dev.off()

metadata = read.csv("/mnt/data/projects/NatMed_HCC/analysis/GSE206325_sample_annots_Liver_Treated_patients.csv", header = T)
patterns <- metadata$sample_ID
replacements <- metadata$treatment_Resp

library(stringr)
Tregs@meta.data$treatment_Resp <- Tregs@meta.data$cell_to_sample_ID
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    Tregs@meta.data$treatment_Resp <- str_replace_all(Tregs@meta.data$treatment_Resp, pattern, replacements[i])
    Tregs@meta.data$treatment_Resp <- gsub("-", "", Tregs@meta.data$treatment_Resp)
}

Tregs_response = rownames(Tregs@meta.data[grep("antiPD1",Tregs@meta.data$treatment_Resp),])
Tregs_response_obj = subset(Tregs, cells = Tregs_response)

pdf(paste0(savedir, "UMAP/treatment_response_splitted.pdf"))
DimPlot(Tregs_response_obj, group.by = "treatment_Resp", split.by = "treatment_Resp", reduction = "harmonyumap", label = FALSE, label.size = 4)
dev.off()

#### Longitudinal analysis
celltypes_ind <- table(Tregs_response_obj@meta.data$cell_to_sample_ID,Tregs_response_obj@meta.data$seurat_clusters)
dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(celltypes_ind, paste0(savedir,"Table/sample_clus_res_0.15.txt"),
            row.names = T, col.names = T, sep = "\t", quote = F)

# celltypes_ind <- read.table(paste("/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/Table/celltype_individual_res_0.8.txt"), header = TRUE, sep = "\t")
n_cells <- celltypes_ind
n_cells_sum <- as.vector(rowSums(celltypes_ind))

### Making an empty dataframe
df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
rownames(df) <- rownames(n_cells)
colnames(df) <- colnames(n_cells)

for (j in 1:nrow(n_cells)) {
  df[j,] <- (n_cells[j,]/n_cells_sum[j])*100
}

library(reshape2)
df$sample <- rownames(df)
df_melted <- melt(df)
colnames(df_melted) <- c("sample","celltype","percentage")
library(ggplot2)

df_melted$sample <- gsub("-","_",df_melted$sample)

library(ggpubr)
library(rstatix)
library(stringr)

metadata = read.csv("/mnt/data/projects/NatMed_HCC/analysis/GSE206325_sample_annots_Liver_Treated_patients.csv")
patterns <- metadata$sample_ID
replacements <- metadata$treatment_Resp

library(stringr)
df_melted$condition <- df_melted$sample
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    df_melted$condition <- str_replace_all(df_melted$condition, pattern, replacements[i])
    df_melted$condition <- gsub("-", "", df_melted$condition)
}

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=sample)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table/sample_cluster_res_0.15.pdf"),width =20, height = 7)
p
dev.off()

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=condition)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table//sample_cluster_res_0.15_2.pdf"),width =10, height = 7)
p
dev.off()

stat.test <- df_melted %>%
  group_by(celltype) %>%
  t_test(percentage ~ condition, paired = FALSE, alternative = "two.sided") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
# stat.test$p <- (stat.test$p)/2

bxp <- ggboxplot(
  df_melted, x = "celltype", y = "percentage", 
  color = "condition", palette = c("#00AFBB", "#E7B800")
) 

bxp2 <- bxp +   geom_dotplot(
  aes(fill = condition, color = condition), trim = FALSE,
  binaxis='y', stackdir='center', dotsize = 0.15,
  position = position_dodge(0.8)
)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))


# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "celltype", dodge = 0.8)
bxp3 <- bxp2 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) +   theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
            panel.background = element_rect(fill = 'white', colour = 'white'),
            panel.grid.minor = element_line(colour = "white"),
            panel.grid.major = element_line(colour = "white"))

# Add 10% spaces between the p-value labels and the plot border
bxp4 <- bxp3 + stat_pvalue_manual(
  stat.test,  label = "p.adj", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

pdf(paste(savedir,"Table/t_test_unpaired_all_sample_0.15.pdf",sep = ""), width = 8, height = 6)
bxp4
dev.off()

# Convert any list columns to simple character columns
stat_test_df_clean <- data.frame(lapply(stat_test_df, function(x) {
  if (is.list(x)) sapply(x, paste, collapse = ",") else x
}), stringsAsFactors = FALSE)

# Then write to file
write.table(
  stat_test_df_clean,
  paste0(savedir, "Table/stat_test.txt"),
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)

dir.create(paste0(savedir,"saveRDS"))
saveRDS(Tregs, paste0(savedir,"saveRDS/Treg_removed_7_8.RDS"))

#region Responders
library(Seurat)
savedir = "/mnt/data/projects/NatMed_HCC/analysis/Tregs/removed_clus12_14/"
setwd(savedir)
Tregs = readRDS(paste0(savedir,"saveRDS/Treg_removed_7_8.RDS"))

Tregs_responder = Tregs[,grep("antiPD1_R",Tregs@meta.data$treatment_Resp)]

Tregs_responder <- NormalizeData(Tregs_responder, normalization.method = "LogNormalize", scale.factor = 10000)
Tregs_responder <- FindVariableFeatures(Tregs_responder, selection.method = "vst", nfeatures = 3000)
Tregs_responder <- ScaleData(Tregs_responder)
Tregs_responder <- RunPCA(Tregs_responder, features = VariableFeatures(object = Tregs_responder))
Tregs_responder <- FindNeighbors(Tregs_responder, dims = 1:30)
Tregs_responder <- RunUMAP(Tregs_responder, dims = 1:30)

savedir = "/mnt/data/projects/NatMed_HCC/analysis/Tregs/removed_clus12_14/responder/"
setwd(savedir)
dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "UMAP/Treg_responder_unintegrated.pdf"))
DimPlot(Tregs_responder)
dev.off()

### Integration using Harmony
library(harmony)
Tregs_responder <- RunHarmony(Tregs_responder, group.by.vars = "cell_to_sample_ID")

Tregs_responder <- RunUMAP(
    Tregs_responder,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

pdf(paste0(savedir, "UMAP/Treg_integrated.pdf"))
DimPlot(Tregs_responder, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

pdf(paste0(savedir, "UMAP/Treg_integrated_samples.pdf"))
DimPlot(Tregs_responder, group.by = "cell_to_sample_ID", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

Tregs_responder <- FindNeighbors(Tregs_responder, reduction = "harmony", dims = 1:30)
Tregs_responder <- FindClusters(Tregs_responder, resolution = 0.3)

pdf(paste0(savedir, "UMAP/celltype_integrated_clusters.pdf"))
DimPlot(Tregs_responder, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

markers <- FindAllMarkers(Tregs_responder)

dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(markers, paste0(savedir, "Table/Tregs_responder_markers.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

genes = c("FOXP3", "IKZF2", "TIGIT", "CTLA4", "ICOS", "HLA-DRB1", "TNFRSF9",
"TCF7", "LEF1", "CXCR6", "CXCR3", "IL10", "TBX21", "RORC",
"BCL6", "PDCD1", "CXCR5", "IL21R", "CCR7", "IL2RA", "FCRL3","CXCL13" )

dir.create(paste0(savedir,"vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/gene_vlnplots.pdf"), width = 12, height = 15)
VlnPlot(Tregs_responder, genes, pt.size = 0)
dev.off()

dir.create(paste0(savedir,"saveRDS_obj"), showWarnings = FALSE)
saveRDS(Tregs_responder, paste0(savedir,"saveRDS_obj/Tregs_responder.RDS"))

#region Non-Responders
# library(Seurat)
# savedir = "/mnt/data/projects/NatMed_HCC/analysis/Tregs/removed_clus12_14/"
# setwd(savedir)
# Tregs = readRDS(paste0(savedir,"saveRDS/Treg_removed_7_8.RDS"))

Tregs_nonresponder = Tregs[,grep("antiPD1_NR",Tregs@meta.data$treatment_Resp)]

Tregs_nonresponder <- NormalizeData(Tregs_nonresponder, normalization.method = "LogNormalize", scale.factor = 10000)
Tregs_nonresponder <- FindVariableFeatures(Tregs_nonresponder, selection.method = "vst", nfeatures = 3000)
Tregs_nonresponder <- ScaleData(Tregs_nonresponder)
Tregs_nonresponder <- RunPCA(Tregs_nonresponder, features = VariableFeatures(object = Tregs_nonresponder))
Tregs_nonresponder <- FindNeighbors(Tregs_nonresponder, dims = 1:30)
Tregs_nonresponder <- RunUMAP(Tregs_nonresponder, dims = 1:30)

savedir = "/mnt/data/projects/NatMed_HCC/analysis/Tregs/removed_clus12_14/non_responder/"
dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE, recursive = TRUE)
setwd(savedir)
pdf(paste0(savedir, "UMAP/Treg_responder_unintegrated.pdf"))
DimPlot(Tregs_nonresponder)
dev.off()

### Integration using Harmony
library(harmony)
Tregs_nonresponder <- RunHarmony(Tregs_nonresponder, group.by.vars = "cell_to_sample_ID")

Tregs_nonresponder <- RunUMAP(
    Tregs_nonresponder,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

pdf(paste0(savedir, "UMAP/Treg_integrated.pdf"))
DimPlot(Tregs_nonresponder, group.by = "celltype", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

pdf(paste0(savedir, "UMAP/Treg_integrated_samples.pdf"))
DimPlot(Tregs_nonresponder, group.by = "cell_to_sample_ID", reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

Tregs_nonresponder <- FindNeighbors(Tregs_nonresponder, reduction = "harmony", dims = 1:30)
Tregs_nonresponder <- FindClusters(Tregs_nonresponder, resolution = 0.3)

pdf(paste0(savedir, "UMAP/Tregs_non_responders_celltype_integrated_clusters.pdf"))
DimPlot(Tregs_nonresponder, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

markers <- FindAllMarkers(Tregs_nonresponder)

dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(markers, paste0(savedir, "Table/Tregs_nonresponder_markers.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

genes = c("FOXP3", "IKZF2", "TIGIT", "CTLA4", "ICOS", "HLA-DRB1", "TNFRSF9",
"TCF7", "LEF1", "CXCR6", "CXCR3", "IL10", "TBX21", "RORC",
"BCL6", "PDCD1", "CXCR5", "IL21R", "CCR7", "IL2RA", "FCRL3","CXCL13" )

dir.create(paste0(savedir,"vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/non_responder_gene_vlnplots.pdf"), width = 12, height = 15)
VlnPlot(Tregs_nonresponder, genes, pt.size = 0)
dev.off()

dir.create(paste0(savedir,"saveRDS_obj"), showWarnings = FALSE)
saveRDS(Tregs_nonresponder, paste0(savedir,"saveRDS_obj/Tregs_nonresponder.RDS"))
