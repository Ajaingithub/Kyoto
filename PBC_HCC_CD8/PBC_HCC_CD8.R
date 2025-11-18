library(Seurat)
library(harmony)
PBC_obj <- readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/CD8/analysis/saveRDS/PBC_CD8_Tcells.RDS")
HCC_obj <- readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/CD8_3/saveRDS/HCC_CD8.RDS")

HCC_obj@meta.data$orig.ident <- paste0("Sample_", HCC_obj@meta.data$cell_to_sample_ID)
HCC_obj@assays$RNA <- HCC_obj@assays$originalexp
# DefaultAssay(HCC_obj) <- "RNA"
# HCC_obj <- DietSeurat(HCC_obj, assays = "RNA", layers = "counts")
# # HCC_obj <- UpdateSeuratObject(HCC_obj) ## from V4 to V5
# layers_to_remove <- setdiff(Layers(HCC_obj[["RNA"]]), "counts")
# for (layer in layers_to_remove) {
#     HCC_obj[["RNA"]]@layers[[layer]] <- NULL
# }

# PBC_obj <- DietSeurat(PBC_obj, assays = "RNA", layers = "counts")
# layers_to_remove <- setdiff(Layers(PBC_obj[["RNA"]]), "counts")
# for (layer in layers_to_remove) {
#     PBC_obj[["RNA"]]@layers[[layer]] <- NULL
# }

#### Combining by counts data directly
library(Seurat)
counts_HCC <- GetAssayData(HCC_obj, assay = "RNA", slot = "counts") # v4 object
counts_PBC <- GetAssayData(PBC_obj, assay = "RNA", layer = "counts") # v5 object

common_genes <- intersect(rownames(counts_HCC), rownames(counts_PBC))
counts_HCC_common <- counts_HCC[common_genes, ]
counts_PBC_common <- counts_PBC[common_genes, ]
combined_counts <- cbind(counts_HCC_common, counts_PBC_common)

# Get all unique columns
meta_HCC <- HCC_obj@meta.data
meta_PBC <- PBC_obj@meta.data
all_cols <- union(colnames(meta_HCC), colnames(meta_PBC))

# Add missing columns with NA
for (col in setdiff(all_cols, colnames(meta_HCC))) {
    meta_HCC[[col]] <- NA
}
for (col in setdiff(all_cols, colnames(meta_PBC))) {
    meta_PBC[[col]] <- NA
}

# Reorder columns to match
meta_HCC <- meta_HCC[, all_cols]
meta_PBC <- meta_PBC[, all_cols]

# Add disease label
meta_HCC$disease <- "HCC"
meta_PBC$disease <- "PBC"

combined_meta <- rbind(meta_HCC, meta_PBC)
all(colnames(combined_counts) == rownames(combined_meta))

HCC_obj <- CreateSeuratObject(counts_HCC_common, meta.data = meta_HCC)
PBC_obj <- CreateSeuratObject(counts_PBC_common, meta.data = meta_PBC)

# Step 5: Create a new Seurat object
PBC_HCC <- CreateSeuratObject(counts = combined_counts, meta.data = combined_meta)

# PBC_HCC <- readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/saveRDS/PBC_HCC.RDS")
PBC_HCC[["RNA"]] <- split(PBC_HCC[["RNA"]], f = PBC_HCC$disease)

PBC_HCC <- NormalizeData(PBC_HCC)
PBC_HCC <- FindVariableFeatures(PBC_HCC)
PBC_HCC <- ScaleData(PBC_HCC)
PBC_HCC <- RunPCA(PBC_HCC)
### Without Integration
# PBC_HCC <- RunUMAP(PBC_HCC, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC_CD8/"
# dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE)
# pdf(paste0(savedir, "UMAP/disease_unitegrated.pdf"), width = 6, height = 5)
# DimPlot(PBC_HCC, reduction = "umap.unintegrated", group.by = c("disease"))
# dev.off()

# PBC_HCC <- RunHarmony(PBC_HCC, group.by.vars = "disease")

# PBC_HCC <- RunUMAP(
#     PBC_HCC,
#     reduction.key = "harmonyUMAP_",
#     reduction = "harmony",
#     reduction.name = "harmonyumap",
#     dims = 1:30
# )

# pdf(paste0(savedir, "UMAP/integrated_umap.pdf"))
# DimPlot(PBC_HCC, reduction = "harmonyumap", group.by = c("disease"))
# dev.off()

# savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/"
# pdf(paste0(savedir, "UMAP/harmony_umap_splitted_celltypes_group.pdf"))
# DimPlot(PBC_HCC, reduction = "harmonyumap", group.by = c("celltypes"), split.by = "disease")
# dev.off()

# pdf(paste0(savedir, "UMAP/integrated_umap_celltypes_group.pdf"))
# DimPlot(PBC_HCC, reduction = "harmonyumap", group.by = c("celltype"))
# dev.off()

s
### Integrated CCA
PBC_HCC[["RNA"]] <- split(PBC_HCC[["RNA"]], f = PBC_HCC$disease)

PBC_HCC <- NormalizeData(PBC_HCC)
PBC_HCC <- FindVariableFeatures(PBC_HCC)
PBC_HCC <- ScaleData(PBC_HCC)
PBC_HCC <- RunPCA(PBC_HCC)


PBC_HCC <- IntegrateLayers(
    object = PBC_HCC, method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    verbose = FALSE
)

# re-join layers after integration
PBC_HCC[["RNA"]] <- JoinLayers(PBC_HCC[["RNA"]])
PBC_HCC <- RunUMAP(PBC_HCC, dims = 1:30, reduction = "integrated.rpca")

# Visualization
savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC_CD8/"
pdf(paste0(savedir, "UMAP/RPCA_integrated_umap.pdf"))
DimPlot(PBC_HCC, reduction = "umap", group.by = c("disease"))
dev.off()

pdf(paste0(savedir, "UMAP/RPCA_integrated_umap_splitted_celltypes_group.pdf"))
DimPlot(PBC_HCC, reduction = "umap", group.by = c("celltype"), split.by = "disease")
dev.off()

setwd(savedir)
saveRDS(PBC_HCC, "PBC_HCC_CD8.RDS")

PBC_HCC <- IntegrateLayers(
    object = PBC_HCC, method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE
)

# re-join layers after integration
PBC_HCC[["RNA"]] <- JoinLayers(PBC_HCC[["RNA"]])
PBC_HCC <- RunUMAP(
    PBC_HCC,
    dims = 1:30,
    reduction = "integrated.cca",
    reduction.name = "cca.umap"
)

pdf(paste0(savedir, "UMAP/CCA_integrated_umap.pdf"))
DimPlot(PBC_HCC, reduction = "cca.umap", group.by = c("seurat_clusters"), label = TRUE, label.size = 6)
dev.off()

pdf(paste0(savedir, "UMAP/CCA_integrated_umap_celltypes_group_labelled.pdf"), width = 10, height = 6)
DimPlot(PBC_HCC, reduction = "cca.umap", group.by = c("celltype"), label = TRUE)
dev.off()

pdf(paste0(savedir, "UMAP/CCA_integrated_umap_celltypes_group_label.pdf"), width = 14, height = 6)
DimPlot(PBC_HCC, reduction = "cca.umap", group.by = c("celltype"), split.by = "disease", label = TRUE)
dev.off()

PBC_HCC <- FindNeighbors(PBC_HCC, reduction = "integrated.cca", dims = 1:30)
PBC_HCC <- FindClusters(PBC_HCC, resolution = 0.4)

pdf(paste0(savedir, "UMAP/CCA_integrated_umap_cluster.pdf"))
DimPlot(PBC_HCC, reduction = "cca.umap", group.by = c("seurat_clusters"), label = TRUE, label.size = 6)
dev.off()

pdf(paste0(savedir, "UMAP/CCA_integrated_umap_cluster_group_label.pdf"), width = 14, height = 6)
DimPlot(PBC_HCC, reduction = "cca.umap", group.by = c("seurat_clusters"), split.by = "disease", label = TRUE)
dev.off()

saveRDS(PBC_HCC, "PBC_HCC_CD8_CCA.RDS")

all_markers <- FindAllMarkers(PBC_HCC)
dir.create(paste0(savedir, "Table"))
write.table(all_markers, paste0(savedir, "Table/all_marker.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

genes <- c(
    "TCF7", "LEF1", "IL7R", "GPR183", "LTB", "CD28", "ICOS", "MAF", "TNFSF8", "CXCL13", "IL21",
    "CSF2", "IFNG", "TNF", "FOXP3", "IL2RA", "ENTPD1", "CTLA4", "PDCD1", "TIGIT", "HAVCR2", "GZMK",
    "CXCR5", "CXCR6", "CXCR3", "CCR7"
)

dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir, "vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/PBC_HCC_genes.pdf"), width = 12, height = 12)
VlnPlot(PBC_HCC, genes, pt.size = 0)
dev.off()
