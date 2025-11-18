library(Seurat)
library(harmony)
PBC_obj <- readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/CD4/saveRDS/CD4_obj.RDS")
HCC_obj <- readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/NatMed_HCC/analysis/saveRDS/CD4_5.RDS")

HCC_obj@meta.data$orig.ident <- paste0("Sample_", HCC_obj@meta.data$cell_to_sample_ID)
HCC_obj@assays$RNA <- HCC_obj@assays$originalexp
HCC_obj <- DietSeurat(HCC_obj, assays = "RNA", layers = "counts")
# HCC_obj <- UpdateSeuratObject(HCC_obj) ## from V4 to V5
layers_to_remove <- setdiff(Layers(HCC_obj[["RNA"]]), "counts")
for (layer in layers_to_remove) {
    HCC_obj[["RNA"]]@layers[[layer]] <- NULL
}

PBC_obj <- DietSeurat(PBC_obj, assays = "RNA", layers = "counts")
layers_to_remove <- setdiff(Layers(PBC_obj[["RNA"]]), "counts")
for (layer in layers_to_remove) {
    PBC_obj[["RNA"]]@layers[[layer]] <- NULL
}

#### Combining by counts data directly
library(Seurat)
counts_HCC <- GetAssayData(HCC_obj, assay = "RNA", slot = "counts") # v4 object
counts_PBC <- GetAssayData(PBC_obj, assay = "RNA", layer = "counts") # v5 object

common_genes <- intersect(rownames(counts_HCC), rownames(counts_PBC))
counts_HCC_common <- counts_HCC[common_genes, ]
counts_PBC_common <- counts_PBC[common_genes, ]

# combined_counts <- cbind(counts_HCC_common, counts_PBC_common)
# Get all unique columns
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
# colnames(combined_counts) == rownames(combined_meta)

HCC_obj <- CreateSeuratObject(counts_HCC_common, meta.data = meta_HCC)
PBC_obj <- CreateSeuratObject(counts_PBC_common, meta.data = meta_PBC)

# Step 5: Create a new Seurat object
PBC_HCC <- CreateSeuratObject(counts = combined_counts, meta.data = combined_meta)

PBC_HCC <- readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/saveRDS/PBC_HCC.RDS")
PBC_HCC[["RNA"]] <- split(PBC_HCC[["RNA"]], f = PBC_HCC$disease)

PBC_HCC <- NormalizeData(PBC_HCC)
PBC_HCC <- FindVariableFeatures(PBC_HCC)
PBC_HCC <- ScaleData(PBC_HCC)
PBC_HCC <- RunPCA(PBC_HCC)
### Without Integration
PBC_HCC <- RunUMAP(PBC_HCC, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/"
dir.create(paste0(savedir, "UMAP"), showWarnings = FALSE)
pdf(paste0(savedir, "UMAP/disease_unitegrated.pdf"), width = 6, height = 5)
DimPlot(PBC_HCC, reduction = "umap.unintegrated", group.by = c("disease"))
dev.off()

PBC_HCC <- RunHarmony(PBC_HCC, group.by.vars = "disease")

PBC_HCC <- RunUMAP(
    PBC_HCC,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

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
savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/"
pdf(paste0(savedir, "UMAP/integrated_umap.pdf"))
DimPlot(PBC_HCC, reduction = "umap", group.by = c("disease"))
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/"
pdf(paste0(savedir, "UMAP/integrated_umap_splitted_celltypes_group.pdf"))
DimPlot(PBC_HCC, reduction = "umap", group.by = c("celltypes"), split.by = "disease")
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/"
pdf(paste0(savedir, "UMAP/integrated_umap_celltypes_group.pdf"))
DimPlot(PBC_HCC, reduction = "umap", group.by = c("celltypes"))
dev.off()

### Integrated CCA
PBC_HCC <- IntegrateLayers(
    object = PBC_HCC, method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE
)

# re-join layers after integration
PBC_HCC[["RNA"]] <- JoinLayers(PBC_HCC[["RNA"]])
PBC_HCC <- RunUMAP(PBC_HCC, dims = 1:30, reduction = "integrated.cca")

# Visualization
savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/"
pdf(paste0(savedir, "UMAP/CCA_integrated_umap.pdf"))
DimPlot(PBC_HCC, reduction = "umap", group.by = c("disease"))
dev.off()

pdf(paste0(savedir, "UMAP/CCA_integrated_umap_splitted_celltypes_group.pdf"))
DimPlot(PBC_HCC, reduction = "umap", group.by = c("celltypes"), split.by = "disease")
dev.off()

pdf(paste0(savedir, "UMAP/CCA_integrated_umap.pdf"))
DimPlot(PBC_HCC_req, reduction = "umap", group.by = c("seurat_clusters"), label = TRUE, label.size = 6)
dev.off()

pdf(paste0(savedir, "UMAP/CCA_integrated_umap_celltypes_group_labelled.pdf"), width = 10, height = 6)
DimPlot(PBC_HCC, reduction = "umap", group.by = c("celltypes"), label = TRUE)
dev.off()

pdf(paste0(savedir, "UMAP/CCA_integrated_umap_celltypes_group_label.pdf"), width = 14, height = 6)
DimPlot(PBC_HCC, reduction = "umap", group.by = c("celltypes"), split.by = "disease", label = TRUE)
dev.off()

PBC_HCC <- FindNeighbors(PBC_HCC, reduction = "integrated.cca", dims = 1:30)
PBC_HCC <- FindClusters(PBC_HCC, resolution = 1)

PBC_HCC <- FindClusters(PBC_HCC, resolution = 0.3, clus)

pdf(paste0(savedir, "UMAP/CCA_integrated_seurat_cluster_2.pdf"))
DimPlot(PBC_HCC, reduction = "umap", split.by = "disease", label = TRUE)
dev.off()

Idents(PBC_HCC) <- "seurat_clusters"
stem_marker_conserved <- FindConservedMarkers(PBC_HCC, ident.1 = c("4", "6", "10", "18"), grouping.var = "disease", verbose = FALSE)

dir.create(paste0(savedir, "Table"), showWarnings = FALSE)
write.table(stem_marker_conserved, paste0(savedir, "Table/stem_marker_conserved_PBC_HCC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

# PBC_HCC$disease <- factor(PBC_HCC$disease, levels = c("PBC", "HCC"))
# Idents(PBC_HCC) <- "seurat_clusters"
# stem_marker <- FindMarkers(PBC_HCC, ident.1 = c("4", "6", "10", "18"), grouping.var = "disease", verbose = FALSE)
# write.table(stem_marker, paste0(savedir, "Table/stem_marker_diff_PBC_vs_HCC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

# PBC_HCC@meta.data$celltypes <- gsub("Stem-like CD4 T cells", "Stem-like CD4 T cell", PBC_HCC@meta.data$celltypes)

# Idents(PBC_HCC) <- "celltypes"
# stem_marker <- FindMarkers(PBC_HCC, ident.1 = c("Stem-like CD4 T cell"), grouping.var = "disease", verbose = FALSE)
# write.table(stem_marker, paste0(savedir, "Table/Stem_like_CD4_T_cell_marker_diff_PBC_vs_HCC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

# GZMK_Th1 <- FindMarkers(PBC_HCC, ident.1 = c("GZMK Th1"), grouping.var = "disease", verbose = TRUE)
# write.table(GZMK_Th1, paste0(savedir, "Table/GZMK_Th1_marker_diff_PBC_vs_HCC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

# ## Confirming
# Th1_cellnames <- rownames(PBC_HCC@meta.data[grep("GZMK Th1", PBC_HCC@meta.data$celltypes), ])
# PBC_HCC_Th1 <- subset(PBC_HCC, cells = Th1_cellnames)

# stem_marker_genes <- c("DUSP4", "CST7", "ID2", "CCR7", "TCF7", "SELL")
# pdf(paste0(savedir, "vlnplots/SM_marker_disease_nopt.pdf"))
# VlnPlot(PBC_HCC_SM, pt.size = 0, features = dual_producer_genes, group.by = "disease", ncol = 3)
# dev.off()

# Tph <- FindMarkers(PBC_HCC, ident.1 = c("Tph"), grouping.var = "disease", verbose = TRUE)
# write.table(Tph, paste0(savedir, "Table/Tph_marker_diff_PBC_vs_HCC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

# PBC_HCC@meta.data$celltypes <- gsub("Dual producing Th1", "Dual producer", PBC_HCC@meta.data$celltypes)
# Dual_producer <- FindMarkers(PBC_HCC, ident.1 = c("Dual producer"), grouping.var = "disease", verbose = TRUE)
# write.table(Dual_producer, paste0(savedir, "Table/Dual_producer_marker_diff_PBC_vs_HCC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

### The above differential does not work as we expected it is performing differential between the cluster not within
SM_cellnames <- rownames(PBC_HCC@meta.data[grep("Stem-like CD4 T cell", PBC_HCC@meta.data$celltypes), ])
PBC_HCC_SM <- subset(PBC_HCC, cells = SM_cellnames)

PBC_HCC_SM <- NormalizeData(PBC_HCC_SM)
SM_marker_subset <- FindMarkers(PBC_HCC_SM, ident.1 = "PBC", ident.2 = "HCC", group.by = "disease")
write.table(SM_marker_subset, paste0(savedir, "Table/SM_marker_subset_marker_diff_PBC_vs_HCC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

stem_marker_genes <- c("TCF7", "LEF1", "BACH2", "BCL11B", "KLF12", "ZBTB20", "FOXN3", "TNFSF8", "PDCD1", "TIGIT")
pdf(paste0(savedir, "vlnplots/SM_marker_subset_disease_nopt.pdf"), height = 7)
VlnPlot(PBC_HCC_SM, pt.size = 0, features = stem_marker_genes, group.by = "disease", ncol = 4)
dev.off()

library(dplyr)
library(ggplot2)
library(ggrepel)
library(Rmpfr)
savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/"
PBC_HCC_SM <- read.table(paste(savedir, "Table/SM_marker_subset_marker_diff_PBC_vs_HCC.txt", sep = ""), sep = "\t")

tops_plot <- select(PBC_HCC_SM, avg_log2FC, p_val_adj, p_val)
p_floor <- 1e-320 # or e.g. 1e-320, just above double precision limit
tops_plot$neg_log10_pval_adj <- -log10(pmax(tops_plot$p_val_adj, p_floor))

tops_plot[, "-log10(adj_pvalue)"] <- -log10(tops_plot$p_val_adj)
tops_plot <- mutate(tops_plot,
    sig = ifelse(tops_plot[, "p_val_adj"] < 0.01, ### adjusted pvalue < 0.01
        "Sig", "Not Sig"
    )
)

tops_plot$significant <- tops_plot$sig
tops_plot_sig <- tops_plot[grep("^Sig$", tops_plot$significant), ]
tops_plot_sig_up <- tops_plot_sig[which(tops_plot_sig$avg_log2FC > 0), ]
tops_plot_sig_down <- tops_plot_sig[which(tops_plot_sig$avg_log2FC < 0), ]

tops_plot[which(rownames(tops_plot) %in% rownames(tops_plot_sig_up)), "significant"] <- "Up_sig"
tops_plot[which(rownames(tops_plot) %in% rownames(tops_plot_sig_down)), "significant"] <- "Down_sig"

tops_plot$significant <- gsub("^Sig$", "Not Sig", tops_plot$significant)
volc <- ggplot(tops_plot, aes(avg_log2FC, neg_log10_pval_adj)) +
    geom_point(aes(col = significant), size = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    # geom_vline(xintercept = -2, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("blue", "grey", "red")) +
    annotate("text", x = 7, y = 300, label = paste0("Up=", nrow(tops_plot_sig_up))) +
    annotate("text", x = -7, y = 300, label = paste0("Down=", nrow(tops_plot_sig_down))) +
    # ggtitle(colnames(cm)) +
    theme(
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black")
    )
# panel.grid.minor = element_line(colour = "grey"),
# panel.grid.major = element_line(colour = "grey"))
tops_plot$gene <- rownames(tops_plot)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/Table/"
pdf(paste(savedir, "SM_PBC_vs_HCC_Volcano_plots.pdf", sep = ""), width = 5.5, height = 8)
volc # + geom_text_repel(data = head(tops_plot[order(tops_plot$p_val_adj), ], 100), aes(label = gene))
dev.off()


### TpH
Tph_cellnames <- rownames(PBC_HCC@meta.data[grep("Tph", PBC_HCC@meta.data$celltypes), ])
PBC_HCC_Tph <- subset(PBC_HCC, cells = Tph_cellnames)

Tph_marker_subset <- FindMarkers(PBC_HCC_Tph, ident.1 = "PBC", ident.2 = "HCC", group.by = "disease")
write.table(Tph_marker_subset, paste0(savedir, "Table/Tph_marker_subset_marker_diff_PBC_vs_HCC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

### confirming
Tph_genes <- c(
    "ICOS", "MAF", "RORA", "TOX", "TOX2", "PDCD1", "LAG3", "TIGIT", "HAVCR2", "IL21",
    "CXCL13", "TCF7", "TNFSF8"
)
pdf(paste0(savedir, "vlnplots/Tph_marker_subset_disease_nopt.pdf"), height = 10)
VlnPlot(PBC_HCC_Tph, pt.size = 0, features = Tph_genes, group.by = "disease", ncol = 4)
dev.off()

### DP Testing
DP_cellnames <- rownames(PBC_HCC@meta.data[grep("Dual producer", PBC_HCC@meta.data$celltypes), ])
PBC_HCC_DP <- subset(PBC_HCC, cells = DP_cellnames)

DP_marker_subset <- FindMarkers(PBC_HCC_DP, ident.1 = "PBC", ident.2 = "HCC", group.by = "disease")
write.table(DP_marker_subset, paste0(savedir, "Table/DP_marker_subset_marker_diff_PBC_vs_HCC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

### confirming
dual_producer_genes <- c("HSPA1B", "DNAJB1", "IFNG", "TNF", "CCL5", "PDCD1")
pdf(paste0(savedir, "vlnplots/DP_marker_subset_disease_nopt.pdf"), height = 4.5)
VlnPlot(PBC_HCC_DP, pt.size = 0, features = dual_producer_genes, group.by = "disease", ncol = 4)
dev.off()

# GZMK_Th1 <- FindMarkers(PBC_HCC, ident.1 = c("GZMK Th1"), grouping.var = "disease", verbose = TRUE)
# write.table(GZMK_Th1, paste0(savedir, "Table/GZMK_Th1_marker_diff_PBC_vs_HCC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

## GZMK Th1
Th1_cellnames <- rownames(PBC_HCC@meta.data[grep("GZMK Th1", PBC_HCC@meta.data$celltypes), ])
PBC_HCC_Th1 <- subset(PBC_HCC, cells = Th1_cellnames)

Th1_marker_subset <- FindMarkers(PBC_HCC_Th1, ident.1 = "PBC", ident.2 = "HCC", group.by = "disease")
write.table(Th1_marker_subset, paste0(savedir, "Table/Th1_marker_subset_marker_diff_PBC_vs_HCC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

### confirming
Th1_producer_genes <- c("GZMK", "IFNG", "TNF", "TOX", "PDCD1", "TIGIT", "LAG3")
pdf(paste0(savedir, "vlnplots/Th1_marker_subset_disease_nopt.pdf"), height = 4.5)
VlnPlot(PBC_HCC_Th1, pt.size = 0, features = Th1_producer_genes, group.by = "disease", ncol = 4)
dev.off()

# Idents(PBC_HCC) <- "res"
PBC_vs_HCC <- FindMarkers(PBC_HCC, ident.1 = "PBC", ident.2 = "HCC", group.by = "disease", verbose = FALSE)
write.table(PBC_vs_HCC, paste0(savedir, "Table/marker_diff_PBC_vs_HCC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

marker <- FindAllMarkers(PBC_HCC)
write.table(marker, paste0(savedir, "Table/cluster_markers.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

PBC_HCC <- FindClusters(PBC_HCC, resolution = 0.4, cluster.name = "res0.4")
Idents(PBC_HCC) <- "res0.4"
pdf(paste0(savedir, "UMAP/CCA_integrated_seurat_cluster_res0.4.pdf"), width = 14, height = 6)
DimPlot(PBC_HCC, reduction = "umap", split.by = "disease", label = TRUE)
dev.off()

PBC_HCC_req <- FindClusters(PBC_HCC_req, resolution = 0.3, cluster.name = "res0.3")
pdf(paste0(savedir, "UMAP/CCA_integrated_seurat_cluster_res0.3_0.25.pdf"), width = 8, height = 6)
DimPlot(PBC_HCC_req, reduction = "umap", group.by = c("res0.3", "res0.25"), label = TRUE, label.size = 6)
dev.off()

marker <- FindAllMarkers(PBC_HCC)
write.table(marker, paste0(savedir, "Table/cluster_markers_res0.4.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

### removing cluster 10 to 13 due to contamination from resolution 0.4
req_cells <- rownames(PBC_HCC@meta.data[grep("10|11|12|13", PBC_HCC@meta.data$res0.4, invert = TRUE), ])
PBC_HCC_req <- subset(PBC_HCC, cells = req_cells)

pdf(paste0(savedir, "UMAP/CCA_integrated_seurat_cluster_res0.4_removed_clus_2.pdf"), width = 10, height = 5)
DimPlot(PBC_HCC_req, reduction = "umap", split.by = "disease", label = TRUE, group.by = "res0.4")
dev.off()

pdf(paste0(savedir, "UMAP/CCA_integrated_seurat_cluster_res0.4_removed_clus_nosplit.pdf"), height = 5, width = 5)
DimPlot(PBC_HCC_req, reduction = "umap", label = TRUE, group.by = "res0.4")
dev.off()

pdf(paste0(savedir, "UMAP/CCA_integrated_seurat_cluster_res0.4_removed_clus_nosplit.pdf"), height = 5, width = 5)
DimPlot(PBC_HCC_req, reduction = "umap", label = TRUE, group.by = "res0.4")
dev.off()

PBC_cells <- rownames(PBC_HCC_req@meta.data[grep("PBC", PBC_HCC_req@meta.data$disease), ])
HCC_cells <- rownames(PBC_HCC_req@meta.data[grep("HCC", PBC_HCC_req@meta.data$disease), ])

PBC_req_obj <- subset(PBC_HCC_req, cells = PBC_cells)
HCC_req_obj <- subset(PBC_HCC_req, cells = HCC_cells)

pdf(paste0(savedir, "UMAP/PBC_req_obj.pdf"), height = 5, width = 5)
DimPlot(PBC_req_obj, reduction = "umap", label = TRUE, group.by = "res0.4")
dev.off()

pdf(paste0(savedir, "UMAP/HCC_req_obj.pdf"), height = 5, width = 5)
DimPlot(HCC_req_obj, reduction = "umap", label = TRUE, group.by = "res0.4")
dev.off()


write.table(table(PBC_HCC_req@meta.data$res0.4, PBC_HCC_req@meta.data$disease),
    paste0(savedir, "Table/PBC_HCC_cluster_disease.txt"),
    sep = "\t", row.names = T, col.names = T, quote = F
)

saveRDS(PBC_HCC_req, paste0(savedir, "saveRDS/PBC_HCC_CCA.RDS"))

Idents(PBC_HCC) <- (PBC_HCC@meta.data$res0.4)
PBC_HCC <- NormalizeData(PBC_HCC)
PBC_HCC_markers <- FindAllMarkers(PBC_HCC)
savedir = "/mnt/data/projects/PBC_HCC/"
write.table(PBC_HCC_markers, paste0(savedir,"Table/PBC_HCC_markers_res0.4_no_clus10_to_13.txt"),  quote = F, col.names = T, row.names = T, sep = "\t")

#### Violin plots
Set1 <- c("TCF7", "LEF1", "SELL", "CCR7", "FOXP3", "TIGIT")
Set2 <- c("CTLA4", "IFNG", "TNF", "CCL5", "CD40LG", "GZMK")
Set3 <- c("HSPA1A", "HSPA1B", "HSPD1", "NR4A1", "HSPA6", "BAG3")
Set4 <- c("CXCL13", "PDCD1", "TNFSF8", "IL21", "TOX2", "MAF", "GPR183")
Set5 <- c("ZEB2", "ID2", "XCL1", "XCL2", "CCL4", "GZMA")

sets <- c("Set1","Set2","Set3","Set4","Set5")
savedir = "/mnt/data/projects/PBC_HCC/"

for(i in 1:length(sets)){
    # pdf(paste0(savedir,"vlnplots/",sets[i],"_PBC_HCC_genes.pdf"))
    # print(VlnPlot(PBC_HCC,get(sets[i])))
    # dev.off()
    pdf(paste0(savedir,"vlnplots/",sets[i],"_PBC_HCC_genes_no_points.pdf"))
    print(VlnPlot(PBC_HCC,get(sets[i]), pt.size = 0))
    dev.off()
}

# region CD4 Trajectory
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

PBC_HCC <- readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/saveRDS/PBC_HCC_CCA.RDS")

req_cells <- rownames(PBC_HCC@meta.data[grep("10|11|12|13", PBC_HCC@meta.data$res0.4, invert = TRUE), ])
PBC_HCC_req <- subset(PBC_HCC, cells = req_cells)

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
# Save the objects as separate matrices for input in slingshot
## convert back to singleCellExperiment
# Random subsample, e.g. 5000 cells
set.seed(123)
subset_cells <- sample(rownames(PBC_HCC_req@meta.data), 25000)
PBC_HCC_req2 <- subset(PBC_HCC_req, cells = subset_cells)

DefaultAssay(PBC_HCC_req2) <- "RNA"
sce <- as.SingleCellExperiment(PBC_HCC_req2, assay = "RNA")


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
savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/Trajectory/"
cell_embedding <- PBC_HCC_req2@reductions$umap@cell.embeddings
cell_embedding <- as.data.frame(cell_embedding)
cell_embedding$seurat_clusters <- PBC_HCC_req2@meta.data$seurat_clusters
cell_embedding_6clus <- cell_embedding[grep("6", cell_embedding$seurat_clusters), ]
remove_cells <- rownames(cell_embedding_6clus[cell_embedding_6clus$umap_1 > 0, ])
remain_cells <- grep(paste0(remove_cells, collapse = "|"), rownames(PBC_HCC_req2@meta.data), invert = TRUE, value = TRUE)
PBC_HCC_req3 <- subset(PBC_HCC_req2, cells = remain_cells)

dir.create(savedir, showWarnings = FALSE)
dir.create(paste(savedir, "UMAP/", sep = ""), showWarnings = FALSE)
p <- DimPlot(PBC_HCC_req2, reduction = "umap", label = TRUE, group.by = "seurat_clusters")
pdf(paste(savedir, "UMAP/combine.pdf", sep = ""))
p
dev.off()

p <- DimPlot(PBC_HCC_req2, reduction = "umap", label = TRUE, group.by = "res0.8")
pdf(paste(savedir, "UMAP/res0.8.pdf", sep = ""))
p
dev.off()

## highlight cells
clus6 <- WhichCells(PBC_HCC_req3, idents = c("6"))
pdf(paste(savedir, "UMAP/cluster6_UMAP_2.pdf", sep = ""))
DimPlot(PBC_HCC_req3, label = T, group.by = "seurat_clusters", reduction = "umap", cells.highlight = clus6, cols.highlight = c("darkred"), cols = "grey")
dev.off()

# CD4_sce <- slingshot(sce,
#     clusterLabels = "seurat_clusters2",
#     reducedDim = "UMAP", approx_points = 100,
#     omega = TRUE,
#     omega_scale = 1.5
# )

DefaultAssay(PBC_HCC_req3) <- "RNA"
sce <- as.SingleCellExperiment(PBC_HCC_req3, assay = "RNA")

CD4_sce_naive <- slingshot(sce,
    clusterLabels = colData(sce)$res0.4,
    # clusterLabels = colData(sce)$celltypes,
    reducedDim = "UMAP",
    start.clus = c("7"),
    # start.clus = c("Stem-like CD4 T cells"),
    # end.clus = c("7"),
    approx_points = 150
    # omega = TRUE,omega_scale=1.5
)

# library(ArchR)
# library(scater)
# embedded_orig <- embedCurves(CD4_sce, "UMAP")
# rm(plot_list)
# plot_list <- list()

# for (i in 1:8) {
#   embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
#   embedded <- data.frame(embedded$s[embedded$ord, ])
#   g <- plotUMAP(CD4_sce, colour_by = paste("slingPseudotime_", i, sep = ""))
#   stopifnot(all(rownames(combine@reductions$umap@cell.embeddings) == rownames(g$data)))
#   data <- merge(combine@reductions$umap@cell.embeddings, g$data, by = "row.names")
#   colnames(data) <- c("cellname", "UMAP_1", "UMAP_2", "X", "Y", paste("Lineage", i, sep = ""))
#   p <- ggplot(data, aes_string("UMAP_1", "UMAP_2", color = paste("Lineage", i, sep = ""))) +
#     geom_point(size = 0.01) +
#     scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
#     geom_path(data = embedded, aes(x = UMAP_1, y = UMAP_2), color = "black", size = 1.2) +
#     theme_bw()
#   plot_list[[i]] <- p
# }

# require(gridExtra)
# pdf(paste(savedir, "UMAP/sce_UMAP_splitted_unbias.pdf", sep = ""), width = 12, height = 9)
# grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]], nrow = 3, ncol = 3)
# dev.off()

library(ArchR)
library(scater)
embedded_orig <- embedCurves(CD4_sce_naive, "UMAP")

rm(plot_list)
plot_list <- list()
for (i in 1:4) {
    embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
    embedded <- data.frame(embedded$s[embedded$ord, ])
    g <- plotUMAP(CD4_sce_naive, colour_by = paste("slingPseudotime_", i, sep = ""), dimred = "UMAP", text_by = "celltypes")
    stopifnot(all(rownames(PBC_HCC_req3@reductions$umap@cell.embeddings) == rownames(g$data)))
    # data <- merge(embedded_orig@reductions$umap@cell.embeddings, g$data, by = "row.names")
    # colnames(data) <- c("cellname", "HARMONYUMAP_1", "HARMONYUMAP_2", "X", "Y", paste("Lineage", i, sep = ""))
    colnames(g$data) <- c("UMAP_1", "UMAP_2", paste("Lineage", i, sep = ""))
    p <- ggplot(g$data, aes_string("UMAP_1", "UMAP_2", color = paste("Lineage", i, sep = ""))) +
        geom_point(size = 0.01) +
        scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
        geom_path(data = embedded, aes(x = umap_1, y = umap_2), color = "black", size = 1.2) +
        theme_bw()
    plot_list[[i]] <- p
}


require(gridExtra)
pdf(paste(savedir, "UMAP/PBC_HCC_CD4_cluster7_res0.4_remove_cells.pdf", sep = ""), width = 8.5, height = 8)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], nrow = 2, ncol = 2)
dev.off()

pdf(paste(savedir, "UMAP/PBC_HCC_CD4_UMAP_cluster7_lineage_res0.4_remove_cells.pdf", sep = ""), width = 4.5, height = 4)
plot_list
dev.off()

# Extract plot data
# plot_data <- g$data
# plot_data$celltypes <- colData(CD4_sce_naive)$celltypes
# colnames(plot_data) <- c("UMAP_1", "UMAP_2", paste("Lineage", i, sep = ""), "order", "celltypes")

# # Compute cluster label positions (centroids)
# label_positions <- plot_data %>%
#     group_by(celltypes) %>%
#     summarise(
#         UMAP_1 = median(UMAP_1),
#         UMAP_2 = median(UMAP_2)
#     )

# # Build your ggplot with pseudotime
# p <- ggplot(plot_data, aes_string("UMAP_1", "UMAP_2", color = paste0("Lineage", i))) +
#     geom_point(size = 0.01) +
#     scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
#     geom_path(data = embedded, aes(x = umap_1, y = umap_2), color = "black", size = 1.2) +
#     ggrepel::geom_text_repel(
#         data = label_positions,
#         aes(x = UMAP_1, y = UMAP_2, label = celltypes),
#         size = 3,
#         color = "black"
#     ) +
#     theme_bw()

# plot_list[[i]] <- p

# pdf(paste(savedir, "UMAP/test_label.pdf", sep = ""), width = 4.5, height = 4)
# p
# dev.off()


# Example: Read data from a CSV or TSV file if needed
# gene_data <- read.csv("gene_list.csv", header = TRUE)

# Loop through each column
setwd("/diazlab/data3/.abhinav/.immune/Kyoto/resource/CD8_genes/")
CD8_genes <- read.delim("/diazlab/data3/.abhinav/.immune/Kyoto/resource/CD8_gene_list_pancer_Nat_Med")
for (col in colnames(CD8_genes)) {
    # Get non-empty genes
    genes <- na.omit(CD8_genes[[col]])
    genes <- genes[genes != ""] # Remove empty strings if any

    # Save to file named after the column
    writeLines(genes, paste0(gsub("\\.", "_", col), ".txt"))
}

setwd("/diazlab/data3/.abhinav/.immune/Kyoto/resource/CD4_genes/")
CD4_genes <- read.delim("/diazlab/data3/.abhinav/.immune/Kyoto/resource/CD4_gene_list_pancer_Nat_Med")
for (col in colnames(CD4_genes)) {
    # Get non-empty genes
    genes <- na.omit(CD4_genes[[col]])
    genes <- genes[genes != ""] # Remove empty strings if any

    # Save to file named after the column
    writeLines(genes, paste0(gsub("\\.", "_", col), ".txt"))
}

DefaultAssay(PBC_HCC) <- "RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/.immune/Kyoto/resource/CD4_genes/", pattern = "txt", full.names = TRUE)
filename <- gsub(".txt", "", basename(files))
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[, 1]
    Tcellsubset <- rownames(PBC_HCC)[match(Tcellsubset, rownames(PBC_HCC), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

PBC_HCC <- AddModuleScore(PBC_HCC, features = markers, slot = "data")
colnames(PBC_HCC@meta.data)[44:60] <- paste0("CD4_", filename)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/"
dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir, "vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/CD4_pancancer_module_score_genes.pdf"), width = 12, height = 12)
VlnPlot(PBC_HCC_req, paste0("CD4_", filename), pt.size = 0)
dev.off()

dir.create(paste0(savedir, "featureplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/CD4_pancancer_module_score_genes.pdf"), width = 12, height = 12)
FeaturePlot(PBC_HCC_req, paste0("CD4_", filename), reduction = "umap", label = TRUE)
dev.off()

### CD8
DefaultAssay(PBC_HCC) <- "RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/.immune/Kyoto/resource/CD8_genes/", pattern = "txt", full.names = TRUE)
filename <- gsub(".txt", "", basename(files))
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[, 1]
    Tcellsubset <- rownames(PBC_HCC)[match(Tcellsubset, rownames(PBC_HCC), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

PBC_HCC <- AddModuleScore(PBC_HCC, features = markers, slot = "data")
colnames(PBC_HCC@meta.data)[61:79] <- paste0("CD8_", filename)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/"
dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir, "vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/CD8_pancancer_module_score_genes.pdf"), width = 12, height = 12)
VlnPlot(PBC_HCC_req, paste0("CD8_", filename), pt.size = 0)
dev.off()

dir.create(paste0(savedir, "featureplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/CD8_pancancer_module_score_genes.pdf"), width = 12, height = 12)
FeaturePlot(PBC_HCC_req, paste0("CD8_", filename), reduction = "umap", label = TRUE)
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/"
pdf(paste0(savedir, "vlnplots/CD4_CXCL13_module_score_genes.pdf"), width = 4.4, height = 4.3)
VlnPlot(PBC_HCC_req, "CD4_CXCL13", pt.size = 0)
dev.off()

DefaultAssay(PBC_HCC_req) <- "RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/.immune/Kyoto/liver/resources/genelist/", pattern = "CD", full.names = TRUE)
filename <- paste0(basename(files), "_2")
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[, 1]
    Tcellsubset <- rownames(PBC_HCC_req)[match(Tcellsubset, rownames(PBC_HCC_req), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

PBC_HCC_req <- AddModuleScore(PBC_HCC_req, features = markers, slot = "data")
colnames(PBC_HCC_req@meta.data)[82:92] <- filename

dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir, "vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/Cell_2017_violinplot_PBC_HCC.pdf"), width = 12, height = 9)
VlnPlot(PBC_HCC_req, filename, pt.size = 0, group.by = "res0.4")
dev.off()

pdf(paste0(savedir, "featureplots/Cell_2017_featureplot_PBC_HCC.pdf"), width = 12, height = 12)
FeaturePlot(PBC_HCC_req, filename, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

pdf(paste0(savedir, "UMAP/dimplot_seurat_clusters2.pdf"), width = 4.5, height = 4.5)
DimPlot(PBC_HCC_req2, group.by = "res0.4", reduction = "umap")
dev.off()

#region RUnning TFs
PBC_HCC@meta.data$celltypes = gsub("Stem-like CD4 T cells","Stem-like CD4 T cell",PBC_HCC@meta.data$celltypes)
stem_cells = rownames(PBC_HCC@meta.data[grep("Stem-like CD4 T cell",PBC_HCC@meta.data$celltypes),])

PBC_HCC_SM = subset(PBC_HCC, cells = stem_cells)
PBC_HCC_SM = NormalizeData(PBC_HCC_SM)
write.table(PBC_HCC_SM@assays$RNA$data, paste0("/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/Table/PBC_HCC_SM_normalize.txt"), col.names=T, row.names=T, quote=F, sep="\t")

DefaultAssay(PBC_HCC) <- "RNA"
PBC_HCC <- NormalizeData(PBC_HCC)
write.table(PBC_HCC@assays$RNA$data, paste0("/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/Table/PBC_HCC_normalize.txt"), col.names=T, row.names=T, quote=F, sep="\t")
