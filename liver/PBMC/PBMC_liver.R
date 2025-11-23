#### After aggregation Performing QC and combining and integration
### Loading the data and performing the QC on the object
library(Seurat)
samples <- c(
    "pbc_pbmc1_2_aggregated", "pbc_pbmc2_2_aggregated", "pbc_pbmc3_2_aggregated", "pbc_pbmc4_2_aggregated", "pbc_pbmc5_2_aggregated",
    "control_pbmc1_2_aggregated", "control_pbmc2_2_aggregated", "control_pbmc3_2_aggregated", "control_pbmc4_2_aggregated", "control_pbmc5_2_aggregated"
)

source("/diazlab/data3/.abhinav/resources/all_scripts/R/scRNA_QC2.R")
for (i in 1:length(samples)) {
    sample_obj <- scRNA_QC(
        Dir = "/diazlab/data3/.abhinav/.immune/Kyoto/liver/preprocessing/", Sample = samples[i],
        saveDir = "/diazlab/data3/.abhinav/.immune/Kyoto/liver/PBMC/analysis/"
    )
    assign(paste0(samples[i], "_obj"), sample_obj)
}

objname <- ls(pattern = "_aggregated_obj")

## We have to combine the RNA singlet object to perform the normalization for the unwanted the effects.
# combined <- merge(
#     x = get(objname[1]), y = c(
#         get(objname[2]), get(objname[3]), get(objname[4]), get(objname[5]), get(objname[6]),
#         get(objname[7]), get(objname[8]), get(objname[9]), get(objname[10])
#     ),
#     add.cell.ids = objname, project = "combined"
# )

combined <- merge(
    x = get(objname[1]), y = c(
        get(objname[2]), get(objname[3]), get(objname[4]), get(objname[5]), get(objname[6]),
        get(objname[7]), get(objname[8]), get(objname[9]), get(objname[10])
    ),
    add.cell.ids = objname, project = "combined"
)

dir.create("/diazlab/data3/.abhinav/.immune/Kyoto/liver/PBMC/analysis/saveRDS", showWarnings = FALSE)
saveRDS(combined, "/diazlab/data3/.abhinav/.immune/Kyoto/liver/PBMC/analysis/saveRDS/combined.RDS")

### Since now we have moved to the Kyoto Google cloud
library(Seurat)
combined = readRDS("/mnt/data/projects/liver/PBMC/analysis/saveRDS/combined.RDS")


##### Below this code has not been run for PBMC

# source("/diazlab/data3/.abhinav/resources/all_scripts/R/scCITESeq_sctransform_V2.R")
# objname = "CD4"
# Assay = "RNA"
# process = "sctransform"

# combined_sctransformed <- sctransform_V2_integration(obj = combined, saveDir = savedir, ngenes = 4000,
#                                                     regress = c("nCount_RNA"),
#                                                     dims = 30,
#                                                     Assay = Assay, process = process, objname = objname,
#                                                     split_by = "orig.ident",
#                                                     reference = NULL,
#                                                     sample_tree = NULL)

### Performing the Harmony integration. Much faster and easy to run
# combined <- subset(obj, cells = cells_to_keep)  # Subset the Seurat object

library(harmony)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunHarmony(combined, group.by.vars = "orig.ident")

combined <- RunUMAP(
    combined,
    reduction.key = "harmonyUMAP_",
    reduction = "harmony",
    reduction.name = "harmonyumap",
    dims = 1:30
)

savedir <- "/mnt/data/projects/liver/PBMC/"
dir.create(paste0(savedir, "analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/combined_samples_2.pdf"), width = 6, height = 5)
DimPlot(combined, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.4)

dir.create(paste0(savedir, "analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/combined_cluster_2.pdf"), width = 6, height = 5)
DimPlot(combined, reduction = "harmonyumap", group.by = "seurat_clusters", label = TRUE, label.size = 4)
dev.off()

dir.create(paste0(savedir, "analysis/featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/featureplots/Tcell_Bcells_2.pdf"), width = 10)
FeaturePlot(combined, c("CD3E", "CD3D", "CD3G", "CD4", "CD8A", "CD8B", "CD19", "CD20", "CD79A", "CD22", "CD27"), reduction = "harmonyumap")
dev.off()

dir.create(paste0(savedir, "analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Tcell_Bcells_nopoint_2.pdf"), width = 12, height = 5)
VlnPlot(combined, c("CD3E", "CD3D", "CD3G", "CD4", "CD8A", "CD8B", "CD19", "CD20", "CD79A", "CD22", "CD27"), pt.size = 0)
dev.off()

genes <- c(
    "TCF7", "LEF1", "CCR7", "SELL",
    "GZMK", "CCL5", "IL7R", "GZMB", "PRF1", "GNLY", "NKG7",
    "PDCD1", "TIGIT", "LAG3", "HAVCR2", "FOXP3", "CTLA4",
    "SDC1", "PRDM1", "TBX21", "ITGAX", "FCRL5"
)

dir.create(paste0(savedir, "analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Tcellstates_nopoint_2.pdf"), width = 12, height = 10)
VlnPlot(combined, genes, pt.size = 0)
dev.off()

### Running with the join layers
# Idents(combined) <- combined@meta.data$RNA_snn_res.0.4
combined <- JoinLayers(combined)
combined<- NormalizeData(combined)
markers <- FindAllMarkers(combined)
dir.create(paste0(savedir,"analysis/Table"), showWarnings = FALSE)
write.table(markers, paste0(savedir,"analysis/Table/marker_harmony.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

saveRDS(combined, paste0(savedir,"analysis/saveRDS/combined_integrated.RDS"))
# endregion

#region reference
reference = readRDS("/mnt/data/resources/pbmc_multimodal_2023.rds")
# gc()

pdf(paste0(savedir,"analysis/UMAP/reference_umap.pdf"))
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
dev.off()

library(Seurat)
combined = readRDS("/mnt/data/projects/liver/PBMC/analysis/saveRDS/combined_integrated.RDS")
combined <- SCTransform(combined, verbose = TRUE, vst.flavor='v2')

reference = readRDS("/mnt/data/resources/pbmc_multimodal_2023.rds")

anchors <- FindTransferAnchors(
  reference = reference,
  query = combined,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

combined <- MapQuery(
  anchorset = anchors,
  query = combined,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

p1 = DimPlot(combined, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(combined, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()

savedir = "/mnt/data/projects/liver/PBMC/"
pdf(paste0(savedir,"analysis/UMAP/reference_umap.pdf"), width = 12, height = 5)
p1 + p2
dev.off()

p1 = DimPlot(combined, reduction = "harmonyumap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(combined, reduction = "harmonyumap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()

savedir = "/mnt/data/projects/liver/PBMC/"
pdf(paste0(savedir,"analysis/UMAP/PBMC_trasnfer_label_umap_harmony.pdf"), width = 12, height = 5)
p1 + p2
dev.off()

saveRDS(combined, paste0(savedir,"analysis/saveRDS/combined_integrated_mapped.RDS"))

combined@meta.data$status <- gsub("_pbmc.*.","",combined@meta.data$orig.ident)

p1 = DimPlot(combined, reduction = "harmonyumap", group.by = "predicted.celltype.l2", split.by = "status", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(combined, reduction = "ref.umap", group.by = "predicted.celltype.l2", split.by = "status", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()

savedir = "/mnt/data/projects/liver/PBMC/"
pdf(paste0(savedir,"analysis/UMAP/PBMC_splitted_PBC_control_celltypes_label2.pdf"),width = 12, height = 12)
p1 + p2
dev.off()

pdf(paste0(savedir,"analysis/UMAP/PBMC_splitted_PBC_control.pdf"),width = 10, height = 5)
DimPlot(combined, reduction = "harmonyumap", split.by = "status", group.by = "status") + NoLegend()
dev.off()

combined@meta.data$status <- gsub("_.*.","",combined@meta.data$orig.ident)
pdf(paste0(savedir,"UMAP/PBMC_PBC_control.pdf"),width = 5.5, height = 5)
DimPlot(combined, reduction = "harmonyumap", group.by = "status")
dev.off()

combined_status_celltypel2 <- t(table(combined@meta.data$status, combined@meta.data$predicted.celltype.l2))
combined_status_celltypel1 <- t(table(combined@meta.data$status, combined@meta.data$predicted.celltype.l1))
combined_status_clusters <- t(table(combined@meta.data$status, combined@meta.data$seurat_clusters))

write.table(combined_status_celltypel2, paste0(savedir,"analysis/Table/combined_status_celltypel2.txt"), sep = "\t", col.names = T, row.names = T, quote =F)

#region TCF7
# Identify TCF7 gene expression difference between PBC and Control
library(Seurat)
library(dplyr)
library(ggplot2)

# Assuming your Seurat object is called 'pbmc'
# and has metadata columns:
# 1, Healthy vs PBC naive CD4 T cells?
# 2, Healthy vs PBC naive CD8 T cells?
# 3, Naive CD4 vs Naive CD8 in healthy?
# 4, Naive CD4 vs Naive CD8 in PBC?

library(Seurat)
savedir = "/mnt/data/projects/liver/PBMC/analysis/"
combined <- readRDS("/mnt/data/projects/liver/PBMC/analysis/saveRDS/combined_integrated_mapped.RDS")
DefaultAssay(combined) <- "RNA"

CD4_naive_cells = rownames(combined@meta.data[grep("CD4 Naive",combined@meta.data$predicted.celltype.l2),])
Naive_CD4 <- subset(combined, cells = CD4_naive_cells)
naive_CD4_mat = AggregateExpression(Naive_CD4, group.by="orig.ident")
naive_CD4_RNA <- naive_CD4_mat$RNA
naive_CD4_RNA_TCF7 <- as.data.frame(naive_CD4_RNA[grep("^TCF7$",rownames(naive_CD4_RNA)),])
colnames(naive_CD4_RNA_TCF7) <- "TCF7_expr"
columnsum <- as.vector(colSums(naive_CD4_RNA))
naive_CD4_RNA_TCF7$Total_UMI <- as.vector(colSums(naive_CD4_RNA))
naive_CD4_RNA_TCF7$TCF7_logexpr = log1p((naive_CD4_RNA_TCF7$TCF7_expr / columnsum)*1e6)
naive_CD4_RNA_TCF7$condition <- gsub('-pbmc.*.',"",rownames(naive_CD4_RNA_TCF7))
ttesting <- t.test(naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "control","TCF7_logexpr"], naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "pbc","TCF7_logexpr"])
pvalue <- ttesting$p.value
log2FC <- log2(mean(naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "pbc","TCF7_logexpr"])/mean(naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "control","TCF7_logexpr"]))

pdf(paste0(savedir,"analysis/Table/naive_CD4_RNA_TCF7_condition.pdf"))
ggplot(naive_CD4_RNA_TCF7, aes(x = condition, y = TCF7_logexpr, fill = condition)) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  labs(
    title = paste("Naive CD4 p=",pvalue,"FC=",log2FC),
    x = "Condition",
    y = "log(Expression)"
  ) +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("control" = "#1f77b4", "pbc" = "#ff7f0e"))
dev.off()

## saving the file
write.table(naive_CD4_RNA_TCF7, paste0(savedir,"Table/naive_CD4_RNA_TCF7_condition.txt"), sep = "\t", row.names = T, col.names = T, quote = F)

#region CD8
CD4_naive_cells = rownames(combined@meta.data[grep("CD8 Naive",combined@meta.data$predicted.celltype.l2),])
Naive_CD4 <- subset(combined, cells = CD4_naive_cells)
naive_CD4_mat = AggregateExpression(Naive_CD4, group.by="orig.ident")
naive_CD4_RNA <- naive_CD4_mat$RNA
naive_CD4_RNA_TCF7 <- as.data.frame(naive_CD4_RNA[grep("^TCF7$",rownames(naive_CD4_RNA)),])
colnames(naive_CD4_RNA_TCF7) <- "TCF7_expr"
columnsum <- as.vector(colSums(naive_CD4_RNA))
naive_CD4_RNA_TCF7$Total_UMI <- as.vector(colSums(naive_CD4_RNA))
naive_CD4_RNA_TCF7$TCF7_logexpr = log1p((naive_CD4_RNA_TCF7$TCF7_expr / columnsum)*1e6)
naive_CD4_RNA_TCF7$condition <- gsub('-pbmc.*.',"",rownames(naive_CD4_RNA_TCF7))
ttesting <- t.test(naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "control","TCF7_logexpr"], naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "pbc","TCF7_logexpr"])
pvalue <- ttesting$p.value
log2FC <- log2(mean(naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "pbc","TCF7_logexpr"])/mean(naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "control","TCF7_logexpr"]))

pdf(paste0(savedir,"analysis/Table/naive_CD8_RNA_TCF7_condition.pdf"))
ggplot(naive_CD4_RNA_TCF7, aes(x = condition, y = TCF7_logexpr, fill = condition)) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  labs(
    title = paste("Naive CD8 p=",pvalue,"FC=",log2FC),
    x = "Condition",
    y = "log(Expression)"
  ) +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("control" = "#1f77b4", "pbc" = "#ff7f0e"))
dev.off()

## saving the file
write.table(naive_CD4_RNA_TCF7, paste0(savedir,"Table/naive_CD8_RNA_TCF7_condition.txt"), sep = "\t", row.names = T, col.names = T, quote = F)


#region CD4 and CD8
CD4_CD8_naive_cells = rownames(combined@meta.data[grep("CD8 Naive|CD4 Naive",combined@meta.data$predicted.celltype.l2),])
Naive_CD4 <- subset(combined, cells = CD4_CD8_naive_cells)

# Control and PBC
Naive_CD4@meta.data$celltype_samples <- paste(Naive_CD4@meta.data$predicted.celltype.l2, Naive_CD4@meta.data$orig.ident, sep = "_") %>% gsub(" Naive","",.)
control_CD4_CD8 <- rownames(Naive_CD4@meta.data[grep("_control",Naive_CD4@meta.data$celltype_samples),])

Naive_control_CD4_CD8 <- subset(Naive_CD4, cells = control_CD4_CD8)

naive_CD4_mat = AggregateExpression(Naive_control_CD4_CD8, group.by="celltype_samples")
naive_CD4_RNA <- naive_CD4_mat$RNA
naive_CD4_RNA_TCF7 <- as.data.frame(naive_CD4_RNA[grep("^TCF7$",rownames(naive_CD4_RNA)),])
colnames(naive_CD4_RNA_TCF7) <- "TCF7_expr"
columnsum <- as.vector(colSums(naive_CD4_RNA))
naive_CD4_RNA_TCF7$Total_UMI <- as.vector(colSums(naive_CD4_RNA))
naive_CD4_RNA_TCF7$TCF7_logexpr = log1p((naive_CD4_RNA_TCF7$TCF7_expr / columnsum)*1e6)
naive_CD4_RNA_TCF7$condition <- gsub('-control-pbmc.*.',"",rownames(naive_CD4_RNA_TCF7))
ttesting <- t.test(naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "CD4","TCF7_logexpr"], naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "CD8","TCF7_logexpr"])
pvalue <- ttesting$p.value
log2FC <- log2(mean(naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "CD4","TCF7_logexpr"])/mean(naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "CD8","TCF7_logexpr"]))

pdf(paste0(savedir,"analysis/Table/naive_CD4_vs_CD8_RNA_TCF7_condition.pdf"))
ggplot(naive_CD4_RNA_TCF7, aes(x = condition, y = TCF7_logexpr, fill = condition)) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  labs(
    title = paste("Naive CD4 vs CD8 p=",round(pvalue,3),"FC=",round(log2FC,3)),
    x = "Condition",
    y = "log(Expression)"
  ) +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("CD4" = "#1f77b4", "CD8" = "#ff7f0e"))
dev.off()

## saving the file
write.table(naive_CD4_RNA_TCF7, paste0(savedir,"Table/naive_Control_RNA_TCF7_CD4_CD8.txt"), sep = "\t", row.names = T, col.names = T, quote = F)

#region PBC TCF7
CD4_naive_cells = rownames(combined@meta.data[grep("CD8 Naive|CD4 Naive",combined@meta.data$predicted.celltype.l2),])
Naive_CD4 = subset(combined, cells = CD4_naive_cells)
Naive_CD4@meta.data$celltype_samples <- paste(Naive_CD4@meta.data$predicted.celltype.l2, Naive_CD4@meta.data$orig.ident, sep = "_") %>% gsub(" Naive","",.)
control_CD4_CD8 <- rownames(Naive_CD4@meta.data[grep("_pbc",Naive_CD4@meta.data$celltype_samples),])
Naive_control_CD4_CD8 <- subset(Naive_CD4, cells = control_CD4_CD8)

naive_CD4_mat = AggregateExpression(Naive_control_CD4_CD8, group.by="celltype_samples")
naive_CD4_RNA <- naive_CD4_mat$RNA
naive_CD4_RNA_TCF7 <- as.data.frame(naive_CD4_RNA[grep("^TCF7$",rownames(naive_CD4_RNA)),])
colnames(naive_CD4_RNA_TCF7) <- "TCF7_expr"
columnsum <- as.vector(colSums(naive_CD4_RNA))
naive_CD4_RNA_TCF7$Total_UMI <- as.vector(colSums(naive_CD4_RNA))
naive_CD4_RNA_TCF7$TCF7_logexpr = log1p((naive_CD4_RNA_TCF7$TCF7_expr / columnsum)*1e6)
naive_CD4_RNA_TCF7$condition <- gsub('-pbc-pbmc.*.',"",rownames(naive_CD4_RNA_TCF7))
ttesting <- t.test(naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "CD4","TCF7_logexpr"], naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "CD8","TCF7_logexpr"])
pvalue <- ttesting$p.value
log2FC <- log2(mean(naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "CD4","TCF7_logexpr"])/mean(naive_CD4_RNA_TCF7[naive_CD4_RNA_TCF7$condition == "CD8","TCF7_logexpr"]))

pdf(paste0(savedir,"analysis/Table/naive_PBC_CD4_vs_CD8_RNA_TCF7_condition.pdf"))
ggplot(naive_CD4_RNA_TCF7, aes(x = condition, y = TCF7_logexpr, fill = condition)) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  labs(
    title = paste("Naive PBC CD4 vs CD8 p=",round(pvalue,3),"FC=",round(log2FC,3)),
    x = "Condition",
    y = "log(Expression)"
  ) +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("CD4" = "#1f77b4", "CD8" = "#ff7f0e"))
dev.off()

write.table(naive_CD4_RNA_TCF7, paste0(savedir,"Table/naive_PBC_RNA_TCF7_CD4_CD8.txt"), sep = "\t", row.names = T, col.names = T, quote = F)

#endregion

#region RAM score
PBMC <- readRDS("/mnt/data/projects/Kyoto/liver/PBMC/analysis/saveRDS/combined_integrated_mapped.RDS")

savedir = "/mnt/data/projects/Kyoto/liver/PBMC/analysis/naive/"
dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir,"UMAP/PBMC_celltypes.pdf"), width = 10, height = 5.5)
DimPlot(PBMC, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = T)
dev.off()

naive_cellnames <- rownames(PBMC@meta.data[grep("CD4 Naive", PBMC@meta.data$predicted.celltype.l2), ])
PBMC_naive <- subset(PBMC, cells = naive_cellnames)

pdf(paste0(savedir,"UMAP/PBMC_subset_celltypes.pdf"))
DimPlot(PBMC_naive, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = T)
dev.off()

### Adding the gene list
oldvsyoung= read.table("/mnt/data/projects/resource/old_vs_young_Claires_Nature_each_celltype.txt", header = TRUE, sep = "\t")
oldvsyoung_naive <- oldvsyoung[grep("Core naive CD4 T cell",oldvsyoung$AIFI_L3),]
naive_aging_high <- oldvsyoung_naive[oldvsyoung_naive$padj < 0.05 & oldvsyoung_naive$log2fc < 0,] ## taking only high older adult for aging
rm(genelist)
genelist <- list()
genelist[[1]] <- naive_aging_high$gene

PBMC_naive[["RNA"]] <- split(PBMC_naive[["RNA"]], f = PBMC_naive$orig.ident)
PBMC_naive <- NormalizeData(PBMC_naive)
PBMC_naive <- AddModuleScore(PBMC_naive, genelist)
# PBC_HCC_naive_req@meta.data <- PBC_HCC_naive_req@meta.data[,grep("Cluster",colnames(PBC_HCC_naive_req@meta.data),invert = TRUE)]

req_index = grep("Cluster",colnames(PBMC_naive@meta.data))
colnames(PBMC_naive@meta.data)[req_index] = "RAM_coreCD4_naive_score_lfc_0.25_2"

# PBMC_naive@meta.data$condition <- gsub("_.*.","",PBMC_naive@meta.data$orig.ident)

library(ggplot2)
dir.create(paste0(savedir,"vlnplot"), showWarnings = FALSE)
pdf(paste0(savedir,"vlnplot/RAMcoreCD4naive_aging_boxplot_lfc_0.25_2.pdf"))
VlnPlot(PBMC_naive, "RAM_coreCD4_naive_score", group.by = "condition", pt.size =0) + geom_boxplot()
dev.off()

pdf(paste0(savedir,"vlnplot/RAMcoreCD4naive_aging_score_point_lfc_0.25_2.pdf"))
VlnPlot(PBMC_naive, "RAM_coreCD4_naive_score", group.by = "condition") + geom_boxplot()
dev.off()

### For young
naive_aging_low <- oldvsyoung_naive[oldvsyoung_naive$padj < 0.05 & oldvsyoung_naive$log2fc < -0.5,] ## taking only high older adult for aging
rm(genelist)
genelist <- list()
genelist[[1]] <- naive_aging_low$gene

# PBMC_naive <- NormalizeData(PBMC_naive)
PBMC_naive <- AddModuleScore(PBMC_naive, genelist)
# PBC_HCC_naive_req@meta.data <- PBC_HCC_naive_req@meta.data[,grep("Cluster",colnames(PBC_HCC_naive_req@meta.data),invert = TRUE)]

req_index = grep("Cluster",colnames(PBMC_naive@meta.data))
colnames(PBMC_naive@meta.data)[req_index] = "RAM_coreCD4_naive_young_score"

# PBMC_naive@meta.data$condition <- gsub("_.*.","",PBMC_naive@meta.data$orig.ident)

library(ggplot2)
dir.create(paste0(savedir,"vlnplot"), showWarnings = FALSE)
pdf(paste0(savedir,"vlnplot/RAM_coreCD4_naive_young_score_boxplot.pdf"))
VlnPlot(PBMC_naive, "RAM_coreCD4_naive_young_score", group.by = "condition", pt.size =0) + geom_boxplot()
dev.off()

pdf(paste0(savedir,"vlnplot/RAM_coreCD4_naive_young_score_point.pdf"))
VlnPlot(PBMC_naive, "RAM_coreCD4_naive_young_score", group.by = "condition") + geom_boxplot()
dev.off()


# require_cell = rownames(PBC_HCC_naive_req@meta.data[grep("control",PBC_HCC_naive_req@meta.data$com_condition,invert=TRUE),])
# PBC_HCC_req2 = subset(PBC_HCC_naive_req, cells = require_cell)

anova_res <- aov(Tn_score ~ treatment_Resp, data = PBC_HCC_naive_req@meta.data)
summary(anova_res)

dir.create(paste0(savedir,"saveRDS"), showWarnings = FALSE)
saveRDS(PBC_HCC_naive_req, paste0(savedir,"saveRDS/PBC_HCC_Tfh.RDS"))

#endregion

