library(Seurat)
CD4_liver = readRDS("/mnt/data/projects/liver/Tcells/remove_more_mps/CD4/saveRDS/CD4_obj.RDS")
PBMC_liver = readRDS("/mnt/data/projects/liver/PBMC/analysis/saveRDS/combined_integrated_mapped.RDS")

PBMC_liver = DietSeurat(PBMC_liver,layers = "counts",assays = "RNA",dimreducs = "harmonyumap")
CD4_liver = DietSeurat(CD4_liver,layers = "counts",assays = "RNA",dimreducs = "harmonyumap")

### Subsetting CD4 Naive
CD4_naive_cells = rownames(PBMC_liver@meta.data[grep("CD4 Naive",PBMC_liver@meta.data$predicted.celltype.l2),])
PBMC_liver_CD4 <- subset(PBMC_liver, cells = CD4_naive_cells)
liver_CD4_naive = rownames(CD4_liver@meta.data[grep("C3",CD4_liver@meta.data$seurat_clusters_new),])
naive_CD4_liver <- subset(CD4_liver, cells = liver_CD4_naive)

### Common genes
common_genes <- intersect(rownames(naive_CD4_liver), rownames(PBMC_liver_CD4))
naive_CD4_liver <- subset(naive_CD4_liver, features = common_genes) 
PBMC_liver_CD4 <- subset(PBMC_liver_CD4, features = common_genes)

naive_CD4_liver$sample_origin <- "naive_CD4_liver"
PBMC_liver_CD4$sample_origin <- "naive_CD4_PBMC"

combined_CD4 <- merge(naive_CD4_liver, y = PBMC_liver_CD4)

combined_CD4 <- JoinLayers(combined_CD4)
combined_CD4 <- NormalizeData(combined_CD4)

savedir = "/mnt/data/projects/liver/PBMC_vs_Tissue/"
setwd(savedir)
dir.create(paste0(savedir,"Table"), showWarnings = FALSE)

Idents(combined_CD4) <- combined_CD4@meta.data$sample_origin
naive_CD4_liver_vs_PBMC = FindMarkers(combined_CD4, ident.1 = "naive_CD4_liver", ident.2 = "naive_CD4_PBMC")
write.table(naive_CD4_liver_vs_PBMC, paste0(savedir,"Table/naive_CD4_liver_vs_PBMC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

genes = c(
    "TCF7","LEF1","CCR7","TGFB1","SMAD3","ITGB1","ITGB2","ITGB1BP1",
    "ETS1","KLF2","ID3","IRF1","BATF","MAF","IKZF1","FOXO1","CXCR4"
)

dir.create(paste0(savedir,"vlnplots/"), showWarnings = FALSE)
pdf(paste0(savedir,"vlnplots/CD4_genes.pdf"), width = 12, height = 15)
VlnPlot(combined_CD4,genes, ncol = 5)
dev.off()

pdf(paste0(savedir,"vlnplots/CD4_genes_no_points.pdf"), width = 12, height = 15)
VlnPlot(combined_CD4,genes, ncol = 5, pt.size = 0)
dev.off()

dir.create(paste0(savedir,"saveRDS"))
saveRDS(combined_CD4,paste0(savedir,"saveRDS/combined_CD4.RDS"))

### Performing Healthy and Normal samples separately
library(Seurat)
savedir = "/mnt/data/projects/liver/PBMC_vs_Tissue/"
combined_CD4 = readRDS(paste0(savedir,"saveRDS/combined_CD4.RDS"))

# control
control_CD4=rownames(combined_CD4@meta.data[grep("control",combined_CD4@meta.data$orig.ident),])
control_CD4_obj = subset(combined_CD4, cells = control_CD4)

Idents(control_CD4_obj) <- control_CD4_obj@meta.data$sample_origin
naive_CD4_liver_vs_PBMC = FindMarkers(control_CD4_obj, ident.1 = "naive_CD4_liver", ident.2 = "naive_CD4_PBMC")
write.table(naive_CD4_liver_vs_PBMC, paste0(savedir,"Table/naive_control_CD4_liver_vs_PBMC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

genes = c(
    "TCF7","LEF1","CCR7","TGFB1","SMAD3","ITGB1","ITGB2","ITGB1BP1",
    "ETS1","KLF2","ID3","IRF1","BATF","MAF","IKZF1","FOXO1","CXCR4"
)

# dir.create(paste0(savedir,"vlnplots/"))
pdf(paste0(savedir,"vlnplots/CD4_control_genes.pdf"), width = 12, height = 15)
VlnPlot(control_CD4_obj,genes, ncol = 5)
dev.off()

pdf(paste0(savedir,"vlnplots/CD4_control_genes_no_points.pdf"), width = 12, height = 15)
VlnPlot(control_CD4_obj,genes, ncol = 5, pt.size = 0)
dev.off()

saveRDS(control_CD4_obj,paste0(savedir,"saveRDS/control_CD4.RDS"))

### Performing Healthy and PBC samples separately
# control
PBC_CD4=rownames(combined_CD4@meta.data[grep("^pbc",combined_CD4@meta.data$orig.ident),])
PBC_CD4_obj = subset(combined_CD4, cells = PBC_CD4)

Idents(PBC_CD4_obj) <- PBC_CD4_obj@meta.data$sample_origin
naive_CD4_liver_vs_PBMC = FindMarkers(PBC_CD4_obj, ident.1 = "naive_CD4_liver", ident.2 = "naive_CD4_PBMC")
write.table(naive_CD4_liver_vs_PBMC, paste0(savedir,"Table/naive_PBC_CD4_liver_vs_PBMC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

genes = c(
    "TCF7","LEF1","CCR7","TGFB1","SMAD3","ITGB1","ITGB2","ITGB1BP1",
    "ETS1","KLF2","ID3","IRF1","BATF","MAF","IKZF1","FOXO1","CXCR4"
)

# dir.create(paste0(savedir,"vlnplots/"))
pdf(paste0(savedir,"vlnplots/CD4_PBC_genes.pdf"), width = 12, height = 15)
VlnPlot(PBC_CD4_obj,genes, ncol = 5)
dev.off()

pdf(paste0(savedir,"vlnplots/CD4_PBC_genes_no_points.pdf"), width = 12, height = 15)
VlnPlot(PBC_CD4_obj,genes, ncol = 5, pt.size = 0)
dev.off()

saveRDS(PBC_CD4_obj,paste0(savedir,"saveRDS/PBC_CD4.RDS"))


#region CD8
library(Seurat)
CD8_liver = readRDS("/mnt/data/projects/liver/Tcells/remove_more_mps/CD8/analysis/saveRDS/PBC_CD8_Tcells.RDS")
PBMC_liver = readRDS("/mnt/data/projects/liver/PBMC/analysis/saveRDS/combined_integrated_mapped.RDS")

PBMC_liver = DietSeurat(PBMC_liver,layers = "counts",assays = "RNA",dimreducs = "harmonyumap")
CD8_liver = DietSeurat(CD8_liver,layers = "counts",assays = "RNA",dimreducs = "harmonyumap")

### Subsetting CD8 Naive
CD8_naive_cells = rownames(PBMC_liver@meta.data[grep("CD8 Naive",PBMC_liver@meta.data$predicted.celltype.l2),])
PBMC_liver_CD8 <- subset(PBMC_liver, cells = CD8_naive_cells)
liver_CD8_naive = rownames(CD8_liver@meta.data[grep("^5$",CD8_liver@meta.data$seurat_clusters),])
naive_CD8_liver <- subset(CD8_liver, cells = liver_CD8_naive)

### Common genes
common_genes <- intersect(rownames(naive_CD8_liver), rownames(PBMC_liver_CD8))
naive_CD8_liver <- subset(naive_CD8_liver, features = common_genes) 
PBMC_liver_CD8 <- subset(PBMC_liver_CD8, features = common_genes)

naive_CD8_liver$sample_origin <- "naive_CD8_liver"
PBMC_liver_CD8$sample_origin <- "naive_CD8_PBMC"

combined <- merge(naive_CD8_liver, y = PBMC_liver_CD8)
combined <- JoinLayers(combined)
combined <- NormalizeData(combined)

savedir = "/mnt/data/projects/liver/PBMC_vs_Tissue/"
setwd(savedir)
dir.create(paste0(savedir,"Table"), showWarnings = FALSE)

Idents(combined) <- combined@meta.data$sample_origin
naive_CD8_liver_vs_PBMC = FindMarkers(combined, ident.1 = "naive_CD8_liver", ident.2 = "naive_CD8_PBMC")
write.table(naive_CD8_liver_vs_PBMC, paste0(savedir,"Table/naive_CD8_liver_vs_PBMC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

genes = c(
    "TCF7","LEF1","CCR7","TGFB1","SMAD3","ITGB1","ITGB2","ITGB1BP1",
    "ETS1","KLF2","ID3","IRF1","BATF","MAF","IKZF1","FOXO1","CXCR4"
)

# dir.create(paste0(savedir,"vlnplots/"))
pdf(paste0(savedir,"vlnplots/CD8_genes.pdf"), width = 12, height = 15)
VlnPlot(combined,genes, ncol = 5)
dev.off()

pdf(paste0(savedir,"vlnplots/CD8_genes_no_points.pdf"), width = 12, height = 15)
VlnPlot(combined,genes, ncol = 5, pt.size = 0)
dev.off()

saveRDS(combined,paste0(savedir,"saveRDS/combined_CD8.RDS"))

### Performing Healthy and Normal samples separately
library(Seurat)
savedir = "/mnt/data/projects/liver/PBMC_vs_Tissue/"
combined = readRDS(paste0(savedir,"saveRDS/combined_CD8.RDS"))

# control
control_CD8=rownames(combined@meta.data[grep("control",combined@meta.data$orig.ident),])
control_CD8_obj = subset(combined, cells = control_CD8)

Idents(control_CD8_obj) <- control_CD8_obj@meta.data$sample_origin
naive_CD8_liver_vs_PBMC = FindMarkers(control_CD8_obj, ident.1 = "naive_CD8_liver", ident.2 = "naive_CD8_PBMC")
write.table(naive_CD8_liver_vs_PBMC, paste0(savedir,"Table/naive_control_CD8_liver_vs_PBMC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

genes = c(
    "TCF7","LEF1","CCR7","TGFB1","SMAD3","ITGB1","ITGB2","ITGB1BP1",
    "ETS1","KLF2","ID3","IRF1","BATF","MAF","IKZF1","FOXO1","CXCR4"
)

# dir.create(paste0(savedir,"vlnplots/"))
pdf(paste0(savedir,"vlnplots/CD8_control_genes.pdf"), width = 12, height = 15)
VlnPlot(control_CD8_obj,genes, ncol = 5)
dev.off()

pdf(paste0(savedir,"vlnplots/CD8_control_genes_no_points.pdf"), width = 12, height = 15)
VlnPlot(control_CD8_obj,genes, ncol = 5, pt.size = 0)
dev.off()

saveRDS(control_CD8_obj,paste0(savedir,"saveRDS/control_CD8.RDS"))

### Performing Healthy and PBC samples separately
# control
PBC_CD8=rownames(combined@meta.data[grep("^pbc",combined@meta.data$orig.ident),])
PBC_CD8_obj = subset(combined, cells = PBC_CD8)

Idents(PBC_CD8_obj) <- PBC_CD8_obj@meta.data$sample_origin
naive_CD8_liver_vs_PBMC = FindMarkers(PBC_CD8_obj, ident.1 = "naive_CD8_liver", ident.2 = "naive_CD8_PBMC")
write.table(naive_CD8_liver_vs_PBMC, paste0(savedir,"Table/naive_PBC_CD8_liver_vs_PBMC.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

genes = c(
    "TCF7","LEF1","CCR7","TGFB1","SMAD3","ITGB1","ITGB2","ITGB1BP1",
    "ETS1","KLF2","ID3","IRF1","BATF","MAF","IKZF1","FOXO1","CXCR4"
)

# dir.create(paste0(savedir,"vlnplots/"))
pdf(paste0(savedir,"vlnplots/CD8_PBC_genes.pdf"), width = 12, height = 15)
VlnPlot(PBC_CD8_obj,genes, ncol = 5)
dev.off()

pdf(paste0(savedir,"vlnplots/CD8_PBC_genes_no_points.pdf"), width = 12, height = 15)
VlnPlot(PBC_CD8_obj,genes, ncol = 5, pt.size = 0)
dev.off()

saveRDS(PBC_CD8_obj,paste0(savedir,"saveRDS/PBC_CD8.RDS"))


