library(Seurat)
library(harmony)
PBC_obj <- readRDS("/mnt/data/projects/liver/Tcells/remove_more_mps/CD4/Tregs/saveRDS/Treg2.RDS")
HCC_obj <- readRDS("/mnt/data/projects/NatMed_HCC/analysis/Tregs/removed_clus12_14/saveRDS/Treg_removed_7_8.RDS")

HCC_obj@meta.data$orig.ident <- paste0("Sample_", HCC_obj@meta.data$cell_to_sample_ID)
HCC_obj@assays$RNA <- HCC_obj@assays$originalexp
DefaultAssay(HCC_obj) <- "RNA"
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

saveRDS(PBC_HCC, paste0(savedir,"saveRDS/PBC_HCC_Tregs.RDS"))

# PBC_HCC <- readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/saveRDS/PBC_HCC.RDS")
PBC_HCC@meta.data$com_condition <- PBC_HCC@meta.data$condition

genelist = list()
genelist[[1]] = c("IL2RA", "CTLA4", "ICOS", "TNFRSF4", "TNFRSF18", "CCR8", "BATF",
'PRDM1', "IRF4", "ENTPD1", "LAYN", 'TIGIT', 'HAVCR2')

PBC_HCC@meta.data[is.na(PBC_HCC@meta.data$com_condition),"com_condition"] <- PBC_HCC@meta.data[is.na(PBC_HCC@meta.data$com_condition),"treatment_Resp"]

cellnames = rownames(PBC_HCC@meta.data[grep("antiPD1_NR|antiPD1_R|control|pbc",PBC_HCC@meta.data$com_condition),])

PBC_HCC_req = subset(PBC_HCC, cells = cellnames)
PBC_HCC_req[["RNA"]] <- split(PBC_HCC_req[["RNA"]], f = PBC_HCC_req$com_condition)
PBC_HCC_req <- NormalizeData(PBC_HCC_req)
PBC_HCC_req <- AddModuleScore(PBC_HCC_req, genelist)

req_index = grep("Cluster",colnames(PBC_HCC_req@meta.data))
colnames(PBC_HCC_req@meta.data)[req_index] = "Tregs_score"

dir.create(paste0(savedir,"vlnplot"), showWarnings = FALSE)
pdf(paste0(savedir,"vlnplot/Treg_score_boxplot.pdf"))
VlnPlot(PBC_HCC_req, "Tregs_score", group.by = "com_condition", pt.size =0) + geom_boxplot()
dev.off()

pdf(paste0(savedir,"vlnplot/Treg_score_point.pdf"))
VlnPlot(PBC_HCC_req, "Tregs_score", group.by = "com_condition") + geom_boxplot()
dev.off()

require_cell = rownames(PBC_HCC_req@meta.data[grep("control",PBC_HCC_req@meta.data$com_condition,invert=TRUE),])
PBC_HCC_req2 = subset(PBC_HCC_req, cells = require_cell)

t.test(Tregs_score ~ com_condition, data = PBC_HCC_req2@meta.data[PBC_HCC_req2@meta.data$com_condition %in% c("antiPD1_NR", "antiPD1_R"), ])
t.test(Tregs_score ~ com_condition, data = PBC_HCC_req2@meta.data[PBC_HCC_req2@meta.data$com_condition %in% c("antiPD1_NR", "pbc"), ])
t.test(Tregs_score ~ com_condition, data = PBC_HCC_req2@meta.data[PBC_HCC_req2@meta.data$com_condition %in% c("antiPD1_R", "pbc"), ])

wilcox.test(Tregs_score ~ com_condition, data = PBC_HCC_req2@meta.data[PBC_HCC_req2@meta.data$com_condition %in% c("antiPD1_NR", "antiPD1_R"), ])
wilcox.test(Tregs_score ~ com_condition, data = PBC_HCC_req2@meta.data[PBC_HCC_req2@meta.data$com_condition %in% c("antiPD1_NR", "pbc"), ])
wilcox.test(Tregs_score ~ com_condition, data = PBC_HCC_req2@meta.data[PBC_HCC_req2@meta.data$com_condition %in% c("antiPD1_R", "pbc"), ])

pdf(paste0(savedir,"vlnplot/Treg_score_point_2.pdf"))
VlnPlot(PBC_HCC_req2, "Tregs_score", group.by = "com_condition") + geom_boxplot()
dev.off()

pdf(paste0(savedir,"vlnplot/Treg_score_boxplot_2.pdf"))
VlnPlot(PBC_HCC_req2, "Tregs_score", group.by = "com_condition", pt.size =0) + geom_boxplot()
dev.off()

pdf(paste0(savedir,"vlnplot/Treg_score_boxplot.pdf"))
VlnPlot(PBC_HCC_req2, "Tregs_score", group.by = "com_condition", pt.size =0)
dev.off()

saveRDS(PBC_HCC_req2, paste0(savedir,"saveRDS/PBC_HCC_Tregs_nocontrol.RDS"))
