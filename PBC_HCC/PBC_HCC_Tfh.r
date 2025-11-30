library(Seurat)
library(harmony)
PBC_HCC = readRDS("/mnt/data/projects/Kyoto/PBC_HCC/saveRDS/PBC_HCC_CCA.RDS")

dir.create("/mnt/data/projects/Kyoto/PBC_HCC/Tfh/UMAP/", showWarnings = FALSE, recursive  = T)

savedir = "/mnt/data/projects/Kyoto/PBC_HCC/Tfh/"

pdf(paste0(savedir,"UMAP/PBC_HCC_celltypes.pdf"))
DimPlot(PBC_HCC, reduction = "umap", group.by = "celltypes", label = T)
dev.off()

pdf(paste0(savedir,"UMAP/PBC_HCC_seurat_clusters.pdf"))
DimPlot(PBC_HCC, reduction = "umap", group.by = "seurat_clusters", label = T)
dev.off()

pdf(paste0(savedir,"UMAP/PBC_HCC_res0.4.pdf"))
DimPlot(PBC_HCC, reduction = "umap", group.by = "res0.4", label = T)
dev.off()

Tph_cellnames <- rownames(PBC_HCC@meta.data[grep("Tph", PBC_HCC@meta.data$celltypes), ])
PBC_HCC_Tph <- subset(PBC_HCC, cells = Tph_cellnames)

### Adding Responder and non-responder
metadata = read.csv("/mnt/data/projects/Kyoto/NatMed_HCC/analysis/GSE206325_sample_annots_Liver_Treated_patients.csv", header = T)
patterns <- metadata$sample_ID
replacements <- metadata$treatment_Resp

library(stringr)
PBC_HCC_Tph@meta.data$treatment_Resp <- PBC_HCC_Tph@meta.data$cell_to_sample_ID
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    PBC_HCC_Tph@meta.data$treatment_Resp <- str_replace_all(PBC_HCC_Tph@meta.data$treatment_Resp, pattern, replacements[i])
    PBC_HCC_Tph@meta.data$treatment_Resp <- gsub("-", "", PBC_HCC_Tph@meta.data$treatment_Resp)
}

PBC_HCC_Tph@meta.data$treatment_Resp[grep("liver",PBC_HCC_Tph@meta.data$orig.ident)] <- gsub("_.*.","",grep("liver",PBC_HCC_Tph@meta.data$orig.ident,value=TRUE))
cellnames <- rownames(PBC_HCC_Tph@meta.data[grep("^[0-9]",PBC_HCC_Tph@meta.data$treatment_Resp, invert = TRUE),])

cellnames <- rownames(PBC_HCC_Tph@meta.data[grep("^[0-9]",PBC_HCC_Tph@meta.data$treatment_Resp, invert = TRUE),])
PBC_HCC_Tph_req = subset(PBC_HCC_Tph, cells = cellnames)

genelist <- list()
genelist[[1]] <- read.table("/mnt/data/projects/resource/CD4_Tfh", header = FALSE)[,1]
PBC_HCC_Tph_req[["RNA"]] <- split(PBC_HCC_Tph_req[["RNA"]], f = PBC_HCC_Tph_req$treatment_Resp)
PBC_HCC_Tph_req <- NormalizeData(PBC_HCC_Tph_req)
PBC_HCC_Tph_req <- AddModuleScore(PBC_HCC_Tph_req, genelist)
# PBC_HCC_Tph_req@meta.data <- PBC_HCC_Tph_req@meta.data[,grep("Cluster",colnames(PBC_HCC_Tph_req@meta.data),invert = TRUE)]

req_index = grep("Cluster",colnames(PBC_HCC_Tph_req@meta.data))
colnames(PBC_HCC_Tph_req@meta.data)[req_index] = "Tfh_score"

library(ggplot2)
dir.create(paste0(savedir,"vlnplot"), showWarnings = FALSE)
pdf(paste0(savedir,"vlnplot/Tfh_score_boxplot.pdf"))
VlnPlot(PBC_HCC_Tph_req, "Tfh_score", group.by = "treatment_Resp", pt.size =0) + geom_boxplot()
dev.off()

pdf(paste0(savedir,"vlnplot/Tfh_score_point.pdf"))
VlnPlot(PBC_HCC_Tph_req, "Tfh_score", group.by = "treatment_Resp") + geom_boxplot()
dev.off()

# require_cell = rownames(PBC_HCC_Tph_req@meta.data[grep("control",PBC_HCC_Tph_req@meta.data$com_condition,invert=TRUE),])
# PBC_HCC_req2 = subset(PBC_HCC_Tph_req, cells = require_cell)

anova_res <- aov(Tfh_score ~ treatment_Resp, data = PBC_HCC_Tph_req@meta.data)
summary(anova_res)
#                   Df Sum Sq Mean Sq F value Pr(>F)
# treatment_Resp     3   29.3   9.754   179.3 <2e-16 ***
# Residuals      16772  912.3   0.054
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(anova_res)
#   Tukey multiple comparisons of means
#     95% family-wise confidence level
# Fit: aov(formula = Tfh_score ~ treatment_Resp, data = PBC_HCC_Tph_req@meta.data)
# $treatment_Resp
#                             diff          lwr         upr     p adj
# antiPD1_R-antiPD1_NR  0.05421619  0.042198691  0.06623370 0.0000000
# control-antiPD1_NR   -0.06049371 -0.084783459 -0.03620396 0.0000000
# pbc-antiPD1_NR       -0.04446840 -0.055919618 -0.03301718 0.0000000
# control-antiPD1_R    -0.11470990 -0.138982608 -0.09043720 0.0000000
# pbc-antiPD1_R        -0.09868459 -0.110099612 -0.08726958 0.0000000
# pbc-control           0.01602531 -0.007972072  0.04002269 0.3153945

dir.create(paste0(savedir,"saveRDS"), showWarnings = FALSE)
saveRDS(PBC_HCC_Tph_req, paste0(savedir,"saveRDS/PBC_HCC_Tfh.RDS"))

savedir = "/mnt/data/projects/Kyoto/PBC_HCC/Tfh/"
PBC_HCC_Tph_req = readRDS(paste0(savedir,"saveRDS/PBC_HCC_Tfh.RDS"))

treatment_levels = c("control","pbc","antiPD1_NR","antiPD1_R")
PBC_HCC_Tph_req@meta.data$treatment_Resp <- factor(PBC_HCC_Tph_req@meta.data$treatment_Resp, levels = treatment_levels)

pdf(paste0(savedir,"vlnplot/Tfh_score_point_2.pdf"))
VlnPlot(PBC_HCC_Tph_req, "Tfh_score", group.by = "treatment_Resp") + geom_boxplot()
dev.off()

pdf(paste0(savedir,"vlnplot/Tfh_score_boxplot_2.pdf"))
VlnPlot(PBC_HCC_Tph_req, "Tfh_score", group.by = "treatment_Resp", pt.size =0) + geom_boxplot()
dev.off()
å