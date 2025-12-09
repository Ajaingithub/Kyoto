#region Combined PBC_HCC for 7 clus
library(Seurat)
PBC_HCC = readRDS("/mnt/data/projects/Kyoto/PBC_HCC/saveRDS/PBC_HCC_CCA_required.RDS")
dir.create(savedir, showWarnings =  FALSE)
setwd(savedir)

dir.create("UMAP", showWarnings =  FALSE)  
pdf(paste0(savedir,"UMAP/PBC_HCC.pdf"))
DimPlot(PBC_HCC, reduction = "umap", group.by = "res0.4", label = T, label.size = 5)
dev.off()

# cyto_HSP with cluster 5
savedir = "/mnt/data/projects/Kyoto/PBC_HCC/cyto_HSP/cyto_HSP_new_cluster/"
dir.create(savedir, showWarnings = FALSE, recursive = TRUE)
cyto_HSP_cellnames = rownames(PBC_HCC@meta.data[grep("^3$",PBC_HCC@meta.data$res0.4),])
cyto_HSP = subset(PBC_HCC, cells = cyto_HSP_cellnames)

cyto_HSP <- NormalizeData(cyto_HSP)

rm(genelist)
genelist = list()
genelist[[1]] = c("GZMA", "GZMB", "GZMK", "PRF1", "NKG7")
genelist[[2]] = c("HSPA1A", "HSPA1B", "HSPH1", "NR4A1", "HSPA6", "BAG3", "HSPE1", "DNAJB4", "DNAJA1")

cyto_HSP <- AddModuleScore(cyto_HSP, genelist)

req_index = grep("Cluster",colnames(cyto_HSP@meta.data))
colnames(cyto_HSP@meta.data)[req_index] = c("cyto_score_clus3","HSPscore_clus3")

cyto_HSP@meta.data$treatment_Resp = factor(cyto_HSP@meta.data$treatment_Resp, levels = c("control","pbc","antiPD1_NR","antiPD1_R"))

library(ggplot2)
dir.create(paste0(savedir,"vlnplot"), showWarnings = FALSE, recursive = TRUE)
genescore = c("cyto_score_clus3","HSPscore_clus3")
for(i in 1:length(genescore)){
    pdf(paste0(savedir,"vlnplot/",genescore[i],"_boxplot_clus3.pdf"))
    print(VlnPlot(cyto_HSP, genescore[i], group.by = "treatment_Resp", pt.size =0) + geom_boxplot())
    dev.off()
    
    pdf(paste0(savedir,"vlnplot/",genescore[i],"_score_point_clus3.pdf"))
    print(VlnPlot(cyto_HSP, genescore[i], group.by = "treatment_Resp") + geom_boxplot())
    dev.off()
}

dir.create(paste0(savedir,"saveRDS_obj"), showWarnings = FALSE)
saveRDS(cyto_HSP, paste0(savedir,"saveRDS_obj/cyto_HSP_new_clus3.RDS"))


