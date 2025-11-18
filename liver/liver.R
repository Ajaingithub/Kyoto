#region All Cells
### We have performed BAM to Fastq
### Data downloaded from here https://ngdc.cncb.ac.cn/gsa-human/browse/HRA008003

# #!/bin/bash
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=1
# #SBATCH --mem=100G
# #SBATCH --job-name=fastq
# #SBATCH --time=1-00:00:00
# #SBATCH --output=%x-%j.out
# #SBATCH --error=%x-%j.err
# #SBATCH --mail-user=abhinav.jain@ucsf.edu
# #SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
# #SBATCH --array=1-18

# #### List the BAM files
# BAM_FILES=(
#     "HRR1849445.bam" "HRR1849447.bam"
#     "HRR1849451.bam" "HRR1849453.bam"
#     "HRR1849455.bam" "HRR1849457.bam"
#     "HRR1849459.bam" "HRR1849461.bam"
#     "HRR1849463.bam" "HRR1849446.bam"
#     "HRR1849448.bam" "HRR1849450.bam"
#     "HRR1849452.bam" "HRR1849454.bam"
#     "HRR1849456.bam" "HRR1849458.bam"
#     "HRR1849460.bam" "HRR1849462.bam"
# )

# # Get the BAM file corresponding to the array task ID
# BAM_FILE=${BAM_FILES[$SLURM_ARRAY_TASK_ID-1]}

# # Output directory based on the BAM file
# OUTPUT_DIR="./${BAM_FILE%.bam}/fastq"

# # Run the bamtofastq command
# /diazlab/data3/.abhinav/tools/bamtofastq_linux "$BAM_FILE" "$OUTPUT_DIR"

### We are Running cellranger on the fastq files
#### Analyzing the liver dataset
filenames <- list.files("/diazlab/data3/.abhinav/.immune/Kyoto/liver", pattern = "-liver", recursive = TRUE, include.dirs = T, full.names = TRUE)
write.table(filenames, "/diazlab/data3/.abhinav/.immune/Kyoto/liver/filenames.txt", sep = "\t", row.names = F, col.names = F, quote = F)

### Running the cell ranger
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --job-name=fastq
#SBATCH --time=1-00:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-user=abhinav.jain@ucsf.edu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --array=1-25

#### List the FASTQ directories
FASTQ_DIRS=(
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849454/fastq/control-liver1_0_1_H5T7CDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849454/fastq/control-liver1_0_1_H5YLGDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849454/fastq/control-liver1_0_1_HGWVWDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849455/fastq/control-liver2_0_1_H5GVWDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849455/fastq/control-liver2_0_1_H5T7CDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849456/fastq/control-liver3_0_1_H5T7CDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849456/fastq/control-liver3_0_1_H5YLGDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849456/fastq/control-liver3_0_1_HGWVWDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849457/fastq/control-liver4_0_1_H5T7CDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849457/fastq/control-liver4_0_1_HGWVWDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849458/fastq/control-liver5_0_1_H5T7CDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849458/fastq/control-liver5_0_1_HGWVWDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849459/fastq/pbc-liver1_0_1_H5T7CDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849459/fastq/pbc-liver1_0_1_H5YLGDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849459/fastq/pbc-liver1_0_1_HGWVWDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849460/fastq/pbc-liver2_0_1_H5T7CDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849460/fastq/pbc-liver2_0_1_HGWVWDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849461/fastq/pbc-liver3_0_1_H5T7CDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849461/fastq/pbc-liver3_0_1_HGWVWDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849462/fastq/pbc-liver4_0_1_H5T7CDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849462/fastq/pbc-liver4_0_1_H725MDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849462/fastq/pbc-liver4_0_1_HGWVWDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849463/fastq/pbc-liver5_0_1_H5T7CDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849463/fastq/pbc-liver5_0_1_H5TWNDSXY"
    "/diazlab/data3/.abhinav/.immune/Kyoto/liver/HRR1849463/fastq/pbc-liver5_0_1_HGWVWDSXY"
)

# Get the FASTQ directory corresponding to the array task ID
FASTQ_DIR=${FASTQ_DIRS[$SLURM_ARRAY_TASK_ID-1]}

# Correct the syntax for basename
idname=$(basename "$FASTQ_DIR")

# Run the cellranger count command
cellranger count --id $idname \
 --transcriptome /diazlab/data3/.abhinav/resources/10x_reference/refdata-gex-GRCh38-2024-A/ \
 --fastqs $FASTQ_DIR \
 --sample bamtofastq \
 --create-bam=false \
 --nosecondary


 ### Running the cellRanger aggregate
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --job-name=CR_aggr
#SBATCH --time=1-20:00:00
#SBATCH --output=./log_files/%x-%A_%a_aggr.out
#SBATCH --error=./log_files/%x-%A_%a_aggr.err
#SBATCH --mail-user=abhinav.jain@ucsf.edu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --array=0-9

# Load Cell Ranger module
module load cellranger

# Define array of CSV files
CSV_FILES=("control_liver1.csv" "control_liver2.csv" "control_liver3.csv" "control_liver4.csv" "control_liver5.csv"
    "pbc_liver1.csv" "pbc_liver2.csv" "pbc_liver3.csv" "pbc_liver4.csv" "pbc_liver5.csv")

# Get the specific CSV file based on the array task ID
CSV_FILE=${CSV_FILES[$SLURM_ARRAY_TASK_ID]}
SAMPLE_NAME=$(basename "$CSV_FILE" .csv)

echo "Running aggregation for $SAMPLE_NAME"

# Run cellranger aggr
cellranger aggr \
    --id="${SAMPLE_NAME}_aggregated" \
    --csv="$CSV_FILE" \
    --nosecondary

#### After aggregation Performing QC and combining and integration
### Loading the data and performing the QC on the object
library(Seurat)
samples <- c("control_liver1_aggregated", "control_liver2_aggregated","control_liver3_aggregated","control_liver4_aggregated", "control_liver5_aggregated",
            "pbc_liver1_aggregated","pbc_liver2_aggregated","pbc_liver3_aggregated","pbc_liver4_aggregated","pbc_liver5_aggregated")

source("/diazlab/data3/.abhinav/resources/all_scripts/R/scRNA_QC2.R")

for(i in 1:length(samples)){
    sample_obj <- scRNA_QC(Dir = "/diazlab/data3/.abhinav/.immune/Kyoto/liver/",Sample = samples[i],saveDir = "/diazlab/data3/.abhinav/.immune/Kyoto/liver/analysis/")
    assign(paste0(samples[i],"_obj"),sample_obj)
}

objname = ls(pattern = "_aggregated_obj")

## We have to combine the RNA singlet object to perform the normalization for the unwanted the effects.
combined <- merge(x=get(objname[1]), y = c(get(objname[2]),get(objname[3]),get(objname[4]),get(objname[5]),get(objname[6]),
                                               get(objname[7]),get(objname[8]),get(objname[9]),get(objname[10])),
                                               add.cell.ids = objname, project = "combined")

saveRDS(combined, "/diazlab/data3/.abhinav/.immune/Kyoto/liver/analysis/saveRDS/combined.RDS")

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
combined <- RunHarmony(combined,group.by.vars = "orig.ident")

combined <- RunUMAP(
  combined,
  reduction.key = "harmonyUMAP_",
  reduction = "harmony",
  reduction.name = "harmonyumap",
  dims = 1:30
)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/liver/"
dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/combined_samples.pdf"), width = 6, height = 5)
DimPlot(combined, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.4)

dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/combined_cluster.pdf"), width = 6, height = 5)
DimPlot(combined, reduction = "harmonyumap", group.by = "seurat_clusters", label = TRUE, label.size = 4)
dev.off()

dir.create(paste0(savedir,"analysis/featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/featureplots/Tcell_Bcells.pdf"), width = 10)
FeaturePlot(combined, c("CD3E","CD3D","CD3G","CD4","CD8A","CD8B","CD19","CD20","CD79A","CD22","CD27"), reduction = "harmonyumap")
dev.off()

dir.create(paste0(savedir,"analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Tcell_Bcells_nopoint.pdf"), width = 12, height = 5)
VlnPlot(combined, c("CD3E","CD3D","CD3G","CD4","CD8A","CD8B","CD19","CD20","CD79A","CD22","CD27"), pt.size = 0)
dev.off()

dir.create(paste0(savedir,"analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Tcell_Bcells_nopoint.pdf"), width = 12, height = 5)
VlnPlot(combined, c("CD3E","CD3D","CD3G","CD4","CD8A","CD8B","CD19","CD20","CD79A","CD22","CD27"), pt.size = 0)
dev.off()s

### Running with the join layers
# https://satijalab.org/seurat/articles/seurat5_integration
combined <- IntegrateLayers(
  object = combined, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE
)

combined <- FindNeighbors(combined, reduction = "integrated.rpca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.8, cluster.name = "rpca_clusters")
combined <- RunUMAP(combined, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")


dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/combined_sample_rpca.pdf"), width = 6, height = 5)
DimPlot(combined, reduction = "umap.rpca", group.by = "orig.ident", label = TRUE, label.size = 4)
dev.off()

pdf(paste0(savedir, "analysis/UMAP/combined_cluster_rpca.pdf"), width = 6, height = 5)
DimPlot(combined, reduction = "umap.rpca", group.by = "rpca_clusters", label = TRUE, label.size = 4)
dev.off()

dir.create(paste0(savedir,"analysis/featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/featureplots/Tcell_Bcells_rpca.pdf"), width = 10)
FeaturePlot(combined, c("CD3E","CD3D","CD3G","CD4","CD8A","CD8B","CD19","CD20","CD79A","CD22","CD27"), reduction = "umap.rpca")
dev.off()

dir.create(paste0(savedir,"analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Tcell_Bcells_nopoint_rpca.pdf"), width = 12, height = 5)
VlnPlot(combined, c("CD3E","CD3D","CD3G","CD4","CD8A","CD8B","CD19","CD20","CD79A","CD22","CD27"), pt.size = 0, group.by = "rpca_clusters")
dev.off()

combined <- JoinLayers(combined)
combined<- NormalizeData(combined)
markers <- FindAllMarkers(combined)
dir.create(paste0(savedir,"analysis/Table"), showWarnings = FALSE)
write.table(markers, paste0(savedir,"analysis/Table/marker_rpca.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

Idents(combined) <- combined@meta.data$RNA_snn_res.0.4
markers <- FindAllMarkers(combined, group.by = "RNA_snn_res.0.4")
dir.create(paste0(savedir,"analysis/Table"), showWarnings = FALSE)
write.table(markers, paste0(savedir,"analysis/Table/marker_harmony.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

saveRDS(combined, "/diazlab/data3/.abhinav/.immune/Kyoto/liver/analysis/saveRDS/combined_integrated.RDS")
#endregion

#region Tcells
# Yuki emails
# T cell:  Using CD3D/E/G, CD2 as T cell markers, then clustering for sorting out T cells
# Subsequent analysis: CD4, CD8A, IL7R, TCF7, LEF1, TBX21, EOMES, BCL6, CD28, FAS, TNFSF8, IL21, IFNG, GZMK
# B cells: CD19, CD79A, MS4A1 for sorting out B cells
# Subsequent analysis: CD27, IGHD, ITGAX, FCRL5,  IGHM, IGHG, CD38, CD21
library(Seurat)
library(Harmony)
combined = readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/liver/analysis/saveRDS/combined_integrated.RDS")
combined@meta.data$seurat_clusters <- (combined@meta.data$RNA_snn_res.0.4)
combined@meta.data$seurat_clusters <- factor(combined@meta.data$seurat_clusters, levels = c(0:20))
Idents(combined) <- factor(Idents(combined), levels = c(0:20))

pdf("/diazlab/data3/.abhinav/.immune/Kyoto/liver/analysis/UMAP/combined_UMAP.pdf")
DimPlot(combined, group.by = "seurat_clusters", label = TRUE) + NoLegend()
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/liver/"
pdf(paste0(savedir, "analysis/vlnplots/Tcells.pdf"), width = 12, height = 5)
VlnPlot(combined, c("CD3E","CD3D","CD3G","CD2"), pt.size = 0)
dev.off()

Tcells_clus = c(0,1,2,4,5,18)
Tcells <- subset(combined, idents = Tcells_clus)

Tcells <- NormalizeData(Tcells)
Tcells <- FindVariableFeatures(Tcells)
Tcells <- ScaleData(Tcells)
Tcells <- RunPCA(Tcells)
Tcells <- RunHarmony(Tcells,group.by.vars = "orig.ident")

Tcells <- RunUMAP(
  Tcells,
  reduction.key = "harmonyUMAP_",
  reduction = "harmony",
  reduction.name = "harmonyumap",
  dims = 1:30
)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/"
dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/combined_samples.pdf"), width = 6, height = 5)
DimPlot(Tcells, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

Tcells <- FindNeighbors(Tcells, reduction = "harmony", dims = 1:30)
Tcells <- FindClusters(Tcells, resolution = 0.8)
Tcell_obj <- FindClusters(Tcell_obj, resolution = 0.4)

dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/Tcells_cluster_0.4_2.pdf"), width = 6, height = 5)
DimPlot(Tcell_obj, reduction = "harmonyumap", group.by = "seurat_clusters", label = TRUE, label.size = 4)
dev.off()

genes <- c("CD3E","CD3G","CD3D","CD2","CD4", "CD8A","CD8B", "IL7R", 'TCF7', "LEF1", 'TBX21', "EOMES", "BCL6", "CD28", "FAS", "TNFSF8", "IL21", "IFNG", "GZMK")
genes <- c("TNF", "FOXP3", "MKI67", "CXCR6", "TYROBP", "LYZ")

dir.create(paste0(savedir,"analysis/featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/featureplots/Tcell_required.pdf"))
FeaturePlot(Tcells, genes, reduction = "harmonyumap",label = TRUE)
dev.off()

dir.create(paste0(savedir,"analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Tcell_nopoint_required.pdf"), width = 12, height = 10)
VlnPlot(Tcells, genes, pt.size = 0, group.by = "seurat_clusters")
dev.off()

Set1 = c("TCF7", "LEF1", "BACH2", "MAF", "TOX", "ZEB2")
Set2 = c("TNF", "IFNG", "CD40LG", "XCL1", "CCL5", "STAT1")
Set3 = c("HSPA1A", "HSPA1B", "HSPH1", "NR4A1", "DNAJA1", "HSPD1")
Set4 = c("JUN", "JUNB", "FOS", "FOSB", "REL", "NFKB1")
Set5 = c("CD4_Naïve", "CD4_Treg_signature", "CD4_Activation_Effector_function", "CD4_Adhesion", "CD4_TCR_signaling", "CD4_Stress_response")  

set = c("Set1","Set2","Set3","Set4","Set5")
savedir = "/mnt/data/projects/liver/"
dir.create(paste0(savedir,"analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
for(i in 1:length(set)){
  pdf(paste0(savedir, "analysis/vlnplots/",set[i],"_nopoint_required.pdf"), width = 12, height = 8)
  print(VlnPlot(CD4_obj, get(set[i]), pt.size = 0, group.by = "seurat_clusters_new"))
  dev.off()
}

Tcellmarkers <- FindAllMarkers(Tcells)
write.table(Tcellmarkers, "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/analysis/Tcellmarker_res0.4.txt", sep = "\t", quote = F, row.names = T, col.names = T)

Tcells <- ScaleData(Tcells, features=rownames(Tcells))

dir.create(paste0(savedir, "analysis/heatmap/"),showWarnings = FALSE)
pdf(paste0(savedir, "analysis/heatmap/Tcell_res0.4.pdf"), height = 4)
DoHeatmap(Tcells, genes)
dev.off()

dir.create(paste0(savedir, "saveRDS"), showWarnings = FALSE)
saveRDS(Tcells, paste0(savedir, "saveRDS/Tcells.RDS"))
#endregion

#region removing macrophages from T cells
Tcell_obj = readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/saveRDS/Tcells.RDS")

mps <- c(3,6,8,10,12)
Tcells_clus <- levels(Tcell_obj@meta.data$seurat_clusters)[grep(paste0("^",mps,"$",collapse = "|"),levels(Tcell_obj@meta.data$seurat_clusters),invert=TRUE)]
Tcell_subset_obj = subset(Tcell_obj, idents = Tcells_clus)

Tcell_subset_obj <- NormalizeData(Tcell_subset_obj)
Tcell_subset_obj <- FindVariableFeatures(Tcell_subset_obj)
Tcell_subset_obj <- ScaleData(Tcell_subset_obj)
Tcell_subset_obj <- RunPCA(Tcell_subset_obj)
Tcell_subset_obj <- RunHarmony(Tcell_subset_obj,group.by.vars = "orig.ident")

Tcell_subset_obj <- RunUMAP(
  Tcell_subset_obj,
  reduction.key = "harmonyUMAP_",
  reduction = "harmony",
  reduction.name = "harmonyumap",
  dims = 1:30
)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_mps/"
dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/Tcells_samples.pdf"), width = 6, height = 5)
DimPlot(Tcell_subset_obj, reduction = "harmonyumap", label = TRUE, label.size = 4, group.by =  "orig.ident")
dev.off()

Tcell_subset_obj <- FindNeighbors(Tcell_subset_obj, reduction = "harmony", dims = 1:30)
# Tcell_subset_obj <- FindClusters(Tcell_subset_obj, resolution = 0.8)
Tcell_subset_obj <- FindClusters(Tcell_subset_obj, resolution = 0.5)

dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/Tcell_subset_obj_cluster_0.4_2.pdf"), width = 6, height = 5)
DimPlot(Tcell_subset_obj, reduction = "harmonyumap", group.by = "seurat_clusters", label = TRUE, label.size = 4) + geom_vline(xintercept = -4, linetype = "dashed", color = "red") 
dev.off()

### cell required
required_cells <- names(Tcell_subset_obj@reductions$harmonyumap@cell.embeddings[,"harmonyUMAP_1"] > -4 )

Tcell_subset_obj_2 = subset(Tcell_subset_obj, cells = required_cells)
Tcell_subset_obj_2 <- NormalizeData(Tcell_subset_obj_2)
Tcell_subset_obj_2 <- FindVariableFeatures(Tcell_subset_obj_2)
Tcell_subset_obj_2 <- ScaleData(Tcell_subset_obj_2)
Tcell_subset_obj_2 <- RunPCA(Tcell_subset_obj_2)
Tcell_subset_obj_2 <- RunHarmony(Tcell_subset_obj_2,group.by.vars = "orig.ident")

Tcell_subset_obj_2 <- RunUMAP(
  Tcell_subset_obj_2,
  reduction.key = "harmonyUMAP_",
  reduction = "harmony",
  reduction.name = "harmonyumap",
  dims = 1:30
)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_mps/"
dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/Tcells_samples_remove_outlier_cells.pdf"), width = 6, height = 5)
DimPlot(Tcell_subset_obj_2, reduction = "harmonyumap", label = TRUE, label.size = 4, group.by =  "orig.ident")
dev.off()

Tcell_subset_obj_2 <- FindNeighbors(Tcell_subset_obj_2, reduction = "harmony", dims = 1:30)
# Tcell_subset_obj_2 <- FindClusters(Tcell_subset_obj_2, resolution = 0.8)
Tcell_subset_obj_2 <- FindClusters(Tcell_subset_obj_2, resolution = 0.5)

dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/Tcell_subset_obj_2_cluster_0.4_remove_outlier_cells.pdf"), width = 6, height = 5)
DimPlot(Tcell_subset_obj_2, reduction = "harmonyumap", group.by = "seurat_clusters", label = TRUE, label.size = 4)
dev.off()


genes <- c("CD3E","CD3G","CD3D","CD2","CD4", "CD8A","CD8B", "IL7R", 'TCF7', "LEF1", 'TBX21', "EOMES", "BCL6", "CD28", "FAS", "TNFSF8", "IL21", "IFNG", "GZMK","TNF", "FOXP3", "MKI67", "CXCR6", "TYROBP", "LYZ")

dir.create(paste0(savedir,"analysis/featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/featureplots/Tcell_required.pdf"), width = 15, height = 15)
FeaturePlot(Tcell_subset_obj_2, genes, reduction = "harmonyumap",label = TRUE)
dev.off()

dir.create(paste0(savedir,"analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Tcell_nopoint_required.pdf"), width = 12, height = 12)
VlnPlot(Tcell_subset_obj_2, genes, pt.size = 0, group.by = "seurat_clusters")
dev.off()

Tcellmarkers <- FindAllMarkers(Tcell_subset_obj_2)
write.table(Tcellmarkers, "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/analysis/Tcellmarker_res0.4.txt", sep = "\t", quote = F, row.names = T, col.names = T)

Tcell_subset_obj_2 <- ScaleData(Tcell_subset_obj_2, features=rownames(Tcell_subset_obj_2))

dir.create(paste0(savedir, "analysis/heatmap/"),showWarnings = FALSE)
pdf(paste0(savedir, "analysis/heatmap/Tcell_res0.4.pdf"), height = 4)
DoHeatmap(Tcell_subset_obj_2, genes)
dev.off()

dir.create(paste0(savedir, "analysis/saveRDS"), showWarnings = FALSE)
saveRDS(Tcell_subset_obj_2, paste0(savedir, "saveRDS/analysis/Tcell_subset_obj_2.RDS"))

#endregion

#region remove clus 3, 10, and 12
savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_mps_clus3_10_12/"
dir.create(savedir, showWarnings = FALSE)

mps <- c(3,10,12)
Tcells_clus <- levels(Tcell_subset_obj_2@meta.data$seurat_clusters)[grep(paste0("^",mps,"$",collapse = "|"),levels(Tcell_subset_obj_2@meta.data$seurat_clusters),invert=TRUE)]
Tcell_subset_obj_3 = subset(Tcell_subset_obj_2, idents = Tcells_clus)

Tcell_subset_obj_3 <- NormalizeData(Tcell_subset_obj_3)
Tcell_subset_obj_3 <- FindVariableFeatures(Tcell_subset_obj_3)
Tcell_subset_obj_3 <- ScaleData(Tcell_subset_obj_3)
Tcell_subset_obj_3 <- RunPCA(Tcell_subset_obj_3)
Tcell_subset_obj_3 <- RunHarmony(Tcell_subset_obj_3,group.by.vars = "orig.ident")

Tcell_subset_obj_3 <- RunUMAP(
  Tcell_subset_obj_3,
  reduction.key = "harmonyUMAP_",
  reduction = "harmony",
  reduction.name = "harmonyumap",
  dims = 1:30
)

dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/Tcells_samples_remove_outlier_cells.pdf"), width = 6, height = 5)
DimPlot(Tcell_subset_obj_3, reduction = "harmonyumap", label = TRUE, label.size = 4, group.by =  "orig.ident")
dev.off()

pdf(paste0(savedir, "analysis/UMAP/Tcells_samples_remove_outlier_cells_splitted.pdf"))
DimPlot(Tcell_subset_obj_3, reduction = "harmonyumap", label = FALSE, label.size = 4, group.by =  "orig.ident", split.by = "orig.ident")
dev.off()


Tcell_subset_obj_3 <- FindNeighbors(Tcell_subset_obj_3, reduction = "harmony", dims = 1:30)
Tcell_subset_obj_3 <- FindClusters(Tcell_subset_obj_3, resolution = 1.2)
# Tcell_subset_obj_3 <- FindClusters(Tcell_subset_obj_3, resolution = 0.5)

dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/Tcell_subset_obj_3_cluster_0.8_remove_outlier_cells.pdf"), width = 6, height = 5)
DimPlot(Tcell_subset_obj_3, reduction = "harmonyumap", group.by = "seurat_clusters", label = TRUE, label.size = 4)
dev.off()


genes <- c("CD3E","CD3G","CD3D","CD2","CD4", "CD8A","CD8B", "IL7R", 'TCF7', "LEF1", 'TBX21', "EOMES", "BCL6", "CD28", "FAS", "TNFSF8", "IL21", "IFNG", "GZMK","TNF", "FOXP3", "MKI67", "CXCR6", "TYROBP", "LYZ")

dir.create(paste0(savedir,"analysis/featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/featureplots/Tcell_required_1.2.pdf"), width = 15, height = 15)
FeaturePlot(Tcell_subset_obj_3, genes, reduction = "harmonyumap",label = TRUE)
dev.off()

dir.create(paste0(savedir,"analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Tcell_nopoint_required.pdf"), width = 15, height = 15)
VlnPlot(Tcell_subset_obj_3, genes, pt.size = 0, group.by = "seurat_clusters")
dev.off()

Tcellmarkers <- FindAllMarkers(Tcell_subset_obj_3)
write.table(Tcellmarkers, "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/analysis/Tcellmarker_res0.4.txt", sep = "\t", quote = F, row.names = T, col.names = T)

Tcell_subset_obj_3 <- ScaleData(Tcell_subset_obj_3, features=rownames(Tcell_subset_obj_3))

dir.create(paste0(savedir, "analysis/heatmap/"),showWarnings = FALSE)
pdf(paste0(savedir, "analysis/heatmap/Tcell_res0.4.pdf"), height = 4)
DoHeatmap(Tcell_subset_obj_3, genes)
dev.off()

dir.create(paste0(savedir, "analysis/saveRDS"), showWarnings = FALSE)
saveRDS(Tcell_subset_obj_3, paste0(savedir, "analysis/saveRDS/Tcell_subset_obj_3.RDS"))
#endregion

#region removing more mps
savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/"
dir.create(savedir, showWarnings = FALSE)

mps <- c(10,13,16,17)
Tcells_clus <- levels(Tcell_subset_obj_3@meta.data$seurat_clusters)[grep(paste0("^",mps,"$",collapse = "|"),levels(Tcell_subset_obj_3@meta.data$seurat_clusters),invert=TRUE)]

Tcell_subset_obj_4 = subset(Tcell_subset_obj_3, idents = Tcells_clus)
Tcell_subset_obj_4 <- NormalizeData(Tcell_subset_obj_4)
Tcell_subset_obj_4 <- FindVariableFeatures(Tcell_subset_obj_4)
Tcell_subset_obj_4 <- ScaleData(Tcell_subset_obj_4)
Tcell_subset_obj_4 <- RunPCA(Tcell_subset_obj_4)
Tcell_subset_obj_4 <- RunHarmony(Tcell_subset_obj_4,group.by.vars = "orig.ident")

Tcell_subset_obj_4 <- RunUMAP(
  Tcell_subset_obj_4,
  reduction.key = "harmonyUMAP_",
  reduction = "harmony",
  reduction.name = "harmonyumap",
  dims = 1:30
)

dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/Tcells_samples_remove_outlier_cells.pdf"), width = 6, height = 5)
DimPlot(Tcell_subset_obj_4, reduction = "harmonyumap", label = TRUE, label.size = 4, group.by =  "orig.ident")
dev.off()

pdf(paste0(savedir, "analysis/UMAP/Tcells_samples_remove_outlier_cells_splitted.pdf"))
DimPlot(Tcell_subset_obj_4, reduction = "harmonyumap", label = FALSE, label.size = 4, group.by =  "orig.ident", split.by = "orig.ident", ncol = 5)
dev.off()


Tcell_subset_obj_4 <- FindNeighbors(Tcell_subset_obj_4, reduction = "harmony", dims = 1:30)
Tcell_subset_obj_4 <- FindClusters(Tcell_subset_obj_4, resolution = 1.2)
# Tcell_subset_obj_4 <- FindClusters(Tcell_subset_obj_4, resolution = 0.5)

dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/Tcell_subset_obj_4_cluster_0.8_remove_outlier_cells.pdf"), width = 6, height = 5)
DimPlot(Tcell_subset_obj_4, reduction = "harmonyumap", group.by = "seurat_clusters", label = TRUE, label.size = 4)
dev.off()


genes <- c("CD3E","CD3G","CD3D","CD2","CD4", "CD8A","CD8B", "IL7R", 'TCF7', "LEF1", 'TBX21', "EOMES", "BCL6", "CD28", "FAS", "TNFSF8", "IL21", "IFNG", "GZMK","TNF", "FOXP3", "MKI67", "CXCR6", "TYROBP", "LYZ")

dir.create(paste0(savedir,"analysis/featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/featureplots/Tcell_required_1.2.pdf"), width = 15, height = 15)
FeaturePlot(Tcell_subset_obj_4, genes, reduction = "harmonyumap",label = TRUE)
dev.off()

dir.create(paste0(savedir,"analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Tcell_nopoint_required.pdf"), width = 15, height = 15)
VlnPlot(Tcell_subset_obj_4, genes, pt.size = 0, group.by = "seurat_clusters")
dev.off()

Tcellmarkers <- FindAllMarkers(Tcell_subset_obj_4)
write.table(Tcellmarkers, paste0(savedir,"analysis/Tcellmarker_res1.2.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

Tcell_subset_obj_4 <- ScaleData(Tcell_subset_obj_4, features=rownames(Tcell_subset_obj_4))

dir.create(paste0(savedir, "analysis/heatmap/"),showWarnings = FALSE)
pdf(paste0(savedir, "analysis/heatmap/Tcell_res0.4.pdf"), height = 4)
DoHeatmap(Tcell_subset_obj_4, genes)
dev.off()

#### Checking the quantitative difference between the cluster pbc and control 
#### Longitudinal analysis
celltypes_ind <- table(Tcell_subset_obj_4@meta.data$orig.ident,Tcell_subset_obj_4@meta.data$seurat_clusters)
dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(celltypes_ind, paste0(savedir,"Table/celltype_individual_res_1.2.txt"),
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

df_melted$condition <- gsub("_.*.","",df_melted$sample)

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=sample)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table//marker_celltype_res_0.8_pts_timepoints.pdf"),width =10, height = 7)
p
dev.off()

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=condition)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table//marker_celltype_res_0.8_condition.pdf"),width =10, height = 7)
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
  stat.test,  label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

pdf(paste(savedir,"Table/t_test_unpaired_all_samples_bxplot_two_sided_clus_res_1.4.pdf",sep = ""), width = 8, height = 6)
bxp4
dev.off()

### Adding genes
DefaultAssay(obj) <- "MAGIC_RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/.immune/Kyoto/liver/resources/genelist/", pattern = "CD", full.names = TRUE)
filename <- basename(files)
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[, 1]
    Tcellsubset <- rownames(obj)[match(Tcellsubset, rownames(obj), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

obj <- AddModuleScore(obj, features = markers, slot = "data")
colnames(obj@meta.data)[15:25] <- filename

savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/analysis/"

dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir,"vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/module_score_genes.pdf"), width = 12, height = 12)
VlnPlot(obj, filename, pt.size = 0)
dev.off()

pdf(paste0(savedir, "featureplots/module_score_genes.pdf"), width = 12, height = 12)
FeaturePlot(obj, filename, reduction = "harmonyumap")
dev.off()

saveRDS(obj, "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/analysis/saveRDS/Tcell_subset_obj_4.RDS")

# PBC_Tcell = readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/analysis/saveRDS/Tcell_subset_obj_4.RDS")
DefaultAssay(PBC_Tcell) <- "RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/.immune/Kyoto/resource/CD4_genes/", pattern = "txt", full.names = TRUE)
filename <- gsub(".txt", "", basename(files))
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[, 1]
    Tcellsubset <- rownames(PBC_Tcell)[match(Tcellsubset, rownames(PBC_Tcell), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

PBC_Tcell <- AddModuleScore(PBC_Tcell, features = markers, slot = "data")

colnames(PBC_Tcell@meta.data)[15:31] <- paste0("CD4_", filename)

# Idents(PBC_Tcell) <- PBC_Tcell@meta.data$seurat_clusters2

# savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/CD4/"
savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/analysis/"
dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir, "vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/CD4_pancancer_PBC_Tcells_module_score_genes_violin.pdf"), width = 12, height = 12)
VlnPlot(PBC_Tcell, paste0("CD4_", filename), pt.size = 0)
dev.off()

dir.create(paste0(savedir, "featureplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/CD4_pancancer_PBC_Tcells_module_score_genes_featureplot.pdf"), width = 12, height = 12)
FeaturePlot(PBC_Tcell, paste0("CD4_", filename), reduction = "harmonyumap", label = TRUE)
dev.off()

### CD8
DefaultAssay(PBC_Tcell) <- "RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/.immune/Kyoto/resource/CD8_genes/", pattern = "txt", full.names = TRUE)
filename <- gsub(".txt", "", basename(files))
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[, 1]
    Tcellsubset <- rownames(PBC_Tcell)[match(Tcellsubset, rownames(PBC_Tcell), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

PBC_Tcell <- AddModuleScore(PBC_Tcell, features = markers, slot = "data")
colnames(PBC_Tcell@meta.data)[32:50] <- paste0("CD8_", filename)

# savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC/"
dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir, "vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/CD8_pancancer_PBC_Tcells_module_score_genes_violin.pdf"), width = 12, height = 12)
VlnPlot(PBC_Tcell, paste0("CD8_", filename), pt.size = 0)
dev.off()

dir.create(paste0(savedir, "featureplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/CD8_pancancer_module_score_genes_featureplot.pdf"), width = 12, height = 12)
FeaturePlot(PBC_Tcell, paste0("CD8_", filename), reduction = "harmonyumap", label = TRUE)
dev.off()

dir.create(paste0(savedir, "UMAP/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "UMAP/Tcell_PBC.pdf"))
DimPlot(PBC_Tcell, reduction = "harmonyumap", label = TRUE)
dev.off()


dir.create(paste0(savedir, "analysis/saveRDS"), showWarnings = FALSE)
saveRDS(Tcell_subset_obj_4, paste0(savedir, "analysis/saveRDS/Tcell_subset_obj_4.RDS"))




#endregion

#region T cell Annotation

clusters = 0:14
celltypes = c(
  "GZMK+CXCR6-CD69- memory CD8", 
      "Th1/Cyto", 
      "CXCR6+ TRM",
      "CX3CR1+ CD8 CTL",
      "Tph-like",
      "Exhausted CD8 ",
      "MAIT cells",
      "post TCF1hi stem CD4",
      "TCF1hi stem CD4",
      "Treg",
      "Dual producing Th1",
      "TCF1hi stem-like CD8",
      "gamma delta",
      "proliferating",
      "CXCR6+ TRM 2"
)

patterns <- clusters
replacements <- celltypes

library(stringr)
obj@meta.data$celltype <- obj@meta.data$seurat_clusters
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    obj@meta.data$celltype <- str_replace_all(obj@meta.data$celltype, pattern, replacements[i])
}

pdf(paste0(savedir,"UMAP/celltype.pdf"), width = 10, height = 5.6)
DimPlot(obj, group.by = "celltype", reduction = "harmonyumap", label = TRUE, 
label.size = 4)
dev.off()

pdf(paste0(savedir,"UMAP/celltype_nolabel.pdf"), width = 10, height = 5.6)
DimPlot(obj, group.by = "celltype", reduction = "harmonyumap", label = FALSE)
dev.off()

#### Making a heatmap
Idents(obj) <- obj@meta.data$celltype
celltype.markers <- FindAllMarkers(obj, only.pos = TRUE, group.by = "celltype", assay = "RNA")
write.table(celltype.markers, paste0(savedir, "Table/celltype_markers.txt"),
quote = F, row.names = F, col.names = T, sep = "\t")

library(dplyr)
celltype.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5)

celltype.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5) %>%
    slice_head(n = 5) %>%
    ungroup() -> top10


dir.create(paste0(savedir, "heatmap"), showWarnings = FALSE)
pdf(paste0(savedir, "heatmap/celltype_markers_top5.pdf"),height = 12, width = 8)
DoHeatmap(obj, features = top10$gene) + NoLegend()
dev.off()

#### Checking the quantitative difference between the cluster pbc and control 
#### Longitudinal analysis
celltypes_ind <- table(obj@meta.data$orig.ident,obj@meta.data$celltype)
dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(celltypes_ind, paste0(savedir,"Table/celltype_individual.txt"),
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

df_melted$condition <- gsub("_.*.","",df_melted$sample)

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=sample)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table//marker_celltype_res_0.8_pts_timepoints.pdf"),width =10, height = 7)
p
dev.off()

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=condition)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table//marker_celltype_res_0.8_condition.pdf"),width =10, height = 7)
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
  stat.test,  label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

pdf(paste(savedir,"Table/t_test_unpaired_all_samples_bxplot.pdf",sep = ""), width = 8, height = 6)
bxp4
dev.off()

#endregion

#region extracting CD4
savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/CD4/"
dir.create(savedir, showWarnings = FALSE)

CD4 <- c(1,4,7,8,9,10)
CD4_cluster <- levels(obj@meta.data$seurat_clusters)[grep(paste0("^",CD4,"$",collapse = "|"),levels(obj@meta.data$seurat_clusters),invert=FALSE)]

Idents(obj) <- obj@meta.data$seurat_clusters
CD4_obj = subset(obj, idents = CD4_cluster)
CD4_obj <- NormalizeData(CD4_obj)
CD4_obj <- FindVariableFeatures(CD4_obj)
CD4_obj <- ScaleData(CD4_obj)
CD4_obj <- RunPCA(CD4_obj)
CD4_obj <- RunHarmony(CD4_obj,group.by.vars = "orig.ident")

CD4_obj <- RunUMAP(
  CD4_obj,
  reduction.key = "harmonyUMAP_",
  reduction = "harmony",
  reduction.name = "harmonyumap",
  dims = 1:30
)

dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "UMAP/CD4_Tcell_umap.pdf"), width = 6, height = 5)
DimPlot(CD4_obj, reduction = "harmonyumap", label = TRUE, label.size = 4, group.by =  "orig.ident")
dev.off()

pdf(paste0(savedir, "UMAP/CD4_Tcell_splitted.pdf"))
DimPlot(CD4_obj, reduction = "harmonyumap", label = FALSE,
        label.size = 4, group.by =  "orig.ident",
        split.by = "orig.ident", ncol = 5)
dev.off()

dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "UMAP/CD4_Tcell_cluster.pdf"), width = 6, height = 5)
DimPlot(CD4_obj, reduction = "harmonyumap", 
group.by = "seurat_clusters", label = TRUE, label.size = 8)
dev.off()

CD4_obj@meta.data$seurat_clusters2 <- CD4_obj@meta.data$seurat_clusters

CD4_obj <- FindNeighbors(CD4_obj, reduction = "harmony", dims = 1:30)
# Tcell_subset_obj_4 <- FindClusters(Tcell_subset_obj_4, resolution = 1.2)
CD4_obj <- FindClusters(CD4_obj, resolution = 0.6)

pdf(paste0(savedir, "UMAP/CD4_Tcell_cluster2.pdf"), width = 6, height = 5)
DimPlot(CD4_obj, reduction = "harmonyumap", 
group.by = "seurat_clusters", label = TRUE, label.size = 8)
dev.off()

genes <- c("TCF7", "LEF1", "IL7R", "GPR183", "LTB", "CD28", "ICOS", "MAF", "TNFSF8", "CXCL13", "IL21", 
"CSF2", "IFNG", "TNF", "FOXP3", "IL2RA", "ENTPD1", "CTLA4", "PDCD1", "TIGIT", "HAVCR2", "GZMK", 
"CXCR5", "CXCR6", "CXCR3", "CCR7")

dir.create(paste0(savedir,"featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/Tcell_required_1.2.pdf"), width = 15, height = 15)
FeaturePlot(CD4_obj, genes, reduction = "harmonyumap",label = TRUE)
dev.off()

dir.create(paste0(savedir,"vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/Tcell_nopoint_required.pdf"), width = 15, height = 15)
VlnPlot(CD4_obj, genes, pt.size = 0, group.by = "seurat_clusters")
dev.off()

Tcellmarkers <- FindAllMarkers(CD4_obj)
write.table(Tcellmarkers, paste0(savedir,"Tcellmarker.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

CD4_obj <- ScaleData(CD4_obj, features=rownames(CD4_obj))

dir.create(paste0(savedir, "heatmap/"),showWarnings = FALSE)
pdf(paste0(savedir, "heatmap/CD4_Tcell.pdf"), height = 12, width = 8)
DoHeatmap(CD4_obj, genes)
dev.off()

savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/CD4/"
dir.create(paste(savedir,"saveRDS", sep = ""), showWarnings = FALSE)

pdf(paste0(savedir, "vlnplots/Tcell_nopoint_required_6_clus.pdf"), width = 15, height = 15)
VlnPlot(CD4_obj, genes, pt.size = 0, group.by = "seurat_clusters2")
dev.off()

### Adding genes
DefaultAssay(CD4_obj) <- "RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/.immune/Kyoto/liver/resources/genelist/", pattern = "CD", full.names = TRUE)
filename <- paste0(basename(files),"_2")
for (i in 1:length(files)) { 
    Tcellsubset <- read.table(files[i])[, 1]
    Tcellsubset <- rownames(obj)[match(Tcellsubset, rownames(obj), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

CD4_obj <- AddModuleScore(CD4_obj, features = markers, slot = "data")
colnames(CD4_obj@meta.data)[29:39] <- filename

dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir,"vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/module_score_genes_6_clusters.pdf"), width = 12, height = 12)
VlnPlot(CD4_obj, filename, pt.size = 0, group.by = "seurat_clusters2")
dev.off()

pdf(paste0(savedir, "featureplots/module_score_genes_6_clusters.pdf"), width = 12, height = 12)
FeaturePlot(CD4_obj, filename, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

pdf(paste0(savedir, "featureplots/Genes_6_clusters.pdf"), width = 12, height = 12)
FeaturePlot(CD4_obj, genes, reduction = "harmonyumap", label = TRUE, label.size = 4)
dev.off()

cluster=c(1,4,7,8,9,10) 
celltypes=c("GZMK Th1","Tph","Post stem-like CD4 T cells 1","Stem-like CD4 T cells","Treg","Dual producing Th1")

patterns <- cluster
replacements <- celltypes

library(stringr)
liver_CD4@meta.data$celltypes <- liver_CD4@meta.data$seurat_clusters2
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    liver_CD4@meta.data$celltypes <- str_replace_all(liver_CD4@meta.data$celltypes, pattern, replacements[i])
}

savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/CD4/"
pdf(paste0(savedir, "UMAP/PBC_celltype_CD4_HCC_Yuki_annot.pdf"))
DimPlot(liver_CD4, group.by = "celltypes", reduction = "harmonyumap", label = TRUE, label.size = 6) + NoLegend()
dev.off()

saveRDS(CD4_obj, paste0(savedir,"saveRDS/CD4_obj.RDS"))
#endregion

#region CD4 Trajectory
CD4_obj=readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/CD4/saveRDS/CD4_obj.RDS")

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
DefaultAssay(CD4_obj) <- "RNA"
sce <- as.SingleCellExperiment(CD4_obj, assay = "RNA")

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
savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/CD4/trajectory/"
dir.create(savedir, showWarnings = FALSE)
dir.create(paste(savedir, "UMAP/", sep = ""), showWarnings = FALSE)
p <- DimPlot(CD4_obj, reduction = "harmonyumap", label = TRUE, group.by = "seurat_clusters2")
pdf(paste(savedir, "UMAP/combine.pdf", sep = ""))
p
dev.off()

CD4_sce <- slingshot(sce,
  clusterLabels = "seurat_clusters2",
  reducedDim = "UMAP", approx_points = 100,
  omega = TRUE,
  omega_scale = 1.5
)

CD4_sce_naive <- slingshot(sce,
  clusterLabels = colData(sce)$seurat_clusters2,
  reducedDim = "HARMONYUMAP",
  start.clus = c("8"),
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
embedded_orig <- embedCurves(CD4_sce_naive, "HARMONYUMAP")

rm(plot_list)
plot_list <- list()
for (i in 1:3) {
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
pdf(paste(savedir, "UMAP/Liver_CD4_clus_3lineages.pdf", sep = ""), width = 12, height = 3.5)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], nrow = 1, ncol = 3)
dev.off()

pdf(paste(savedir, "UMAP/Liver_CD4_UMAP_lineage.pdf", sep = ""), width = 4.5, height = 4)
plot_list
dev.off()

#region Bcells
# B cells: CD19, CD79A, MS4A1 for sorting out B cells
# Subsequent analysis: CD27, IGHD, ITGAX, FCRL5,  IGHM, IGHG, CD38, CD21

library(Seurat)
library(harmony)
combined = readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/liver/analysis/saveRDS/combined_integrated.RDS")
combined@meta.data$seurat_clusters <- (combined@meta.data$RNA_snn_res.0.4)
combined@meta.data$seurat_clusters <- factor(combined@meta.data$seurat_clusters, levels = c(0:20))
Idents(combined) <- factor(Idents(combined), levels = c(0:20))
pdf("/diazlab/data3/.abhinav/.immune/Kyoto/liver/analysis/UMAP/combined_UMAP_Bcells.pdf")
DimPlot(combined, group.by = "seurat_clusters", label = TRUE) + NoLegend()
dev.off()

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Bcells/"
dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir,"analysis/vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Bcells.pdf"), width = 9, height = 3)
VlnPlot(combined, c("CD19", "CD79A", "MS4A1"), pt.size = 0)
dev.off()

dir.create(paste0(savedir,"analysis/featureplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/featureplots/Bcells.pdf"), width = 9, height = 4)
FeaturePlot(combined, c("CD19", "CD79A", "MS4A1"), reduction = "harmonyumap")
dev.off()

Bcells_clus = c(8,12)
Bcells <- subset(combined, idents = Bcells_clus)

Bcells <- NormalizeData(Bcells)
Bcells <- FindVariableFeatures(Bcells)
Bcells <- ScaleData(Bcells)
Bcells <- RunPCA(Bcells)
Bcells <- RunHarmony(Bcells,group.by.vars = "orig.ident")

Bcells <- RunUMAP(
  Bcells,
  reduction.key = "harmonyUMAP_",
  reduction = "harmony",
  reduction.name = "harmonyumap",
  dims = 1:30
)

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Bcells/"
dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/combined_samples.pdf"), width = 6, height = 5)
DimPlot(Bcells, reduction = "harmonyumap", label = TRUE, label.size = 4, group.by = "orig.ident")
dev.off()

pdf(paste0(savedir, "analysis/UMAP/combined_samples_split.pdf"), width = 6, height = 7)
DimPlot(Bcells, reduction = "harmonyumap", label = FALSE, label.size = 4, group.by = "orig.ident", split.by = "orig.ident", ncol = 4)
dev.off()

Bcells <- FindNeighbors(Bcells, reduction = "harmony", dims = 1:30)
Bcells <- FindClusters(Bcells, resolution = 0.4)

dir.create(paste0(savedir,"analysis/UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/UMAP/Bcells_cluster_res0.4.pdf"), width = 6, height = 5)
DimPlot(Bcells, reduction = "harmonyumap", group.by = "seurat_clusters", label = TRUE, label.size = 6)
dev.off()

genes <- c("CD27", "IGHD", "ITGAX", "FCRL5",  "IGHM", "IGHG", "CD38", "CR2")
dir.create(paste0(savedir,"analysis/featureplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/featureplots/Bcell.pdf"), width = 10)
FeaturePlot(Bcells, genes, reduction = "harmonyumap")
dev.off()

dir.create(paste0(savedir,"analysis/vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "analysis/vlnplots/Bcell_nopoint_res0.4.pdf"), width = 12, height = 10)
VlnPlot(Bcells, genes, pt.size = 0, group.by = "seurat_clusters")
dev.off()

Bcellmarkers <- FindAllMarkers(Bcells)
write.table(Bcellmarkers, "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Bcells/analysis/Bcellmarker_res0.2.txt", sep = "\t", quote = F, row.names = T, col.names = T)

dir.create(paste0(savedir, "saveRDS"), showWarnings = FALSE)
saveRDS(Bcells, paste0(savedir, "saveRDS/Bcells.RDS"))

#endregion
PBC <- readRDS("/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/CD4/saveRDS/CD4_obj.RDS")

DefaultAssay(PBC) <- "RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/.immune/Kyoto/resource/CD4_genes/", pattern = "txt", full.names = TRUE)
filename <- gsub(".txt", "", basename(files))
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[, 1]
    Tcellsubset <- rownames(PBC)[match(Tcellsubset, rownames(PBC), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

PBC <- AddModuleScore(PBC, features = markers, slot = "data")
colnames(PBC@meta.data)[29:45] <- paste0("CD4_", filename)

Idents(PBC) <- PBC@meta.data$seurat_clusters2

savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/liver/Tcells/remove_more_mps/CD4/"
dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir, "vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/CD4_pancancer_PBC_module_score_genes_violin.pdf"), width = 12, height = 12)
VlnPlot(PBC, paste0("CD4_", filename), pt.size = 0)
dev.off()

dir.create(paste0(savedir, "featureplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/CD4_pancancer_module_score_genes_featureplot.pdf"), width = 12, height = 12)
FeaturePlot(PBC, paste0("CD4_", filename), reduction = "harmonyumap", label = TRUE)
dev.off()

### CD8
DefaultAssay(PBC) <- "RNA"
rm(markers)
markers <- list()
files <- list.files("/diazlab/data3/.abhinav/.immune/Kyoto/resource/CD8_genes/", pattern = "txt", full.names = TRUE)
filename <- gsub(".txt", "", basename(files))
for (i in 1:length(files)) {
    Tcellsubset <- read.table(files[i])[, 1]
    Tcellsubset <- rownames(PBC)[match(Tcellsubset, rownames(PBC), nomatch = 0)]
    markers[[filename[i]]] <- Tcellsubset
    assign(filename[i], Tcellsubset)
}

PBC <- AddModuleScore(PBC, features = markers, slot = "data")
colnames(PBC@meta.data)[46:64] <- paste0("CD8_", filename)

# savedir <- "/diazlab/data3/.abhinav/.immune/Kyoto/PBC/"
dir.create(savedir, showWarnings = FALSE)
dir.create(paste0(savedir, "vlnplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/CD8_pancancer_module_score_genes_violin.pdf"), width = 12, height = 12)
VlnPlot(PBC, paste0("CD8_", filename), pt.size = 0)
dev.off()

dir.create(paste0(savedir, "featureplots/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "featureplots/CD8_pancancer_module_score_genes_featureplot.pdf"), width = 12, height = 12)
FeaturePlot(PBC, paste0("CD8_", filename), reduction = "harmonyumap", label = TRUE)
dev.off()

### making in the proper order
cluster_order = unique(PBC@meta.data$seurat_clusters2)[order(unique(PBC@meta.data$seurat_clusters2))]
new_order = c(0:5)

patterns <- cluster_order
replacements <- paste0("C",as.character(new_order))

library(stringr)
PBC@meta.data$seurat_clusters_new <- as.character(PBC@meta.data$seurat_clusters2)
for (i in seq_along(patterns)) {
    pattern <- paste0("\\b", patterns[i], "\\b")
    PBC@meta.data$seurat_clusters_new <- str_replace_all(PBC@meta.data$seurat_clusters_new, pattern, replacements[i])
}

pdf(paste0(savedir,"UMAP/PBC_new_cluster.pdf"), width = 4.5, height = 4)
DimPlot(PBC, group.by = "seurat_clusters_new", reduction = "harmonyumap", label = TRUE, 
label.size = 4)
dev.off()

#### Starting to run on google cloud Kytoto
# conda activate r_env
library(Seurat)
CD4_obj = readRDS("/mnt/data/projects/liver/Tcells/remove_more_mps/CD4/saveRDS/CD4_obj.RDS")

TFs= c("TCF7", "LEF1", "FOXP3", "TBX21", "MAF", "TOX", "PRDM1", "ZEB2", "EOMES", "BACH2")
Stem_T= c("IL7R", "CCR7", "GPR183")
Tph_markers= c("IL21", "ICOS", "BCL6", "CXCL13", "TNFSF8", "TOX2")
Cytokines= c("TNF", "IFNG", "CSF2", "LTB")
ICI= c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "CD96")
Cytotoxic= c("NKG7", "GZMB", "GZMK", "PRF1", "XCL1", "CD40LG")

gene_name = c("TFs", "Stem_T","Tph_markers","Cytokines","ICI","Cytotoxic")

savedir = "/mnt/data/projects/liver/Tcells/remove_more_mps/CD4/"
for(i in 1:length(gene_name)){
  pdf(paste0(savedir,"vlnplots/",gene_name[i],"_nodots.pdf"))
  print(VlnPlot(CD4_obj, get(gene_name[i]), group.by="seurat_clusters_new",  pt.size = 0))
  dev.off()
}

#region Tregs
library(Seurat)
library(dplyr)
library(harmony)
CD4_obj = readRDS("/mnt/data/projects/liver/Tcells/remove_more_mps/CD4/saveRDS/CD4_obj.RDS")
Tregs_cells <- rownames(CD4_obj@meta.data[grep("C4",CD4_obj@meta.data$seurat_clusters_new),])
Treg_obj = subset(CD4_obj, cells = Tregs_cells)

Treg_obj <- NormalizeData(Treg_obj)
Treg_obj <- FindVariableFeatures(Treg_obj)
Treg_obj <- ScaleData(Treg_obj)
Treg_obj <- RunPCA(Treg_obj)
Treg_obj <- RunHarmony(Treg_obj,group.by.vars = "orig.ident")

Treg_obj <- RunUMAP(
  Treg_obj,
  reduction.key = "harmonyUMAP_",
  reduction = "harmony",
  reduction.name = "harmonyumap",
  dims = 1:30
)

Treg_obj <- FindNeighbors(Treg_obj, reduction = "harmony", dims = 1:30)

library(Seurat)
savedir <- "/mnt/data/projects/liver/Tcells/remove_more_mps/CD4/Tregs/"
Treg_obj <- readRDS(paste0(savedir,"saveRDS/Treg.RDS"))
Treg_obj <- FindClusters(Treg_obj, resolution = 0.15)

Idents(Treg_obj) = Treg_obj@meta.data$RNA_snn_res.0.2
Tregs_cellnames = rownames(Treg_obj@meta.data[grep("7|8",Treg_obj@meta.data$RNA_snn_res.0.5,invert=T),])

dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "UMAP/Tregs_samples.pdf"), width = 6, height = 5)
DimPlot(Treg_obj, reduction = "harmonyumap", label = TRUE, label.size = 4, group.by =  "orig.ident")
dev.off()

pdf(paste0(savedir, "UMAP/Treg_obj_cluster_0.15.pdf"), width = 6, height = 5)
DimPlot(Treg_obj, reduction = "harmonyumap", group.by = "seurat_clusters", label = TRUE, label.size = 4)
dev.off()

Treg_obj <- JoinLayers(Treg_obj)
Treg_obj<- NormalizeData(Treg_obj)
Treg_markers <- FindAllMarkers(Treg_obj)
dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(Treg_markers, paste0(savedir,"Table/Treg_marker_res0.15.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

dir.create(paste0(savedir,"saveRDS"))
saveRDS(Treg_obj, paste0(savedir,"saveRDS/Treg.RDS"))

#### Longitudinal analysis
celltypes_ind <- table(Treg_obj@meta.data$orig.ident,Treg_obj@meta.data$seurat_clusters)
dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(celltypes_ind, paste0(savedir,"Table/celltype_individual_res_1.2.txt"),
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

df_melted$condition <- gsub("_.*.","",df_melted$sample)

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=sample)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table/marker_celltype_res_0.5_pts_timepoints.pdf"),width =10, height = 7)
p
dev.off()

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=condition)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table//marker_celltype_res_0.5_condition.pdf"),width =10, height = 7)
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
  stat.test,  label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

pdf(paste(savedir,"Table/t_test_unpaired_all_sample_0.5.pdf",sep = ""), width = 8, height = 6)
bxp4
dev.off()

dir.create(paste0(savedir,"saveRDS"))
saveRDS(Treg_obj, paste0(savedir,"saveRDS/Treg.RDS"))

#endregion

#region Treg removing cluster 7 and 8
library(Seurat)
savedir <- "/mnt/data/projects/liver/Tcells/remove_more_mps/CD4/Tregs/"
Treg_obj <- readRDS(paste0(savedir,"saveRDS/Treg.RDS"))
Treg_obj <- FindClusters(Treg_obj, resolution = 0.15)

Idents(Treg_obj) = Treg_obj@meta.data$RNA_snn_res.0.5
Tregs_cellnames = rownames(Treg_obj@meta.data[grep("7|8",Treg_obj@meta.data$RNA_snn_res.0.5,invert=T),])
Treg_obj2 = subset(Treg_obj, cells = Tregs_cellnames)

Treg_obj2@meta.data$seurat_clusters <- Treg_obj2@meta.data$RNA_snn_res.0.5

pdf(paste0(savedir,"UMAP/umap_removed_7_8.pdf"))
DimPlot(Treg_obj2, label = T, label.size = 8)
dev.off()

pdf(paste0(savedir,"UMAP/umap_removed_7_8_splitted.pdf"), width = 15, height = 8)
DimPlot(Treg_obj2, split.by = "orig.ident", label = T, label.size = 8, ncol = 5)
dev.off()

pdf(paste0(savedir,"UMAP/umap_removed_7_8_splitted.pdf"), width = 15, height = 8)
DimPlot(Treg_obj2, split.by = "orig.ident", label = T, label.size = 8, ncol = 5)
dev.off()

Treg_obj2@meta.data$condition = gsub("_.*.","",Treg_obj2@meta.data$orig.ident)

pdf(paste0(savedir,"UMAP/umap_removed_7_8_splitted_condition.pdf"), width = 9, height = 4.5)
DimPlot(Treg_obj2, split.by = "condition", label = T, label.size = 8)
dev.off()

genes = c("FOXP3", "IKZF2", "TIGIT", "CTLA4", "ICOS", "HLA-DRB1", "TNFRSF9",
"TCF7", "LEF1", "CXCR6", "CXCR3", "IL10", "TBX21", "RORC",
"BCL6", "PDCD1", "CXCR5", "IL21R", "CCR7", "IL2RA", "FCRL3" )

dir.create(paste0(savedir,"vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/gene_vlnplots.pdf"), width = 12, height = 15)
VlnPlot(Treg_obj2, genes, pt.size = 0)
dev.off()

#### Longitudinal analysis
celltypes_ind <- table(Treg_obj2@meta.data$orig.ident,Treg_obj2@meta.data$seurat_clusters)
dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(celltypes_ind, paste0(savedir,"Table/celltype_individual_res_0.5_remove7_8.txt"),
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

df_melted$condition <- gsub("_.*.","",df_melted$sample)

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=sample)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table/marker_celltype_res_0.5_pts_timepoints_removed_7_8.pdf"),width =10, height = 7)
p
dev.off()

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=condition)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table//marker_celltype_res_0.5_condition_removed_7_8.pdf"),width =10, height = 7)
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
  stat.test,  label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

pdf(paste(savedir,"Table/t_test_unpaired_all_sample_0.5_removed_7_8.pdf",sep = ""), width = 8, height = 6)
bxp4
dev.off()

dir.create(paste0(savedir,"saveRDS"))
saveRDS(Treg_obj2, paste0(savedir,"saveRDS/Treg2.RDS"))

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
DefaultAssay(Treg_obj2) <- "RNA"
sce <- as.SingleCellExperiment(Treg_obj2, assay = "RNA")

## Starting Cluster
savedir <- paste0(savedir,"trajectory/")
dir.create(savedir, showWarnings = FALSE)
dir.create(paste(savedir, "UMAP/", sep = ""), showWarnings = FALSE)
p <- DimPlot(Treg_obj2, reduction = "harmonyumap", label = TRUE, group.by = "seurat_clusters")
pdf(paste(savedir, "UMAP/Tregs.pdf", sep = ""))
p
dev.off()

# CD4_sce <- slingshot(sce,
#   clusterLabels = "seurat_clusters",
#   reducedDim = "HARMONYUMAP", approx_points = 100,
#   omega = TRUE,
#   omega_scale = 1.5
# )

CD4_sce_naive <- slingshot(sce,
  clusterLabels = colData(sce)$seurat_clusters,
  reducedDim = "HARMONYUMAP",
  start.clus = c("1"),
  # end.clus = c("7"),
  approx_points = 100,
  omega = TRUE,
  omega_scale=1.5
)

library(ArchR)
library(scater)
embedded_orig <- embedCurves(CD4_sce_naive, "HARMONYUMAP")

rm(plot_list)
plot_list <- list()
for (i in 1:4) {
  embedded <- slingCurves(embedded_orig)[[i]] # only 1 path.
  embedded <- data.frame(embedded$s[embedded$ord, ])
  g <- plotUMAP(CD4_sce_naive, colour_by = paste("slingPseudotime_", i, sep = ""), dimred = "HARMONYUMAP")
  stopifnot(all(rownames(Treg_obj2@reductions$umap@cell.embeddings) == rownames(g$data)))
  # data <- merge(embedded_orig@reductions$umap@cell.embeddings, g$data, by = "row.names")
  # colnames(data) <- c("cellname", "HARMONYUMAP_1", "HARMONYUMAP_2", "X", "Y", paste("Lineage", i, sep = ""))
  colnames(g$data) <- c("HARMONYUMAP_1","HARMONYUMAP_2",paste("Lineage", i, sep = ""))
  p <- ggplot(g$data, aes_string("HARMONYUMAP_1", "HARMONYUMAP_2", color = paste("Lineage", i, sep = ""))) +
    geom_point(size = 1) +
    scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
    geom_path(data = embedded, aes(x = harmonyUMAP_1, y = harmonyUMAP_2), color = "black", size = 1.2) +
    theme_bw()
  plot_list[[i]] <- p
}

require(gridExtra)
pdf(paste(savedir, "UMAP/Liver_CD4_clus_1_lineages.pdf", sep = ""), width = 12, height = 10)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]],plot_list[[4]], nrow = 2, ncol = 2)
dev.off()

pdf(paste(savedir, "UMAP/Liver_Treg_UMAP_lineage_clus1.pdf", sep = ""), width = 4.5, height = 4)
plot_list
dev.off()

#region TPH cell extract
library(Seurat)
library(dplyr)
library(harmony)

CD4 = readRDS("/mnt/data/projects/liver/Tcells/remove_more_mps/CD4/saveRDS/CD4_obj.RDS")

Tphs_cells <- rownames(CD4@meta.data[grep("C1",CD4@meta.data$seurat_clusters_new),])
Tph_obj = subset(CD4, cells = Tphs_cells)

Tph_obj <- NormalizeData(Tph_obj)
Tph_obj <- FindVariableFeatures(Tph_obj)
Tph_obj <- ScaleData(Tph_obj)
Tph_obj <- RunPCA(Tph_obj)
Tph_obj <- RunHarmony(Tph_obj,group.by.vars = "orig.ident")

Tph_obj <- RunUMAP(
  Tph_obj,
  reduction.key = "harmonyUMAP_",
  reduction = "harmony",
  reduction.name = "harmonyumap",
  dims = 1:30
)

Tph_obj <- FindNeighbors(Tph_obj, reduction = "harmony", dims = 1:30)

library(Seurat)
savedir <- "/mnt/data/projects/liver/Tcells/remove_more_mps/CD4/Tphs/"
dir.create(savedir,showWarnings = FALSE)
setwd(savedir)
Tph_obj <- FindClusters(Tph_obj, resolution = 0.3)

dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "UMAP/Tphs_samples.pdf"), width = 6, height = 5)
DimPlot(Tph_obj, reduction = "harmonyumap", label = TRUE, label.size = 4, group.by =  "orig.ident")
dev.off()

pdf(paste0(savedir, "UMAP/Tph_obj_cluster.pdf"), width = 6, height = 5)
DimPlot(Tph_obj, reduction = "harmonyumap", group.by = "seurat_clusters", label = TRUE, label.size = 4)
dev.off()

Tph_obj <- JoinLayers(Tph_obj)
Tph_obj<- NormalizeData(Tph_obj)
Tph_markers <- FindAllMarkers(Tph_obj)
dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(Tph_markers, paste0(savedir,"Table/Tph_markers.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

#### Longitudinal analysis
celltypes_ind <- table(Tph_obj@meta.data$orig.ident,Tph_obj@meta.data$seurat_clusters)
dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(celltypes_ind, paste0(savedir,"Table/celltype_individual.txt"),
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

df_melted$condition <- gsub("_.*.","",df_melted$sample)

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=sample)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table/marker_celltype_res_pts_timepoints.pdf"),width =10, height = 7)
p
dev.off()

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=condition)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table//marker_celltype_res_condition.pdf"),width =10, height = 7)
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
  stat.test,  label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

pdf(paste(savedir,"Table/t_test_unpaired_all_sample_0.5.pdf",sep = ""), width = 8, height = 6)
bxp4
dev.off()

dir.create(paste0(savedir,"saveRDS"))
saveRDS(Tph_obj, paste0(savedir,"saveRDS/Tph.RDS"))

genes = c("TCF7", "LEF1", "IFNG", "CXCR5", 'CXCR6',"CCR7", "IL21", "TNFSF8", "CCL5", "GZMK", "GZMA", "RELB")

dir.create(paste0(savedir,"vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/Tph_genes_vlnplots.pdf"), width = 12, height = 10)
VlnPlot(Tph_obj, genes, pt.size = 0)
dev.off()

#endregion

#region Tph extract
library(Seurat)
library(dplyr)
library(harmony)

library(Seurat)
savedir = "/mnt/data/projects/liver/Tcells/remove_more_mps/CD4/Tphs/"
Tph_obj = readRDS(paste0(savedir,"saveRDS/Tph.RDS"))

Tphs_cells <- rownames(Tph_obj@meta.data[grep("^5$",Tph_obj@meta.data$seurat_clusters, invert = TRUE),])
Tph_obj2 = subset(Tph_obj, cells = Tphs_cells)

Tph_obj2 <- NormalizeData(Tph_obj2)
Tph_obj2 <- FindVariableFeatures(Tph_obj2)
Tph_obj2 <- ScaleData(Tph_obj2)
Tph_obj2 <- RunPCA(Tph_obj2)
Tph_obj2 <- RunHarmony(Tph_obj2,group.by.vars = "orig.ident")

Tph_obj2 <- RunUMAP(
  Tph_obj2,
  reduction.key = "harmonyUMAP_",
  reduction = "harmony",
  reduction.name = "harmonyumap",
  dims = 1:30
)

Tph_obj2 <- FindNeighbors(Tph_obj2, reduction = "harmony", dims = 1:30)

library(Seurat)
savedir <- "/mnt/data/projects/liver/Tcells/remove_more_mps/CD4/Tphs/Tphs_remove_clus5/"
dir.create(savedir,showWarnings = FALSE)
setwd(savedir)
Tph_obj2 <- FindClusters(Tph_obj2, resolution = 0.15)

dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "UMAP/Tphs_samples.pdf"), width = 6, height = 5)
DimPlot(Tph_obj2, reduction = "harmonyumap", label = TRUE, label.size = 4, group.by =  "orig.ident")
dev.off()

pdf(paste0(savedir, "UMAP/Tph_obj2_cluster.pdf"), width = 6, height = 5)
DimPlot(Tph_obj2, reduction = "harmonyumap", group.by = "seurat_clusters", label = TRUE, label.size = 4)
dev.off()

Tph_obj2 <- JoinLayers(Tph_obj2)
Tph_obj2<- NormalizeData(Tph_obj2)
Tph_markers <- FindAllMarkers(Tph_obj2)
dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(Tph_markers, paste0(savedir,"Table/Tph_markers.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

#### Longitudinal analysis
celltypes_ind <- table(Tph_obj2@meta.data$orig.ident,Tph_obj2@meta.data$seurat_clusters)
dir.create(paste0(savedir,"Table"), showWarnings = FALSE)
write.table(celltypes_ind, paste0(savedir,"Table/celltype_individual.txt"),
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

df_melted$condition <- gsub("_.*.","",df_melted$sample)

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=sample)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table/marker_celltype_res_pts_timepoints.pdf"),width =10, height = 7)
p
dev.off()

p <- ggplot(df_melted, aes(fill=celltype, y=percentage, x=condition)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c('#1a69a1','#436957','#5fe3a7','#107044',"#6fb9ed",'#aa40fc',
                               '#8c564b','black','#b5bd61','hotpink','#17becf',
                               "brown3","cornsilk3","yellow","brown3")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf(paste0(savedir,"Table//marker_celltype_res_condition.pdf"),width =10, height = 7)
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
  stat.test,  label = "p", tip.length = 0
) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

pdf(paste(savedir,"Table/t_test_unpaired_all_sample_0.5.pdf",sep = ""), width = 8, height = 6)
bxp4
dev.off()

dir.create(paste0(savedir,"saveRDS"))
saveRDS(Tph_obj2, paste0(savedir,"saveRDS/Tph_remove_clus5.RDS"))

genes = c("TCF7", "LEF1", "IFNG", "CXCR5", 'CXCR6',"CCR7", "IL21", "TNFSF8", "CCL5", "GZMK", "GZMA", "RELB")

dir.create(paste0(savedir,"vlnplots"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(savedir, "vlnplots/Tph_genes_vlnplots.pdf"), width = 12, height = 10)
VlnPlot(Tph_obj2, genes, pt.size = 0)
dev.off()

### TpH cells

