# region preprocessing
# Visium HD support in Seurat
# We note that Visium HD data is generated from spatially patterned olignocleotides labeled in 2um x 2um bins. However, since the data from this resolution is sparse,
# adjacent bins are pooled together to create 8um and 16um resolutions. 10x recommends the use of 8um binned data for analysis, but Seurat supports in the simultaneous
# loading of multiple binnings - and stores them in a single object as multiple assays.

# Unsupervised clustering
# Identification of spatial tissue domains
# Subsetting spatial regions
# Integration with scRNA-seq data
# Comparing the spatial localization of different cell types

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(arrow)

localdir <- "/diazlab/data3/.abhinav/.immune/Kyoto/liver/visium_HD/rawdata/"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))

# Setting default assay changes between 8um and 16um binning
Assays(object)
DefaultAssay(object) <- "Spatial.008um"

matrix_dir <- paste0(localdir, "GSE289717_square_008um_filtered_feature_bc_matrix")
spatial_dir <- paste0(localdir, "GSE289717_square_008um_spatial")

expression_matrix <- Read10X(data.dir = matrix_dir)
seurat_obj <- CreateSeuratObject(counts = expression_matrix)

seurat_obj <- Load10X_Spatial(data.dir = spatial_dir, assay = "Spatial", slice = "slice1", seurat.object = seurat_obj, load.image = FALSE)

# Read spatial metadata from Visium HD
positions <- read_parquet(file.path(spatial_dir, "tissue_positions.parquet"))
positions <- positions %>% column_to_rownames("barcode")
position_rq <- positions[match(rownames(seurat_obj@meta.data), rownames(positions), nomatch = 0), ]

stopifnot(all(rownames(position_rq) == rownames(seurat_obj@meta.data)))

seurat_obj@meta.data$pxl_row_in_fullres <- position_rq$pxl_row_in_fullres
seurat_obj@meta.data$pxl_col_in_fullres <- position_rq$pxl_col_in_fullres

# Add spatial coordinates (x = col, y = row)
spatial_coords <- seurat_obj@meta.data %>%
    select(pxl_col_in_fullres, pxl_row_in_fullres) %>%
    as.matrix()

colnames(spatial_coords) <- c("spatial_1", "spatial_2")

# Attach spatial coordinates
seurat_obj[["spatail_coord"]] <- Seurat::CreateDimReducObject(
    embeddings = spatial_coords,
    key = "spatial_",
    assay = DefaultAssay(seurat_obj)
)

# Unsupervised clustering
# While the standard scRNA-seq clustering workflow can also be applied to spatial datasets - we have observed that when working with Visium HD datasets, the Seurat
# v5 sketch clustering workflow exhibits improved performance, especially for identifying rare and spatially restricted groups.
# As described in Hao et al, Nature Biotechnology 2023 and Hie et al, sketch-based analyses aim to ‘subsample’ large datasets in a way that preserves rare populations.
# Here, we sketch the Visium HD dataset, perform clustering on the subsampled cells, and then project the cluster labels back to the full dataset.
# Details of the sketching procedure and workflow are described in Hao et al, Nature Biotechnology 2023 and the Seurat v5 sketch clustering vignette. Since the full
# Visium HD dataset fits in memory, we do not use any of the on-disk capabilities of Seurat v5 in this vignette.

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
# we select 50,0000 cells and create a new 'sketch' assay
seurat_obj <- SketchData(
    object = seurat_obj,
    ncells = 50000,
    method = "LeverageScore",
    sketched.assay = "sketch"
)

# switch analysis to sketched cells
DefaultAssay(seurat_obj) <- "sketch"

# perform clustering workflow
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, assay = "sketch", reduction.name = "pca.sketch")
seurat_obj <- FindNeighbors(seurat_obj, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
seurat_obj <- FindClusters(seurat_obj, cluster.name = "seurat_cluster.sketched", resolution = 3)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)


# Now we can project the cluster labels, and dimensional reductions (PCA and UMAP) that we learned from the 50,000 sketched cells - to the entire dataset, using the ProjectData function.
# In the resulting object, for all cells:
#     cluster labels will be stored in object$seurat_cluster.projected
#     Projected PCA embeddings will be stored in object[["pca.008um"]]
#     Projected UMAP embeddings will be stored in object[["umap.sketch"]]

seurat_obj <- ProjectData(
    object = seurat_obj,
    assay = "RNA",
    full.reduction = "full.pca.sketch",
    sketched.assay = "sketch",
    sketched.reduction = "pca.sketch",
    umap.model = "umap.sketch",
    dims = 1:50,
    refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

# endregion
