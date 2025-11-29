# Kyoto Liver Immunology Research Repository

## Overview

This repository contains comparative single-cell transcriptomics analyses of **liver autoimmune disease (PBC) and cancer (HCC)** in relation to immunotherapy response. The project integrates multimodal scCITE-Seq data from multiple studies to identify immune cell subsets, phenotypic markers, and therapeutic biomarkers.

**Research Focus:**

- **PBC (Primary Biliary Cholangitis)**: Autoimmune liver disease
- **HCC (Hepatocellular Carcinoma)**: Liver cancer with immunotherapy response stratification
- **Immune Cell Analysis**: T cells (CD4, CD8, Tregs), B cells, and their phenotypic heterogeneity

---

## Repository Structure

```
Kyoto/
├── README.md                          # This file
├── Readme.md                          # Original minimal README
├── .gitignore                         # Git configuration
│
├── pipeline_functions/                # Shared reusable analysis functions
│   ├── CR7_scCITESeq_QC.R            # Core QC & filtering
│   ├── cellcycle_mitoscore.R         # Cell cycle & mito scoring
│   ├── scCITEseq_Doublet_Finder.R    # DoubletFinder wrapper
│   ├── NN_clustertree.R              # Nearest neighbor clustering
│   ├── pseudobulk_within_cluster_AJ_no_batch_correction.R  # Bulk RNA aggregation
│   ├── RNAseq_limma_EdgeR.R          # Differential expression (bulk)
│   ├── genescore_custom.R            # Custom gene scoring
│   ├── Volcano_plot.R                # DE visualization
│   ├── MA_plot.R                     # MA plot visualization
│   └── scCITE_QC_FIbro.r             # Fibroblast-specific QC
│
├── NatMed_HCC/                        # Nature Medicine HCC Study
│   ├── HCC.R                          # Main HCC CD4+ T cell analysis
│   ├── HCC_CD8.R                      # HCC CD8+ T cell analysis
│   ├── CD8/
│   │   └── HCC_CD8.R                  # CD8 subset processing
│   ├── CD8_2/                         # Alternative CD8 clustering
│   ├── CD8_3/                         # Further CD8 refinement
│   ├── analysis/
│   │   ├── UMAP/                      # Dimension reduction plots
│   │   ├── Tregs/                     # Treg-specific outputs
│   │   │   └── removed_clus12_14/     # Outlier-removed Treg analysis
│   │   └── [other analysis outputs]
│   └── rawdata/
│       ├── GSE206325_data_SingleCellExperiment_object.Rda
│       └── GSE206325_full_HCC_cluster_annotation.csv
│
├── PBC_HCC/                           # PBC vs HCC Comparative Analysis
│   ├── PBC_HCC_Tregs.r               # Cross-disease Treg integration
│   ├── Tfh/                          # T follicular helper cells
│   ├── Tregs/                        # Regulatory T cells
│   ├── Trajectory/                   # Pseudotime/trajectory analysis
│   ├── stemlike/                     # Stem-like cell identification
│   ├── featureplots/                 # Gene expression overlays
│   ├── vlnplots/                     # Violin plot outputs
│   ├── UMAP/                         # Dimension reduction plots
│   ├── TFs/                          # Transcription factor analysis
│   ├── Table/                        # Statistical tables & summaries
│   └── saveRDS/                      # Processed Seurat objects
│
├── PBC_HCC_CD8/                       # PBC vs HCC CD8+ T Cell Analysis
│   ├── PBC_HCC_CD8.R                 # CD8 integration & analysis
│   ├── UMAP/                         # CD8-specific UMAPs
│   ├── vlnplots/                     # CD8 violin plots
│   ├── Table/                        # CD8 statistics
│   └── saveRDS/                      # CD8 Seurat objects
│
└── liver/                             # Legacy liver analysis files
    ├── Tcells/                        # T cell subset analysis
    │   ├── remove_more_mps/           # Immune subset refinement
    │   │   ├── CD4/
    │   │   │   └── Tregs/
    │   │   │       └── saveRDS/       # Processed Treg objects
    │   │   └── CD8/
    │   └── ...
    └── [other liver-specific analyses]
```

---

## Project Descriptions

### 1. NatMed_HCC (Nature Medicine Study)

**Objective:** Characterize CD4+ and CD8+ T cell subsets in HCC and identify therapy-responsive phenotypes

**Key Analysis Steps:**

1. Load GSE206325 HCC scRNA-seq data (SingleCellExperiment format)
2. Subset to CD4+, Naive, Tregs, Proliferating cells
3. QC filtering and doublet removal
4. Integration using Harmony (cross-sample batch correction)
5. Dimensionality reduction (UMAP)
6. Cell type annotation using provided metadata

**Main Scripts:**

- `HCC.R` - CD4+ T cell analysis
- `HCC_CD8.R` - CD8+ T cell analysis

**Output:** Annotated Seurat objects with UMAP coordinates and cluster assignments

---

### 2. PBC_HCC (Comparative PBC vs HCC Analysis)

**Objective:** Compare immune cell subsets between autoimmune liver disease (PBC) and cancer (HCC)

**Key Features:**

- Cross-disease integration (PBC vs HCC)
- Focused analysis on functionally relevant subsets:
  - **Tregs**: Regulatory T cells (immune tolerance)
  - **Tfh**: T follicular helper cells (adaptive immunity)
  - **CD8 T cells**: Cytotoxic immunity
- Pseudotime analysis: Track immune differentiation trajectories
- Transcription factor analysis: Identify regulatory networks

**Main Scripts:**

- `PBC_HCC_Tregs.r` - Integrates Treg objects from PBC and HCC
  - Handles Seurat v4↔v5 compatibility
  - Aligns gene expression data
  - Adds disease labels for downstream comparison

**Workflow:**

```r
# Load pre-processed Treg objects
PBC_obj <- readRDS("...")  # Seurat v5
HCC_obj <- readRDS("...")  # Seurat v4

# Harmonize object structure
HCC_obj <- DietSeurat(HCC_obj, assays="RNA", layers="counts")

# Combine gene expression matrices
combined_counts <- cbind(counts_HCC, counts_PBC)

# Integrated analysis
```

---

### 3. PBC_HCC_CD8 (CD8+ T Cell Focused Study)

**Objective:** Deep characterization of CD8+ T cells in PBC and HCC

**Analysis Types:**

- CD8+ subset identification (naive, effector, exhausted, memory)
- Phenotypic profiling (TCR usage, exhaustion markers)
- Cross-disease comparison of anti-tumor vs anti-self immunity

**Key Outputs:**

- UMAP plots colored by subset
- Violin plots of marker genes
- Comparative statistics tables

---

## Core Analysis Workflow

### Typical Analysis Pipeline

```
Raw Data (CellRanger/GEO)
    ↓
QC Filtering [CR7_scCITESeq_QC.R]
    ↓
Doublet Removal [scCITEseq_Doublet_Finder.R]
    ↓
Cell Cycle & Mito Scoring [cellcycle_mitoscore.R]
    ↓
Integration (Harmony/SCTransform)
    ↓
Dimensionality Reduction (PCA → UMAP)
    ↓
Clustering [NN_clustertree.R]
    ↓
Pseudo-bulk Aggregation [pseudobulk_within_cluster_AJ_no_batch_correction.R]
    ↓
Differential Expression [RNAseq_limma_EdgeR.R]
    ↓
Visualization [Volcano_plot.R, MA_plot.R, Feature plots]
```

---

## Key Functions Reference

### QC & Preprocessing

| Function              | File                         | Purpose                            |
| --------------------- | ---------------------------- | ---------------------------------- |
| `scCITE_QC()`         | `CR7_scCITESeq_QC.R`         | Core filtering (genes, UMI, mito%) |
| `doublet_scCITEseq()` | `scCITEseq_Doublet_Finder.R` | DoubletFinder workflow             |
| `cellcycle_mito()`    | `cellcycle_mitoscore.R`      | Cell cycle & mito effect scoring   |
| `scCITE_QC_Fibro()`   | `scCITE_QC_FIbro.r`          | Fibroblast-specific QC             |

### Analysis & Visualization

| Function             | File                                                 | Purpose                     |
| -------------------- | ---------------------------------------------------- | --------------------------- |
| `nearest_neigbour()` | `NN_clustertree.R`                                   | Clustering with clustertree |
| `pseudobulk_*()`     | `pseudobulk_within_cluster_AJ_no_batch_correction.R` | Aggregate to pseudo-bulk    |
| `limma_EdgeR_*()`    | `RNAseq_limma_EdgeR.R`                               | Differential expression     |
| `volcano_plot()`     | `Volcano_plot.R`                                     | DE visualization            |
| `MA_plot()`          | `MA_plot.R`                                          | MA plot for DE results      |

---

## Data Format & Requirements

### Input Data Formats

1. **SingleCellExperiment (SCE)** → Convert to Seurat

   ```r
   seurat_obj <- as.Seurat(sce, counts = "counts", data = NULL)
   ```

2. **Seurat v4 Objects** (old format)

   - `.RDS` files with RNA assay
   - May have legacy layer names

3. **Seurat v5 Objects** (new format)
   - Layer structure: `counts`, `data`, `scale.data`
   - More efficient memory usage

### Compatibility Note

This repository handles **v4 ↔ v5 conversion**:

```r
# Load v4, keep only counts layer
obj <- DietSeurat(obj, assays = "RNA", layers = "counts")

# Remove extra layers
layers_to_remove <- setdiff(Layers(obj[["RNA"]]), "counts")
for (layer in layers_to_remove) {
    obj[["RNA"]]@layers[[layer]] <- NULL
}
```

---

## Key Analyses Explained

### 1. Harmony Integration

Corrects for batch effects across samples while preserving biological variation:

```r
library(harmony)
seurat_obj <- RunHarmony(seurat_obj,
                         group.by.vars = "sample_ID",
                         reduction = "pca")
```

### 2. Pseudotime Analysis

Trajectory inference to model immune cell differentiation:

- Orders cells along pseudotime axis
- Identifies transitional states
- Used in `Trajectory/` subdirectory

### 3. Pseudobulk Analysis

Aggregates single-cell data to sample-level counts for bulk RNA-seq methods:

- Improves statistical power
- Enables meta-analysis integration
- Removes cell-type-specific batch effects

### 4. Transcription Factor (TF) Analysis

Identifies regulatory networks controlling immune phenotypes:

- Stored in `TFs/` subdirectory
- TF activity inference (regulon analysis)
- Prioritizes master regulators

---

## Running Analyses

### Example: Run PBC vs HCC Treg Integration

```bash
cd /mnt/data/projects/Kyoto
Rscript PBC_HCC/PBC_HCC_Tregs.r
```

**Expected Output:**

- Merged Seurat object with PBC + HCC Tregs
- Metadata annotations (disease, sample_ID, etc.)
- Saved RDS file in `saveRDS/`

### Example: Differential Expression

```r
library(limma)
source("pipeline_functions/RNAseq_limma_EdgeR.R")

# Run DE on pseudobulk data
results <- limma_EdgeR_analysis(pseudo_bulk_counts,
                                design_matrix,
                                contrast)
```

---

## Data Management & Storage

### File Organization

- **Raw data**: `rawdata/` (GEO downloads, original experiments)
- **Processed objects**: `saveRDS/` (Seurat RDS files after QC)
- **Plots**: `UMAP/`, `vlnplots/`, `featureplots/` (publication-ready figures)
- **Tables**: `Table/` (Statistics, DE results, cell counts)

### Large File Handling

- RDS files can be >1GB (Seurat objects with normalized + scaled data)
- Recommended: Use `DietSeurat()` to reduce size before storing
- Archive analyzed data, version control only code

---

## Dependencies & Installation

### R Packages Required

```r
# Core analysis
install.packages(c("Seurat", "dplyr", "ggplot2"))

# Integration & normalization
install.packages("harmony")
install.packages("sctransform")

# Differential expression
install.packages("limma")
install.packages("edgeR")

# Specialty tools
install.packages("DoubletFinder")      # Doublet detection
install.packages("clustree")           # Clustering resolution
install.packages("monocle3")           # Trajectory (if used)
```

**Check installation:**

```r
library(Seurat)
library(harmony)
packageVersion("Seurat")
```

---

## Common Analyses & Workflows

### 1. Subset & Re-cluster a Cell Type

```r
# Load integrated object
obj <- readRDS("saveRDS/combined.RDS")

# Subset to CD4+ Tregs
treg_cells <- colnames(obj)[obj@meta.data$celltype == "Treg"]
tregs <- subset(obj, cells = treg_cells)

# Re-normalize & cluster
tregs <- NormalizeData(tregs)
tregs <- FindVariableFeatures(tregs)
tregs <- ScaleData(tregs)
tregs <- RunPCA(tregs)
tregs <- RunUMAP(tregs, dims = 1:30)
tregs <- FindClusters(tregs, resolution = 0.5)
```

### 2. Compare Gene Expression Between Diseases

```r
# Subset PBC and HCC cells separately
pbc_cells <- colnames(obj)[obj@meta.data$disease == "PBC"]
hcc_cells <- colnames(obj)[obj@meta.data$disease == "HCC"]

# Marker genes for HCC vs PBC
markers_hcc <- FindMarkers(obj,
                           ident.1 = hcc_cells,
                           ident.2 = pbc_cells)
```

### 3. Generate Publication Figures

```r
# UMAP with custom colors
p1 <- DimPlot(obj, group.by = "celltype",
              cols = my_color_palette,
              reduction = "umap")

# Feature plot for marker gene
p2 <- FeaturePlot(obj, features = "FOXP3",
                  reduction = "umap")

# Violin plot for gene expression by group
p3 <- VlnPlot(obj, features = "FOXP3",
              group.by = "disease",
              split.by = "tissue")
```

---

## Troubleshooting

| Issue                                 | Solution                                             |
| ------------------------------------- | ---------------------------------------------------- |
| "Harmony convergence error"           | Increase `max.iter`, reduce `nclust` in Harmony      |
| "Memory overflow with large objects"  | Use `DietSeurat()` to remove unused layers           |
| "v4 vs v5 compatibility error"        | Re-save object with `saveRDS()` after `DietSeurat()` |
| "DE results all NA"                   | Check pseudobulk aggregation; verify design matrix   |
| "UMAP looks noisy"                    | Increase `n.neighbors`, adjust `min.dist` parameters |
| "Missing marker genes in annotations" | Expand marker gene list, check for aliases           |

---

## Documentation & Resources

- **Seurat Documentation**: https://satijalab.org/seurat/
- **Harmony**: https://github.com/immunogenomics/harmony
- **DoubletFinder**: https://github.com/chris-mcginnis-ucsf/DoubletFinder
- **Limma**: https://bioconductor.org/packages/limma/
- **GEO Database**: https://www.ncbi.nlm.nih.gov/geo/ (raw data source)

---

## Project Timeline

- **NatMed_HCC**: Main HCC immunotherapy study
- **PBC_HCC**: Comparative autoimmune vs cancer analysis
- **PBC_HCC_CD8**: Focused CD8+ characterization
- **liver/**: Legacy analyses and intermediate results

---

## Version Control & Collaboration

This repository uses Git for version control. Latest changes tracked in `.git/`

```bash
# Clone repository
git clone git@github.com:Ajaingithub/Kyoto.git

# View commit history
cd Kyoto
git log --oneline

# Check recent changes
git status
```

---

## Citation & Acknowledgments

If you use code or data from this repository, please cite:

- **Primary study**: Nature Medicine HCC (GSE206325)
- **Repository**: https://github.com/Ajaingithub/Kyoto

---

## Contact

**Repository Owner**: Abhinav Jain (abhinavjj@gmail.com)

For questions about specific analyses, refer to script headers or README in individual project folders.

---

**Last Updated**: November 29, 2025
