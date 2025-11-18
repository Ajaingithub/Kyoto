### Running it on Slurm
#### After aggregation Performing QC and combining and integration
### Loading the data and performing the QC on the object
library(Seurat)
samples <- c(
    "control_liver1_aggregated", "control_liver2_aggregated", "control_liver3_aggregated", "control_liver4_aggregated", "control_liver5_aggregated",
    "pbc_liver1_aggregated", "pbc_liver2_aggregated", "pbc_liver3_aggregated", "pbc_liver4_aggregated", "pbc_liver5_aggregated"
)

source("/diazlab/data3/.abhinav/resources/all_scripts/R/scRNA_QC2.R")

for (i in 1:length(samples)) {
    sample_obj <- scRNA_QC(Dir = "/diazlab/data3/.abhinav/.immune/Kyoto/liver/", Sample = samples[i], saveDir = "/diazlab/data3/.abhinav/.immune/Kyoto/liver/analysis/")
    assign(paste0(samples[i], "_obj"), sample_obj)
}

objname <- ls(pattern = "_aggregated_obj")

## We have to combine the RNA singlet object to perform the normalization for the unwanted the effects.
combined <- merge(
    x = get(objname[1]), y = c(
        get(objname[2]), get(objname[3]), get(objname[4]), get(objname[5]), get(objname[6]),
        get(objname[7]), get(objname[8]), get(objname[9]), get(objname[10])
    ),
    add.cell.ids = objname, project = "combined"
)

saveRDS(combined, "/diazlab/data3/.abhinav/.immune/Kyoto/liver/analysis/saveRDS/combined.RDS")
