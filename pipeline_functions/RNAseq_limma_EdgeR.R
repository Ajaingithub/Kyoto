### RNAseq analysis using Limma and EdgeR

LimmaEdgeR_differential <- function(dds, design0, cm, savedir, logfc = 1.5, p_value_adj = 0.05, grouping_by,group1, group2, column_name_split) {
  # Removing heteroscedascity from count data
  # Heteroskedasticity refers to situations where the variance of the residuals is unequal over a range of measured values.
  # If there is an unequal scatter of residuals, the population used in the regression contains unequal variance, and therefore the analysis results may be invalid.

  # It has been shown that for RNA-seq count data, the variance is not independent of the mean – this is true of raw counts or when transformed
  # to log-CPM values. Methods that model counts using a Negative Binomial distribution assume a quadratic mean-variance relationship. In limma,
  # linear modelling is carried out on the log-CPM values which are assumed to be normally distributed and the mean-variance relationship is
  # accommodated using precision weights calculated by the voom function.
  # When operating on a DGEList-object, voom converts raw counts to log-CPM values by automatically extracting library sizes and normalisation
  # factors from x itself. Additional normalisation to log-CPM values can be specified within voom using the normalize.method argument.

  # des <- DESeq(des) #we will run DESeq to des to get the values for genes like p.values,
  # condition for old and young and Wald statistics
  library(DEFormats)
  desl <- as.DGEList(dds)
  message("Normalizaling Data")
  desl <- calcNormFactors(desl) # This add the effective normalization library size for each sample does not perform normalization
  desl <- estimateGLMCommonDisp(desl, design = design0)

  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4053721/
  # Voom
  # A simple approach to analyzing RNA-seq data would be to input the log-cpm values into the limma lmfit this would be expected to behave well
  # if the counts were all reasonably large, but it ignores the mean-variance trend for lower counts.
  # However, in RNA-seq applications, the count sizes may vary considerably from sample to sample for the same gene. Different samples
  # may be sequenced to different depths, so different count sizes may be quite different even if the cpm values are the same. For this
  # reason, we wish to model the mean-variance trend of the log-cpm values at the individual observation level, instead of applying a
  # gene-level variability estimate to all observations from the same gene.
  # It Control for Type 1 error rate

  # With Quality Weights https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4551905/
  # Removal of high variation samples reduces noise, but at a cost of reducing power, thus limiting our ability to detect biologically
  # meaningful changes. Similarly, retaining these samples in the analysis may not reveal any statistically significant changes due to
  # the higher noise level. A compromise is to use all available data, but to down-weight the observations from more variable samples.
  # We describe a statistical approach that facilitates this by modelling heterogeneity at both the sample and observational levels as
  # part of the differential expression analysis.
  message("Performing Voom and adding weight to the sample")
  desl_Wnorm <- voomWithQualityWeights(desl,
    design0,
    # normalize.method = 'quantile',
    plot = TRUE
  )
  ## No need to remove any technical effect like gender effect

  message("Fitting the model")
  desl_Wnorm_fit <- lmFit(desl_Wnorm, design0,
    # block = desl_Wnorm$targets$batch,
    # correlation = dupcor$consensus,
    method = "robust",
    maxit = 10000
  )
  desl_Wnorm_fit_2 <- contrasts.fit(desl_Wnorm_fit, cm)
  desl_Wnorm_fit_2 <- eBayes(desl_Wnorm_fit_2, proportion = 1 / 10)
  pdf(paste(savedir, "mean_variance_trend.pdf", sep = ""))
  print(plotSA(desl_Wnorm_fit_2, main = "Final model: Mean-variance trend"))
  dev.off()
  results <- decideTests(desl_Wnorm_fit_2, p.value = 0.05) # provide the sample with p-value <0.05

  message("Differential result saved")
  summary(results)
  setwd(savedir)
  write.fit(
    fit = desl_Wnorm_fit_2, results = results,
    file = "All_genes_wnorm_fit.txt",
    adjust = "BH", sep = "\t"
  )
  results.summary <- summary(results)

  message("Making diferential barplot")
  pdf("differential_number_barplot.pdf")
  barplot(results.summary[c(1, 3), ],
    beside = T, las = 3, cex.names = 0.7, ylim = c(0, max(results.summary[c(1, 3), ])),
    legend = rownames(results.summary[c(1, 3), ]), col = c("darkblue", "red")
  )
  dev.off()

  message("Making Volcano plot")
  `%notin%` <- Negate(`%in%`)
  saveDir <- savedir
  dir.create(paste(saveDir, "Differential", sep = ""), showWarnings = FALSE)
  message("saving the top changed genes")
  for (i in 1:length(colnames(cm))) {
    tops <- topTable(desl_Wnorm_fit_2,
      coef = colnames(cm)[i],
      number = Inf,
      sort.by = "none"
    ) # Extract a table of the top-ranked genes from a linear model fit
    all_region <- rownames(tops)
    tops_sig <- subset(tops, adj.P.Val < p_value_adj)
    tops_sig_up <- subset(tops_sig, logFC > logfc)
    write.table(tops_sig_up,
      paste(saveDir, "Differential/", colnames(cm)[i], "_sig_up.txt", sep = ""),
      quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t"
    )
    # write.table(rownames(tops_sig_up),
    #             paste(saveDir,"Differential/",colnames(cm)[i],"_sig_up_genes.txt",sep = ""),
    #             quote = FALSE, row.names = FALSE, col.names = FALSE)

    tops_sig_down <- subset(tops_sig, logFC < -logfc)
    write.table(tops_sig_down,
      paste(saveDir, "Differential/", colnames(cm)[i], "_sig_down.txt", sep = ""),
      quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t"
    )

    # write.table(rownames(tops_sig_down),
    #             paste(saveDir,"Differential/",colnames(cm)[i],"_sig_down_genes.txt",sep = ""),
    #             quote = FALSE, row.names = FALSE, col.names = FALSE)

    write.table(tops,
      paste(saveDir, "Differential/", colnames(cm)[i], "_all_genes_no_filter.txt", sep = ""),
      quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t"
    )

    # all_region_except_differential <- all_region[which(all_region %notin% rownames(tops_sig))]
    # write.table(all_region_except_differential,
    #             paste(saveDir,"Differential/",colnames(cm)[i],"_non_significant_peaks_background.txt",sep = ""),
    #             quote = FALSE, row.names = FALSE, col.names = FALSE)

    #### volcano plot
    source("/mnt/data/projects/pipeline_functions/Volcano_plot.R")
    volcano_plot(
      fit = desl_Wnorm_fit_2, file = colnames(cm)[i], px = 6,
      py = 4, pvalue = p_value_adj, logFC = logfc,
      genes = c(rownames(tops_sig_down), rownames(tops_sig_up)), saveDir = savedir
    )
  }

  norm_counts <- cpm(desl, prior.count = 2, log = TRUE)

  write.table(norm_counts, paste(savedir, "CPM_counts.txt", sep = ""), col.names = T, row.names = T, sep = "\t", quote = F)

  message("Making MA plot")
  source("/mnt/data/projects/pipeline_functions/MA_plot.R")
  for (i in 1:length(colnames(cm))) {
    YF_MA_plot(desl_Wnorm_fit_2, colnames(cm)[i], 12, 2, pvalue = p_value_adj, logfc = logfc, saveDir = savedir)
  }

  save(list = c("desl_Wnorm_fit_2", "desl_Wnorm", "desl", "cm"), file = paste(savedir, "differential_limma.RData", sep = ""))

  # TFs <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/SCENIC_plus/scenicplus/resources/allTFs_hg38.txt")[,1]

  #### Extracting the Transcription factors from the upregulated and downregulated genes
  # tops_sig_up_TF <- tops_sig_up[match(TFs, rownames(tops_sig_up), nomatch = 0),]
  # tops_sig_up_TF_top30 <- tops_sig_up_TF[order(-tops_sig_up_TF[,"logFC"]),] %>% head(30)

  # tops_sig_down_TF <- tops_sig_down[match(TFs, rownames(tops_sig_down), nomatch = 0),]
  # tops_sig_down_TF_top30 <- tops_sig_down_TF[order(tops_sig_down_TF[,"logFC"]),] %>% head(30)

  # ### Combining both the up and downregulated genes
  # tot_TFs <- c(rownames(tops_sig_up_TF_top30), rownames(tops_sig_down_TF_top30))

  # if (length(tot_TFs) > 1) {
  #   rld_TF <- norm_counts[grep(paste("^",tot_TFs,"$",sep = "", collapse="|"), rownames(norm_counts)),]
  #   # rld_TF_req <- rld_TF[,grep("gEVZV",colnames(rld_TF))]
  #   rld_TF_scaled <- t(scale(t(rld_TF)))

  #   Treg_metadata <- data.frame(samples = colnames(rld_TF_scaled))

  #   split_to_dataframe <- function(x) {
  #     split_elements <- strsplit(x, "_")[[1]]
  #     data.frame(t(split_elements))
  #   }

  #   # Convert split elements to dataframe
  #   Treg_metadata_2 <- do.call(rbind, lapply(Treg_metadata$samples, split_to_dataframe))
  #   colnames(Treg_metadata_2) <- c(column_name_split,grouping_by)
  #   Treg_metadata_2$orig.ident <- Treg_metadata$samples

  #   stopifnot(all(Treg_metadata$samples == colnames(rld_TF_scaled))) # if TRUE move forward

  #   library(ComplexHeatmap)
  #   library(circlize)
  #   column_anno <- HeatmapAnnotation(Age = Treg_metadata_2$Age2,
  #                                    col = list(Age = c("P" = "darkolivegreen4", "R" = "firebrick4"))
  #   )

  #   rld_TF_scaled_ordered <- rld_TF_scaled[match(tot_TFs,rownames(rld_TF_scaled)),]
  #   h=Heatmap(rld_TF_scaled_ordered,
  #             name = "Z-score", #title of legend
  #             column_title = "TFs Specific top 60", row_title = "Samples",
  #             row_names_gp = gpar(fontsize = 6), # Text size for row names
  #             column_names_gp = gpar(fontsize = 6),
  #             col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  #             cluster_columns = FALSE,
  #             cluster_rows = FALSE,
  #             top_annotation = column_anno
  #   )

  #   dir.create(paste(savedir,"heatmap",sep = ""), showWarnings = FALSE)
  #   write.table(rld_TF_scaled_ordered, paste(savedir,"heatmap/scaled_ordered.txt",sep = ""), sep = "\t", row.names = T, col.names = T, quote = F)
  #   write.table(rbind(tops_sig_up_TF_top30, tops_sig_down_TF_top30), paste(savedir,"heatmap/TFs_fold_change.txt",sep = ""), sep = "\t", row.names = T, col.names = T, quote = F)


  #   pdf(paste(savedir,"heatmap/top_TFs.pdf",sep = ""))
  #   print(h)
  #   dev.off()
  # }

  #### Heatmap for top 50 gene up and down
  message("Making Heatmap")
  tops_sig_up_top50 <- tops_sig_up[order(-tops_sig_up$logFC), ] %>% head(50)
  tops_sig_down_top50 <- tops_sig_down[order(tops_sig_down$logFC), ] %>% head(50)

  ### Combining both the up and downregulated genes
  tot_genes <- c(rownames(tops_sig_up_top50), rownames(tops_sig_down_top50))

  if (length(tot_genes) > 1) {
    rld_gene <- norm_counts[grep(paste("^", tot_genes, "$", sep = "", collapse = "|"), rownames(norm_counts)), ]
    # rld_gene_req <- rld_gene[,grep("gEVZV",colnames(rld_gene))]
    rld_gene_scaled <- t(scale(t(rld_gene)))

    Treg_metadata <- data.frame(samples = colnames(rld_gene_scaled))

    split_to_dataframe <- function(x) {
      split_elements <- strsplit(x, "-")[[1]]
      data.frame(t(split_elements))
    }

    # Convert split elements to dataframe
    Treg_metadata_2 <- do.call(rbind, lapply(Treg_metadata$samples, split_to_dataframe))
    colnames(Treg_metadata_2) <- c(column_name_split, grouping_by)
    Treg_metadata_2$orig.ident <- Treg_metadata$samples

    stopifnot(all(Treg_metadata$samples == colnames(rld_gene_scaled))) # if TRUE move forward

    library(ComplexHeatmap)
    library(circlize)
    column_anno <- HeatmapAnnotation(
      condition = Treg_metadata_2[,grouping_by],
      col = list(Age = c(group1 = "darkolivegreen4", group2 = "firebrick4"))
    )

    rld_gene_scaled_ordered <- rld_gene_scaled[match(tot_genes, rownames(rld_gene_scaled)), ]
    h <- Heatmap(rld_gene_scaled_ordered,
      name = "Z-score", # title of legend
      column_title = "genes top 100", row_title = "Samples",
      row_names_gp = gpar(fontsize = 6), # Text size for row names
      column_names_gp = gpar(fontsize = 6),
      col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
      cluster_columns = TRUE,
      cluster_rows = FALSE,
      top_annotation = column_anno
    )


    dir.create(paste(savedir, "heatmap", sep = ""), showWarnings = FALSE)
    write.table(rld_gene_scaled_ordered, paste(savedir, "heatmap/scaled_ordered.txt", sep = ""), sep = "\t", row.names = T, col.names = T, quote = F)
    write.table(rbind(tops_sig_up_top50, tops_sig_down_top50), paste(savedir, "heatmap/genes_fold_change.txt", sep = ""), sep = "\t", row.names = T, col.names = T, quote = F)

    pdf(paste(savedir, "heatmap/top_genes.pdf", sep = ""), height = 12, width = 10)
    print(h)
    dev.off()

    h <- Heatmap(rld_gene_scaled_ordered,
      name = "Z-score", # title of legend
      column_title = "genes top 100", row_title = "Samples",
      row_names_gp = gpar(fontsize = 6), # Text size for row names
      column_names_gp = gpar(fontsize = 15),
      col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
      cluster_columns = TRUE,
      cluster_rows = TRUE,
      top_annotation = column_anno
    )

    pdf(paste(savedir, "heatmap/top_genes_clustered.pdf", sep = ""), height = 12, width = 10)
    print(h)
    dev.off()

  }

  message("Please check the ", savedir, " for volcano plot MA plot DIfferential genes")
  return(desl)
}
