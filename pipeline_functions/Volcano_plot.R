## this is for the volcano plot
# fit= desl_wnrom_fit
# file = colname of the confusion matrix
# px and py are the position of x-axis and y-axis where we need to add the up and down regulated genes.
# genes = genes to be marked on the 
volcano_plot <- function(fit,file,px=6,py=20,pvalue=0.05,logFC=1.5,genes=NULL,saveDir){
  library(dplyr)
  library(ggrepel)
  library(ggplot2)
  tops<-topTable(fit,coef = file, number = Inf, sort.by = "none") #Extract a table of the top-ranked genes from a linear model fit
  tops_plot <- dplyr::select(tops,logFC, adj.P.Val, P.Value) 
  tops_plot <- mutate(tops_plot,sig=ifelse(tops_plot$adj.P.Val<pvalue,"Sig", "Not Sig"))
  tops_plot$significant = tops_plot$sig
  tops_plot_sig = tops_plot[grep("^Sig$",tops_plot$significant),]
  tops_plot_sig_up = tops_plot_sig[which(tops_plot_sig$logFC>logFC),]
  tops_plot_sig_down = tops_plot_sig[which(tops_plot_sig$logFC < -logFC),]
  
  tops_plot[which(rownames(tops_plot) %in% rownames(tops_plot_sig_up)),"significant"]  <- "Up_sig"
  tops_plot[which(rownames(tops_plot) %in% rownames(tops_plot_sig_down)),"significant"]  <- "Down_sig"
  
  tops_plot$significant <- gsub("^Sig$","Not Sig",tops_plot$significant)
  volc = ggplot(tops_plot, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(col=significant), size=1) +
    geom_vline(xintercept = logFC, linetype="dashed", color="red") + 
    geom_vline(xintercept = -logFC, linetype="dashed", color="red") + 
    geom_hline(yintercept = -log10(pvalue), linetype="dashed", color="red") + 
    scale_color_manual(values=c("blue", "grey","red")) + 
    annotate("text", x = px, y = py, label = paste0("Up=",nrow(tops_plot_sig_up))) + 
    annotate("text", x = -px, y = py, label = paste0("Down=",nrow(tops_plot_sig_down))) +
    ggtitle(colnames(cm)) + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = 'white', colour = 'black'))
  # panel.grid.minor = element_line(colour = "grey"),
  # panel.grid.major = element_line(colour = "grey"))
  
  dir.create(paste(saveDir,"volcano",sep = ""), showWarnings = FALSE)
  if(length(genes) == 0){
    pdf(paste(saveDir,"volcano/",file,"_volcano_plot.pdf",sep = ""),width = 10, height = 10)
    print(volc)
    dev.off()  
  } else{
    x_limits <- c(px, NA)
    genes_label<- tops_plot[which(rownames(tops_plot) %in% toupper(genes)),]
    p2 <- volc+geom_text_repel(data=tops_plot[which(rownames(tops_plot) %in% toupper(genes)),],
                               aes(label=rownames(genes_label)), size = 3)
                               # , xlim = x_limits)
    pdf(paste(saveDir,"volcano/",file,"_volcano_plot.pdf",sep = ""),width = 10, height = 10)
    print(p2)
    dev.off()  
  }
}


