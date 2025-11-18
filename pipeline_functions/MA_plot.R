# This is for the MA plot
# fit = desl_Wnorm_fit 
# file = comparison i.e. colname of the confusion matrix
# px and py are the position of x-axis and y-axis where we need to add the up and down regulated genes. 
YF_MA_plot <- function(fit,file,px,py,pvalue=0.05,logfc=1.5,saveDir){
  # xloc <- 12   #position of the text for x axis
  # yloc <- 2   #position of the text for y axis
  dir.create(paste(saveDir,"MA_plots",sep=""),showWarnings = FALSE)
  pdf(paste(saveDir,"MA_plots/",file,"_MA_plot.pdf",sep = ""),width = 10,height = 8)
  tops<-topTable(fit, 
                 coef = file, 
                 number = Inf, 
                 sort.by = "none") #Extract a table of the top-ranked genes from a linear model fit
  tops_sig<-subset(tops,adj.P.Val<pvalue)
  tops_sig_up<-subset(tops_sig,logFC>logfc)
  tops_sig_down<-subset(tops_sig,logFC < -logfc)
  tops_sig_mod<-tops_sig_up
  tops_sig_mod<-rbind(tops_sig_up,tops_sig_down)
  status <- row.names(tops)%in%row.names(tops_sig_mod)
  values <- c("FALSE","TRUE")
  col <- c("black","red")
  attr(status,"values") <- values
  attr(status,"col") <- col
  limma::plotMA(fit, coef = file, 
                status = status,
                cex=0.8,main = file,
                legend=FALSE)
  message("Up=",nrow(tops_sig_up))
  message("Down=",nrow(tops_sig_down))
  text(px,py,paste0("Up=",nrow(tops_sig_up)))
  text(px,-py,paste0("Down=",nrow(tops_sig_down)))
  dev.off()
}
