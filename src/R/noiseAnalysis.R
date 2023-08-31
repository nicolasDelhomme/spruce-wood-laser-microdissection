#' ---
#' title: "Tracheid and Xylem Laser Microdissection Project, 'denoising'"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/akaerkoenen/laser-capture-microdissection-olga")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/akaerkoenen/laser-capture-microdissection-olga")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))

#' Source helpers
source("~/Git/UPSCb/src/R/expressionSpecificityUtility.R")
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' Load the data
mat <- read.csv("analysis/HTSeq/raw-unormalised_merged-technical-replicates_data.csv",
                row.names = 1)

#' Plotting parameters
mar <- par("mar")
pal <- c(1,2,brewer.pal(8,"Dark2"))

#' # Raw data analysis
par(mar=c(8.1,4.1,3.1,0.1))
barplot(colSums(mat),las=2,main="number of reads")
boxplot(log10(mat+1),las=2,ylab="log10(counts+1)",main="gene expression")

#' # Noise analysis

#' ## Dispersion estimation
dge <- DGEList(mat,group = sub("S[1-4]\\.","",colnames(mat)))
dge <- estimateCommonDisp(dge,verbose=TRUE)
plotBCV(dge)

#' ## Normalise
conditions <- colnames(mat)
dds <- DESeqDataSetFromMatrix(
  countData = mat,
  colData = data.frame(condition=conditions),
  design = ~ condition)

dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(mat)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

rld <- rlog(dds)
rld <- assay(rld)
colnames(rld) <- colnames(mat)
rld <- rld[!rowSums(rld==0)==ncol(rld),]
vsd <- varianceStabilizingTransformation(dds)
vsd <- assay(vsd)
colnames(vsd) <- colnames(mat)

dev.null <- sapply(1:ncol(vsd),function(i,r,v){
  heatscatter(r[,i],v[match(rownames(r),rownames(v)),i],
              xlab="rlog",ylab="vst",main=colnames(v)[i])
},rld,vsd)

plot.multidensity(as.data.frame(rld),
                  col=pal[as.integer(factor(colnames(rld)))],
                  legend.x="topright",
                  legend=levels(factor(colnames(rld))),
                  legend.col=pal,
                  legend.lwd=2,
                  main="sample rlog expression distribution",
                  xlab="per gene rlog expression (appx. log2)")

plot.multidensity(as.data.frame(vsd),
                  col=pal[as.integer(factor(colnames(rld)))],
                  legend.x="topright",
                  legend=levels(factor(colnames(rld))),
                  legend.col=pal,
                  legend.lwd=2,
                  main="sample vst expression distribution",
                  xlab="per gene vst expression (appx. log2)")

pc <- prcomp(t(rld))
percent <- round(summary(pc)$importance[2,]*100)

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(colnames(rld)))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",pch=19,
       col=pal,
       legend=levels(factor(colnames(rld))))

#' ### The 2nd and 3rd dims
#' The 2nd dimension (accounting 17% of the overall variance) nor the third 
#' dimension seem to have a biological explanation.
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(factor(colnames(rld)))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("center",pch=19,
       col=pal,
       legend=levels(factor(colnames(rld))))

plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(factor(sub("\\..*","",colnames(rld))))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

#' ## Expression specificity
xsp <- expressionSpecificity(vsd-min(vsd),sub("S[1-4]\\.","",colnames(rld)),output = "complete")
tab <- table(xsp[,"n"])
pander(tab)
barplot(tab,ylab="genes",xlab="number of tissues", main="gene expression observed accross tissues")
hist(xsp[,"score"],breaks = seq(0,1,0.01),xlab="expression specificity",main="")
hist(xsp[xsp[,"n"] == 1,"score"],breaks = seq(0,1,0.01),xlab="expression specificity",main="")
hist(xsp[xsp[,"n"] == 2,"score"],breaks = seq(0,1,0.01),xlab="expression specificity",main="")
hist(xsp[xsp[,"n"] == 3,"score"],breaks = seq(0,1,0.01),xlab="expression specificity",main="")

plot.multidensity(split(xsp[xsp[,"n"]>0,"maxn"],xsp[xsp[,"n"]>0,"n"]))

boxplot(split(xsp[xsp[,"n"]>0,"maxn"],xsp[xsp[,"n"]>0,"n"]))

sapply(lapply(split.data.frame(xsp[xsp[,"n"]>1,c("maxn","score")],xsp[xsp[,"n"]>1,"n"]),as.matrix),function(o){
  heatscatter(o[,1],o[,2])})


pander(colSums(!is.na(xsp[xsp[,"n"] == 1,2:4])))


#' ## A cutoff
vsd <- vsd - min(vsd)
mat <- matrix(apply(expand.grid(1:5,1:9),1,function(ro){sum(rowSums(vsd>ro[1])>ro[2])}),
       nrow=5,ncol=9)

heatmap.2(mat/sum(rowSums(vsd)>0))

pander(mat/sum(rowSums(vsd)>0))

## TODO make the selector a parameter to plot the PCA and a heatmap
dev.null <- apply(expand.grid(1:5,1:9),1,function(ro,sds){
  sel <- rowSums(vsd>ro[1])>ro[2]

  pc <- prcomp(t(vsd[sel,]))
  percent <- round(summary(pc)$importance[2,]*100)
  
  plot(pc$x[,1],
       pc$x[,2],
       xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
       ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
       col=pal[as.integer(factor(sub("S[1-4]\\.","",colnames(rld))))],
       pch=19,
       main=paste0("Principal Component Analysis (",ro[1],",",ro[2],")"),
       sub="variance stabilized counts")
  
  legend("top",pch=19,
         col=pal,
         legend=levels(factor(sub("S[1-4]\\.","",colnames(rld)))))
  
  #' ### The 2nd and 3rd dims
  #' The 2nd dimension (accounting 17% of the overall variance) nor the third 
  #' dimension seem to have a biological explanation.
  plot(pc$x[,2],
       pc$x[,3],
       xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
       ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
       col=pal[as.integer(factor(sub("S[1-4]\\.","",colnames(rld))))],
       pch=19,
       main="Principal Component Analysis",sub="variance stabilized counts")
  
  legend("center",pch=19,
         col=pal,
         legend=levels(factor(sub("S[1-4]\\.","",colnames(rld)))))
  
  
  heatmap.2(vsd[sel,][order(sds[sel],decreasing=TRUE)[1:1000],],
            labRow = NA,trace = "none",cexCol = 0.6 )
  
},apply(vsd,1,sd))
  