#' ---
#' title: "Tracheid and Xylem Laser Microdissection Project, combined biological QA"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/akaerkoenen/laser-capture-microdissection-olga/htseq")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/akaerkoenen/laser-capture-microdissection-olga/htseq")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(vsn))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' Create a palette
pal <- c(1,2,brewer.pal(8,"Dark2"))

#' Register the default plot margin
mar <- par("mar")

#' Read the sample information
samples <- read.csv("~/Git/UPSCb/projects/spruce-LMD/doc/Lignifying-xylem-of-Norway-spruce_Samples.csv")

#' Read the HTSeq files in a matrix
res <- mclapply(dir(".",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=max(mcaffinity()))

#' and rename the samples
names(res) <- sapply(dir(".",pattern="*.txt"),function(n,s){
  gsub(" ","-",s[grep(sub("_sortmerna_trimmomatic_STAR.txt","",n),s$FileName)[1],"SampleName"])
},samples)

#' # Raw Data QC analysis
addInfo <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",2))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]
dir.create(file.path("..","analysis","HTSeq"),recursive = TRUE, showWarnings = FALSE)
write.csv(count.table,"../analysis/HTSeq/raw-unormalised-data.csv")

#' ## HTSeq stats
#' Extract the HTSeq stat lines
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

#' as percentages
#' 
#' The S4 whole sections sample seem to be an outlier. It has a lot of reads mapping
#' outside known features. All other samples have comparatively less "aligned"
#' reads than the previous sequencing runs
#' 
pander(apply(count.stats,2,function(co){round(co*100/sum(co))}))

#' as a plot
#' 
#' Apart from the S4 whole sections sample, the proportion looks equivalent. That every
#' technical replicate (resulting of the sequencing lane split) is very similar to
#' one another indicates that the issue is not a sequencing technicality. 
#' 
col <- pal[1:nrow(count.stats)]
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=col,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats) + c(0,4e+6),cex.names=.6)
legend("top",fill=col,legend=gsub("_"," ",rownames(count.stats)),bty="n",cex=0.8,horiz=TRUE)
par(mar=mar)

#' ## Raw expression summary
#' ### Percentage of aligned reads
#' The median percentage of aligned reads is 50% but the distribution is
#' fairly large.
boxplot(unlist(count.stats["aligned",]/colSums(count.stats)),
        main="aligned reads",ylab="percent aligned",ylim=c(0,1))

#' ### Amount of un-expressed
sel <- rowSums(count.table) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(count.table),digits=1),
        sum(sel),
        nrow(count.table))

#' ### Overall average expression
#' Display the per-gene mean expression
#' 
#' i.e. the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The curve looks usual but is slightly more enriched for lower counts.
#' It has a peak around 1.5 on a log10 scale, so an average gene coverage of 30X
#' per technical replicates. If the replicates can be merged - which seem to be
#' the case - that would mean an average 120X average coverage which is ideal for 
#' a Differential Expression (DE) analysis.
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' ### Per sample expression
#' The same as above is done for the individual
#' samples colored by sample type
#' 
#' The observed distribution is strickingly identical between technical 
#' replicates. The samples however have different coverage - i.e. 
#' different effective sequencing depth.
plot.multidensity(log10(count.table),
                  col=pal[as.integer(factor(colnames(count.table)))],
                  legend.x="topright",
                  legend=levels(factor(colnames(count.table))),
                  legend.col=pal[1:6],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' # Data normalisation 
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample type and sample experiment
#'
conditions <- colnames(count.table)
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions),
  design = ~ condition)

#' ## Library size Factor (a.k.a. scaling factor)
#' Here, we check the size factors (i.e. the sequencing library size effect).
#' There is some large variation, making a VST transformation borderline. A rlog
#' transformation should be prefered, but for QC purposes, a VST is sufficient.
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(count.table)
#' The size factors calculated by DESeq2 display the same information as the
#' raw expression density plots above; which is expected.
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation (VST)
#' The data is transformed to be homoscedastic, and shifted so that non-expressed
#' genes have a value of 0. The resulting table, containing all replicates, is 
#' saved.
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)
vst <- vst - min(vst)
write.csv(vst,"../analysis/HTSeq/library-size-normalized_variance-stabilized_data.csv")

#' ## VST validation
#' 
#' We visualize the corrected mean - sd relationship to assess whether 
#' homoscedasticity can be assumed. The red line should ideally be linear, 
#' and the standard deviation (sd) should be contained below 0.5.
#' 
#' As expected the normalisation is not optimal, the homoscedasticity is reached
#' only for highy expressed genes.
meanSdPlot(vst[rowSums(count.table)>0,])

#' # QC on the normalised data
#' 
#' ## Principal Component Analysis (PCA)
#' We perform a PCA of the data to do a quick quality assessment; 
#' i.e. replicate should cluster together and the first 2-3 dimensions should 
#' be explainable by biological means.
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' ### PCA 3 first dimensions
#'
#' Here we do plot all technical rrplicates, but as can be seen,
#' they are undistinguishable from one another. Sample S1 ray cells is an 
#' obvious outlier.
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(factor(colnames(count.table)))],
              pch=19)
legend("top",pch=19,
       col=pal,
       legend=levels(factor(colnames(count.table))))
par(mar=mar)

#' ### The first two dimensions
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(colnames(count.table)))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",pch=19,
       col=pal,
       legend=levels(factor(colnames(count.table))))

#' ### The 2nd and 3rd dims
#' The 2nd dimension (accounting 17% of the overall variance) nor the third 
#' dimension seem to have a biological explanation.
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(factor(colnames(count.table)))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("center",pch=19,
       col=pal,
       legend=levels(factor(colnames(count.table))))

#' ## Heatmap
#' The 1000 most variable genes are selected and plotted as a heatmap. The 
#' results recapitulate what we have observed in the PCA.
#' 
dimnames(vst) <- dimnames(count.table)
sel <- order(apply(vst,1,sd),decreasing=TRUE)[1:1000]
heatmap.2(vst[sel,],labRow = NA,trace = "none",cexCol = 0.6 )

#' Looking at the PCA of the same genes does not reveal a clearer pattern.
pc <- prcomp(t(vst[sel,]))
percent <- round(summary(pc)$importance[2,]*100)

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(colnames(count.table)))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 1 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[3],"%)",sep=""),
     col=pal[as.integer(factor(colnames(count.table)))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",pch=19,
       col=pal,
       legend=levels(factor(colnames(count.table))))

plot.multidensity(as.data.frame(vst[rowSums(vst)>0,]),
                  col=pal[as.integer(factor(colnames(count.table)))],
                  legend.x="topright",
                  legend=levels(factor(colnames(count.table))),
                  legend.col=pal[1:6],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (vst)")

#' # Conclusion
#' From the QA, it is evident that the technical replicates are identical and 
#' hence that we can merge them. However there are technicalities which prevent
#' an immediate use of the data for differential expression. Presumably a lot of
#' the technical variability has been introduced by the amplification protocol
#' used (very low amount of starting material was available).
#' 
m.count.table <- sapply(unique(colnames(count.table)),function(n,m){
  rowSums(m[,colnames(m)==n])
},count.table)
write.csv(m.count.table,"../analysis/HTSeq/raw-unormalised_merged-technical-replicates_data.csv")
#' ```{r empty, eval=FALSE,echo=FALSE}
#' ```
#' 
#' # Follow up
#' 
#' I will next try different approaches to identify an adequate cutoff for 
#' reducing the noise, including looking at estimates of dispersion, relative
#' mean difference, expression specificity (by sample and tissues).
#' 
#' ```{r empty, eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r sessionInfo, echo=FALSE}
#' sessionInfo()
#' ```
#' 
