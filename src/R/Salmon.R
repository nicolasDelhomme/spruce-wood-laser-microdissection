#' ---
#' title: "LMD transcripts (protein coding) expression analysis"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/facility/trinity/LMD-trinity")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/facility/trinity/LMD-trinity")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' Create a palette
pal <- brewer.pal(8,"Dark2")

#' Register the default plot margin
mar <- par("mar")

#' # Raw data
#' ## Loading
#' Read the sample information
samples <- read.csv("~/Git/UPSCb/projects/spruce-LMD/doc/ENA/Lignifying-xylem-of-Norway-spruce_Samples.csv")

#' Read the results
files=dir("Salmon",recursive=TRUE,pattern="quant.sf",full.names = TRUE)
names(files) <- gsub("Salmon/|_sortmerna.*","",files)
txi <- tximport(files, type="salmon", txOut=TRUE,
                countsFromAbundance="no")


#' combine the tech reps
cts <- txi$counts

counts <- sapply(split.data.frame(t(cts),sub("_L00[1-4]$","",colnames(cts))),colSums)
head(counts)

#cts <- cts[rowSums(cts) > 0,]

#' Display the per-gene mean expression
#' 
#' i.e. the mean raw count of every 
#' gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative coverage is as expected, around 200X
plot(density(log10(rowMeans(cts))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' What do we want from here
# tissue specificity
# heatmap
# some filtering
# split between known genes and the rest
# use the GMAP to tx2gene

#' The same is done for the individual
#' samples colored by condition 
#' 
#' The observed distribution is strikingly similar. The re-sequencing of the two
#' undersampled samples (13 and 14), which are both biological replicates brought
#' them to a similar depth as the others. 
plot.multidensity(log10(count.table),
                  col=pal[as.integer(samples$Description)],
                  legend.x="topright",
                  legend=levels(samples$Description),
                  legend.col=pal[1:nlevels(samples$Description)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' # Data normalisation 
#' ## Preparation
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate

#' Create the dds object, without giving any prior on the design
conditions <- colnames(count.table)
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions,
                       genotype=samples$Genotype,
                       age=samples$Age.in.days,
                       temp=samples$Temp),
  design = ~ condition)

#' Save the object for later analysis
dir.create(file.path("analysis","DESeq2"),recursive=TRUE,showWarnings=FALSE)
save(dds,file = file.path("analysis/DESeq2/","DESeq2-object.rda"))

#' Check the size factors (i.e. the sequencing library size effect)
#' There is some known extreme variation (2 samples under and on over-sampled),
#' however, for the sole purpose of quality assessment, a Variance Stabilizing 
#' Transformation can be used (over a Relative Log Transformation)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(count.table)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)
write.csv(vst,"analysis/HTSeq/library-size-normalized_variance-stabilized_data.csv")

#' Validate the VST 
#' 
#' Visualize the corrected mean - sd relationship. It is almost linear,
#' meaning we can assume homoscedasticity.
#' The slight initial trend / bump is due to genes having few counts in
#' a few subset of the samples and hence having a higher variability. This is
#' expected.
meanSdPlot(vst[rowSums(count.table)>0,])

#' # QC on the normalised data
#' 
#' ## PCA
#' 
#' First perform a Principal Component Analysis (PCA) of the data
#'to do a quick quality assessment; i.e. replicate should cluster
#' and the first 2-3 dimensions should be explainable by biological means.
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' ### 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples$Description)],
              pch=19)
legend("topleft",pch=19,
       col=pal[1:nlevels(samples$Description)],
       legend=levels(samples$Description))
par(mar=mar)

#' ### first two dimensions
#' 
#' The biological replicates cluster together.
#' The first dimension separates the 2 "final" temperature, whereas the second dimension
#' separates (for the lower final temperature) samples that were grown in the cold 
#' from samples that were shifted to the cold. 2 Col-0 grown 9 days at 16 degree are behind
#' the legend, next to the their third replicate
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples$Description)],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("top",pch=19,
       col=pal[1:nlevels(samples$Description)],
       legend=levels(samples$Description))


#' Just to make the point clearer, plotting the temperature categorical data
#' (C,W,CW,WC)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(samples$Temp))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",pch=19,
       col=pal[1:nlevels(factor(samples$Temp))],
       legend=levels(factor(samples$Temp)))

#' And the same for the temperature at collection
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(samples$Temperature.at.collection))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",pch=19,
       col=pal[1:nlevels(factor(samples$Temperature.at.collection))],
       legend=levels(factor(samples$Temperature.at.collection)))

#' And the same for the age
#' Age is possibly a confounding factor
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(samples$Age.in.days))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",pch=19,
       col=pal[1:nlevels(factor(samples$Age.in.days))],
       legend=levels(factor(samples$Age.in.days)))

#' ### 2nd and 3rd dims
#' 
#' Which to some extend separates the pcp from the Col-0 genotypes
#' 
#' 3 Col-0 samples grown at 23 degrees for 9 days are behind the legend
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(samples$Description)],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",pch=19,
       col=pal[1:nlevels(samples$Description)],
       legend=levels(samples$Description))

#' And the same plotting the genotype
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(factor(samples$Genotype))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",pch=19,
       col=pal[1:nlevels(factor(samples$Genotype))],
       legend=levels(factor(samples$Genotype)))

#' And finally the same for the age
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(factor(samples$Age.in.days))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",pch=19,
       col=pal[1:nlevels(factor(samples$Age.in.days))],
       legend=levels(factor(samples$Age.in.days)))

#' ## Heatmap
#' The 1000 most variable genes are selected and plotted as a heatmap
dimnames(vst) <- dimnames(count.table)
sel <- order(apply(vst,1,sd),decreasing=TRUE)[1:1000]

#' Labelled with sample name and using expression values (approx. log2)
heatmap.2(vst[sel,],labRow = NA,trace = "none",cexCol = 0.6 )

#' Labelled with genotype and temperature and using per gene expression z-scores.
heatmap.2(vst[sel,],labRow = NA,trace = "none",cexCol = 0.6,scale = "row",
          labCol = paste(samples$Genotype,samples$Temp))

#' # Conclusion
#' The raw quality of the data appears good. The amount of read mapping uniquely 
#' to coding genes is low (36 percent on average), but this could be due to the 
#' mRNA selection protocol, which would have kept all non-polyA transcripts. Re-sampling the
#' 2 originally under-sampled samples has not changed the QA conclusion, but has clarified
#' the PCA and heatmaps. This presumes of an increased power for the DE analysis.
#' 
#' The data normalisation gives satisfying results (as far as the VST fit is concerned). 
#' The PCA and the heatmap identifies no samples as possible outliers and the three first
#' dimensions can be explained by the experimental design. 
#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
