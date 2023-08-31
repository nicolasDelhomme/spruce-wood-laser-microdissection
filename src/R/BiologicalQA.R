#' ---
#' title: "Tracheid and Xylem Laser Microdissection Project, Biological QA"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/laser-capture-microdissection-olga/non-strand-spec/htseq")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/laser-capture-microdissection-olga/non-strand-spec/htseq")
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
# source("~/Git/UPSCb/src/R/plotMA.R")
# source("~/Git/UPSCb/src/R/volcanoPlot.R")

#' Create a palette
pal <- brewer.pal(8,"Dark2")

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
dir.create(file.path("analysis","HTSeq"),recursive = TRUE, showWarnings = FALSE)
write.csv(count.table,"analysis/HTSeq/raw-unormalised-data.csv")

#' ## HTSeq stats
#' Extract the HTSeq stat lines
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

#' as percentages
#' 
#' The S1 ray cells sample seem to be an outlier. It has a lot of reads mapping
#' outside known features. A quick look at the alignments in a genome browser
#' revealed that many non genic regions have been deeply sequenced. It suggests
#' that some DNA might have made it in the RNA extraction and has been amplified
#' during the library prep. I suggest DNA contamination, as in the handful cases
#' I have investigated, there is no other feature of interest (such as transposable
#' elements) that could reasonably be transcribed. The genomic coverage and the 
#' GC content of these region is also standard.
#' 
pander(apply(count.stats,2,function(co){round(co*100/sum(co))}))

#' as a plot
#' 
#' Apart from the S1 ray cells sample, the proportion looks adequate. That every
#' technical replicate (resulting of the sequencing lane split) is very similar to
#' one another indicates that the issue is not a sequencing technicality. The
#' fact that the S1 whole sections have less aligned reads - around 4M per tech.
#' rep. is due to the fact that that library had a 20% 28S rRNA "contamination".
col <- pal[1:nrow(count.stats)]
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=col,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats) + c(0,4e+6),cex.names=.6)
legend("top",fill=col,legend=gsub("_"," ",rownames(count.stats)),bty="n",cex=0.8)
par(mar=mar)

#' ## Raw expression summary
#' ### Percentage of aligned reads
#' The median percentage of aligned reads is 64%, but the S1 ray cells sample 
#' stands out dramatically, as already observed.
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
#' i.e. the mean raw count of every 
#' gene across samples is calculated
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
#' replicates. The samples however have different coverage, S1 ray cells and 
#' whole sections being "undersampled" whereas S2 xylem tracheids is slightly
#' "oversampled" in comparison to the remaining 3 samples.
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
#' There is some large variation, so a Variance Stabilizing Transformation (VST) is
#' here not optimal (over a Regularised Log Transformation), but it is significantly
#' faster. As the results probably won't be significantly affected and as it is
#' only for QA purposes, it is OK to VST over RLT.
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
write.csv(vst,"analysis/HTSeq/library-size-normalized_variance-stabilized_data.csv")

#' ## VST validation
#' 
#' We visualize the corrected mean - sd relationship to assess whether homoscedasticity
#' can be assumed. The red line should ideally be linear, and the standard 
#' deviation (sd) should be contained around 0.5.
#' 
#' Here, it is surely not linear, the influence of the S1 ray cells 
#' sample is evident. It is very likely undersampled in comparison to the 
#' others and is most likely responsible for the peak around 32,000.
#' We cannot assume homoscedasticity, however, the overall sd is controlled
#' below 1.0 and that shoud be enough for QA purposes.
meanSdPlot(vst[rowSums(count.table)>0,], ylim = c(0,2.5))

#' # QC on the normalised data
#' 
#' ## Principal Component Analysis (PCA)
#' We perform a PCA of the data
#'  to do a quick quality assessment; 
#' i.e. replicate should cluster together
#' and the first 2-3 dimensions should 
#' be explainable by biological means.
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' ### PCA 3 first dimensions
#' 
#' Here we do plot the 24 technical samples, but as can be seen,
#' the technical samples are undistinguishable from one another.
#' On the other hand, the 2 samples identified earlier: S1 ray cells and
#' S2 xylem tracheids are the ones that explain the first (54%) and third (10%)
#' dimension of the data - i.e. how much they account for the overall variability.
#' The second dimension (21%) is better seen in the next plot (in 2 dimensions).
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
       col=pal[1:6],
       legend=levels(factor(colnames(count.table))))
par(mar=mar)

#' ### The first two dimensions
#' Clearly the S1 ray cells sample obfuscate the PCA, it is responsible for 
#' the first dimension. And the S1 tracheids is responsible for the second
#' dimension.
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(colnames(count.table)))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",pch=19,
       col=pal[1:6],
       legend=levels(factor(colnames(count.table))))

#' ### The 2nd and 3rd dims
#' The S1 and S2 Xylem tracheids samples also seems to be outliers, 
#' responsible for the second and third dimensions
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(factor(colnames(count.table)))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",pch=19,
       col=pal[1:6],
       legend=levels(factor(colnames(count.table))))

#' Only further dimensions, e.g. from the 4th on - here 4th and 5th -
#' may not be due to technicalities. It is definitely not what we would 
#' expect to see and is probably indicative of sampling prepapration issues
#' rather than true biology.
plot(pc$x[,4],
     pc$x[,5],
     xlab=paste("Comp. 4 (",percent[4],"%)",sep=""),
     ylab=paste("Comp. 5 (",percent[5],"%)",sep=""),
     col=pal[as.integer(factor(colnames(count.table)))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",pch=19,
       col=pal[1:6],
       legend=levels(factor(colnames(count.table))))


#' ## Heatmap
#' The 1000 most variable genes are selected from the 3 samples not identified
#' as "outliers" in the above - i.e. both whole sections and S2 ray cells - and 
#' plotted as a heatmap. This is of course probably still biased, but it seems to
#' take away a lot of the techinal artefacts; i.e. the ray cells cluster
#' together, as do the whole sections. Sadly, the xylem tracheids are behaving
#' differently, the S1 being closer to ray cells, the S2 closer to whole 
#' sections. Could this be an artefact of the microdissection?
dimnames(vst) <- dimnames(count.table)
sel <- order(apply(vst[,!grepl("S1-ray|tracheids",colnames(vst))],1,sd),decreasing=TRUE)[1:1000]
heatmap.2(vst[sel,],labRow = NA,trace = "none",cexCol = 0.6 )

#' # Post-processing
#' From the QA, it is evident that the technical replicates are identical and hence that we cabn merge them.
#' The resulting table is exported as a csv file.
m.vst <- sapply(unique(colnames(vst)),function(n,m){
  rowSums(m[,colnames(m)==n])
},vst)
write.csv(m.vst,"analysis/HTSeq/library-size-normalized_variance-stabilized_merged-technical-replicates_data.csv")

#' # Conclusion
#' The raw quality of the data has limitations. Samples have different sequencing
#' depth that affect their expression profile (see the raw data expression density
#' plots). The technical replicates are on the other hand perfectly identical 
#' and hence can be merged. The sample S1 ray cells is definitely undersampled,
#' and this is due to a low alignment rates to known features (i.e. genes) 
#' presumably due to a DNA contamination during the library preparation.
#' 
#' The data normalisation does not gives satisfying results (as far as the VST 
#' fit is concerned), as it is impeded by the bias specific to several samples. 
#' Hence the PCA major components do not reveal any significant biological factor.
#' On the contrary, the PCA reveals sample specific variation rather than biologically
#' relevant information; i.e. tissue do not cluster together. It appears that in 
#' addition to the S1 ray cells sample, both the tracheids samples suffer from 
#' some technical biases.
#' 
#' However, the heatmap, when using the 3 samples that do not seem affected by technicalities
#' to select for the 1000 most variable genes shows a more interesting picture, where
#' tissues group together, at the exception of the tracheids.
#' 
#' Doing a differential analysis on the dataset at present would be
#' severely limited. The noise from the outliers would mask most of the real 
#' effects and the lack of a 3 replicate would prevent an adequate estimation
#' of the parameters of the negative binomial distribution - the distribution that
#' is used to model RNA-Seq data - and lead to an increase of the false positive
#' rate and a decrease of the detection power. I.e. there would be less differential 
#' expression candidate genes predicted, a larger proportion of which would be 
#' false calls. The data as is is known can however already be used to assess
#' expression of the genes of interest (lignin biosynthesis), but the technical
#' biases should be kept in mind while doing this.
#' ```{r empty, eval=FALSE,echo=FALSE}
#' ```
#' 
#' # Follow up
#' 
#' I would suggest - if money and material permits - to try to get a 3rd and 
#' ideally a 4th replicate (since for example both tracheid samples seem to be
#' somewhat technically biased). Meanwhile, I will generate the regularised log
#' transformed data, including and excluding the S1 ray cells sample, to provide
#' more robust data for the "manual" investigation of the genes of interest.
#' 
#' ```{r empty, eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r sessionInfo, echo=FALSE}
#' sessionInfo()
#' ```
#' 
