#' ---
#' title: "Tracheid and Xylem Laser Microdissection Project, DESeq2 analysis"
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
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(vsn))

#' Source Helpers
source("~/Git/UPSCb/src/R/plotMA.R")
source("~/Git/UPSCb/src/R/volcanoPlot.R")
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' palette
pal <- brewer.pal(8,"Dark2")

#' Read in the data
count.table <- read.csv("analysis/HTSeq/raw-unormalised_merged-technical-replicates_data.csv",
                        row.names = 1)
tissues <- factor(sub("S[1-4]\\.","",colnames(count.table)))

#' # DESeq2
#' ## Instantiate
#' ### With the S1 ray cells

dds <- DESeqDataSetFromMatrix(countData = count.table,
                              colData=DataFrame(
                                tissue=tissues),
                              design=~tissue)

dds <- DESeq(dds)

plotDispEsts(dds)

alpha=0.01

r.vs.w <- results(dds,contrast = c("tissue","ray.cells","whole.sections"))
volcanoPlot(r.vs.w, alpha = alpha)
g.r.vs.w <- rownames(r.vs.w[!is.na(r.vs.w$padj) & r.vs.w$padj <= alpha,])

r.vs.x <- results(dds,contrast = c("tissue","ray.cells","xylem.tracheids"))
volcanoPlot(r.vs.x, alpha = alpha)
g.r.vs.x <- rownames(r.vs.x[!is.na(r.vs.x$padj) & r.vs.x$padj <= alpha,])

x.vs.w <- results(dds,contrast = c("tissue","xylem.tracheids","whole.sections"))
volcanoPlot(x.vs.w, alpha = alpha)
g.x.vs.w <- rownames(x.vs.w[!is.na(x.vs.w$padj) & x.vs.w$padj <= alpha,])

#' #### Venn Diagram
plot.new()
grid.draw(venn.diagram(
  list(g.r.vs.w,g.r.vs.x,g.x.vs.w),
  filename = NULL,col=pal[1:3],
  category.names = c("ray-whole","ray-tracheid","tracheid-whole")))

#' ### Without
ddsW <- DESeqDataSetFromMatrix(countData = count.table[,-1],
                               colData=DataFrame(
                                 tissue=tissues[-1]),
                               design=~tissue)

#' Do the variance stabilisation
vsd <- varianceStabilizingTransformation(ddsW, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)[-1]
vst <- vst - min(vst)
meanSdPlot(vst[rowSums(count.table[,-1])>0,])
write.csv(vst,"analysis/HTSeq/library-size-normalized_variance-stabilized_merged-technical-replicates_noS1ray_data.csv")

#' Do the differential expression
ddsW <- DESeq(ddsW)

plotDispEsts(ddsW)

alpha=0.01

r.vs.w.W <- results(ddsW,contrast = c("tissue","ray.cells","whole.sections"))
volcanoPlot(r.vs.w.W, alpha = alpha)
g.r.vs.w.W <- rownames(r.vs.w.W[!is.na(r.vs.w.W$padj) & r.vs.w.W$padj <= alpha,])

r.vs.x.W <- results(ddsW,contrast = c("tissue","ray.cells","xylem.tracheids"))
volcanoPlot(r.vs.x.W, alpha = alpha)
g.r.vs.x.W <- rownames(r.vs.x.W[!is.na(r.vs.x.W$padj) & r.vs.x.W$padj <= alpha,])

x.vs.w.W <- results(ddsW,contrast = c("tissue","xylem.tracheids","whole.sections"))
volcanoPlot(x.vs.w.W, alpha = alpha)
g.x.vs.w.W <- rownames(x.vs.w.W[!is.na(x.vs.w.W$padj) & x.vs.w.W$padj <= alpha,])

#' #### Venn Diagram
plot.new()
grid.draw(venn.diagram(
  list(g.r.vs.w.W,g.r.vs.x.W,g.x.vs.w.W),
  filename = NULL,col=pal[1:3],
  category.names = c("ray-whole","ray-tracheid","tracheid-whole")))

#' ### Ray cells vs. whole sections
plot.new()
grid.draw(venn.diagram(
  list(g.r.vs.w,g.r.vs.w.W),
  filename = NULL,col=pal[1:2],
  category.names = c("with S1","without S1")))

#' ### Ray cells vs. tracheids
plot.new()
grid.draw(venn.diagram(
  list(g.r.vs.x,g.r.vs.x.W),
  filename = NULL,col=pal[1:2],
  category.names = c("with S1","without S1")))

#' ### Tracheids vs. whole section
plot.new()
grid.draw(venn.diagram(
  list(g.x.vs.w,g.x.vs.w.W),
  filename = NULL,col=pal[1:2],
  category.names = c("with S1","without S1")))

#' ## S1 ray cell sample further assessment
#' ### Fold changes: ray vs. whole
d.r.vs.w.W <- r.vs.w.W[!is.na(r.vs.w.W$padj) & r.vs.w.W$padj <= alpha,]
d.r.vs.w <- r.vs.w[!is.na(r.vs.w$padj) & r.vs.w$padj <= alpha,]

g <- intersect(rownames(d.r.vs.w.W),rownames(d.r.vs.w))

l2fc.W <- d.r.vs.w.W[! rownames(d.r.vs.w.W) %in% g,"log2FoldChange"]
l2fc.w <- d.r.vs.w[! rownames(d.r.vs.w) %in% g,"log2FoldChange"]
l2fc.W.c <- d.r.vs.w.W[rownames(d.r.vs.w.W) %in% g,"log2FoldChange"]
l2fc.w.c <- d.r.vs.w[rownames(d.r.vs.w) %in% g,"log2FoldChange"]
plot.multidensity(list(l2fc.W,l2fc.w,l2fc.W.c,l2fc.w.c,
                       c(l2fc.W,l2fc.W.c),
                       c(l2fc.w,l2fc.w.c)),xlab="log2FC",
                  legend = c("spec. no S1","spec S1",
                             "common no S1", "common S1",
                             "all no S1","all S1"),
                  legend.x="topleft")

#' ### Fold changes: ray vs. tracheids
d.r.vs.x.W <- r.vs.x.W[!is.na(r.vs.x.W$padj) & r.vs.x.W$padj <= alpha,]
d.r.vs.x <- r.vs.x[!is.na(r.vs.x$padj) & r.vs.x$padj <= alpha,]

g <- intersect(rownames(d.r.vs.x.W),rownames(d.r.vs.x))

l2fc.W <- d.r.vs.x.W[! rownames(d.r.vs.x.W) %in% g,"log2FoldChange"]
l2fc.w <- d.r.vs.x[! rownames(d.r.vs.x) %in% g,"log2FoldChange"]
l2fc.W.c <- d.r.vs.x.W[rownames(d.r.vs.x.W) %in% g,"log2FoldChange"]
l2fc.w.c <- d.r.vs.x[rownames(d.r.vs.x) %in% g,"log2FoldChange"]
plot.multidensity(list(l2fc.W,l2fc.w,l2fc.W.c,l2fc.w.c,
                       c(l2fc.W,l2fc.W.c),
                       c(l2fc.w,l2fc.w.c)),xlab="log2FC",
                  legend = c("spec. no S1","spec S1",
                             "common no S1", "common S1",
                             "all no S1","all S1"),
                  legend.x="topleft")

#' ### Fold changes: tracheids vs. whole
d.x.vs.w.W <- x.vs.w.W[!is.na(x.vs.w.W$padj) & x.vs.w.W$padj <= alpha,]
d.x.vs.w <- x.vs.w[!is.na(x.vs.w$padj) & x.vs.w$padj <= alpha,]

g <- intersect(rownames(d.x.vs.w.W),rownames(d.x.vs.w))

l2fc.W <- d.x.vs.w.W[! rownames(d.x.vs.w.W) %in% g,"log2FoldChange"]
l2fc.w <- d.x.vs.w[! rownames(d.x.vs.w) %in% g,"log2FoldChange"]
l2fc.W.c <- d.x.vs.w.W[rownames(d.x.vs.w.W) %in% g,"log2FoldChange"]
l2fc.w.c <- d.x.vs.w[rownames(d.x.vs.w) %in% g,"log2FoldChange"]
plot.multidensity(list(l2fc.W,l2fc.w,l2fc.W.c,l2fc.w.c,
                       c(l2fc.W,l2fc.W.c),
                       c(l2fc.w,l2fc.w.c)),xlab="log2FC",
                  col=pal[1:6],
                  legend = c("spec. no S1","spec S1",
                             "common no S1", "common S1",
                             "all no S1","all S1"),
                  legend.x="topleft")

#' ## Export
#' We export the results from the analysis ignoring sample S1-ray-cells.
#' Indeed, as exemplified in the multi-density plots above, the sample S1-ray
#' introduces a bias in the analysis (favoring positive logFC in the comparisons).
#' The last multidensity plot compares the tracheids with the whole sample and
#' as such does not include the S1 sample and as one can see, all densities (but the
#' S1 specific that has only 4 genes) are very balanced around 0.
dir.create("analysis/DESeq2",showWarnings = FALSE)
write.csv(r.vs.w.W[order(r.vs.w.W$padj),],file = "analysis/DESeq2/ray-cells_vs_whole-section_noS1ray_DE-results.csv")
write.csv(r.vs.x.W[order(r.vs.x.W$padj),],file = "analysis/DESeq2/ray-cells_vs_xylem-tracheids_noS1ray_DE-results.csv")
write.csv(x.vs.w.W[order(x.vs.w.W$padj),],file = "analysis/DESeq2/xylem-tracheids_vs_whole-section_noS1ray_DE-results.csv")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
#' 