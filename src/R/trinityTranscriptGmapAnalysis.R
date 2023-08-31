#' ---
#' title: "Trinity transcripts GMAP alignment analysis"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/facility/trinity/LMD-trinity/GMAP")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/facility/trinity/LMD-trinity/GMAP")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(genomeIntervals))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(VennDiagram))

#' Source helpers
source("~/Git/UPSCb/src/R/gff3Utilities.R")
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' Graphic params
mar <- par("mar")
pal <- brewer.pal(8,"Dark2")

#' The gene annotation
annot <- read.delim("/mnt/picea/projects/spruce/genome/v1.0/annotation/transcriptome/january2013/augustus-eugene-january2013.tsv",as.is=TRUE)

#' The Trinity fasta file
trinity.seq <- readDNAStringSet("../TransDecoder/Trinity.fasta.transdecoder.cds")
names(trinity.seq) <- sub(" .*","",names(trinity.seq))

#' # Analysis
#' Reading in all the gmap gff3 files
gff3s <- mclapply(dir(pattern="Pabies1.0-Trinity.2\\.[m,t,u][^\\.]+$",
                      full.names = TRUE),
                  readGff3,mc.cores=3L)
names(gff3s) <- gsub("Pabies1.0-Trinity.2\\.|\\.gff3","",dir(pattern="Pabies1.0-Trinity.2\\.[m,t,u][^\\.]+$"))

#' ## GMAP results
#' Create a list of IDs
ID.list <- lapply(mclapply(mclapply(gff3s,function(g){g[g$type=="gene"]},mc.cores=length(gff3s)),
                                    getGffAttribute,"ID",mc.cores=length(gff3s)),as.vector)
names(ID.list) <- sub(".*\\.","",names(gff3s))

Parent.ID.list <- mclapply(gff3s,
                           getGffAttribute,c("ID","Parent"),mc.cores=length(gff3s))
names(Parent.ID.list) <- names(ID.list )

#' * Barplot the paths
barplot(elementNROWS(ID.list),main="GMAP alignment (gene)",xlab="alignment type",ylab="occurrence")

#' * Check the GMAP "genes"
barplot(elementNROWS(sapply(lapply(ID.list,function(IDl){sub("\\.path.*","",IDl)}),unique)),
        main="GMAP alignment (unique gene)",xlab="alignment type",ylab="occurrence")

#' * How many exon per gene? There is one mRNA per gene
sel <- grep("exon[0-9]+$",Parent.ID.list[[1]][,"ID"])
barplot(table(elementNROWS(split(Parent.ID.list[[1]][sel,"ID"],
                                 sub("\\.mrna[0-9]+","",Parent.ID.list[[1]][sel,"Parent"])))),
        xlab="exon per mRNA",las=2,ylab="occurrence",main="GMAP alignment number of exons")

#' * Overlap between the alignment types
#' 
#' Expected, but nonetheless :-)
grid.newpage()
grid.draw(venn.diagram(lapply(ID.list,sub,pattern="\\.path[0-9]+$",replacement=""),filename = NULL,fill=pal[1:3]))

#' * Quantify the paths
#' 
#' The transloc involve only 2 parts (possibly a limitation from gmap)
#' 
#' The multiple are mostly found across 2 genes, but 5 (max setting of the gmap run)
#' still does not exhaust the possibilities!
pander(lapply(ID.list,function(ls){table(sub(".*.path","",ls))}))

#' * Get the sub gffs
sgff3 <- mcmapply(extractFromGff3UsingGeneIDs,gff3s,ID.list,mc.cores=3L)

cov.idt <- mclapply(sgff3,function(gff){
  t(apply(getGffAttribute(gff[gff$type=="mRNA",],c("coverage","identity")),1,as.numeric))
},mc.cores=3L)

names(cov.idt) <- names(ID.list)

#' The unique alignment as well as the fragmented ones have a very good 
#' identity distribution. The unique ones can go as down as 80% identity,
#' which may indicate that we may have either wrongly sequenced/assembled reads, 
#' or the converse, badly assembled genomic loci. The identity can get far worse
#' for the translocation which is probably more indicative of genomic assembly 
#' issues (e.g. collapsed gene domains or collapsed gene and pseudogene...)
#' 
#' The multiple have a more variable distribution, which spans however mostly the
#' same range as the unique one, indicative that we may be looking at a real
#' gene family. The peak is still at a hundred percent, so there may be recent
#' duplication or haplotypic sequences in the genome.
plot.multidensity(lapply(cov.idt,function(m){m[,2]}),
                  xlab="identity",lwd=2,
                  legend = names(cov.idt),
                  legend.lwd = 2)

#' * Boxplot
#' Most of the unique path have a full coverage
#' 
#' The transloc are as expected around 50%
#' 
#' The multiple seem to also contain a majority of many-to-many translocations
boxplot(lapply(cov.idt,function(m){m[,1]}),ylab="coverage",
        main="Coverage by alignment type")

#' * Cumulative coverage of transloc
#' 
#' It nicely sum up to a 100% for the vast majority
plot.multidensity(list(cumulative=sapply(split(cov.idt[["transloc"]][,1],
                                               sub("\\.path[0-9]+$","",ID.list[["transloc"]])),sum),
                       individual=cov.idt[["transloc"]][,1]),xlab="coverage",
                  main="GMAP transloc coverage")

#' * Cumulative coverage of the mult
#' 
#' They nicely sum up to multiple copies (2x to 5X)
d.list <- list(cumulative=sapply(split(cov.idt[["mult"]][,1],
                                       sub("\\.path[0-9]+$","",ID.list[["mult"]])),sum),
               individual=cov.idt[["mult"]][,1])
plot.multidensity(d.list,xlab="coverage",
                  main="GMAP multiple coverage")

boxplot(d.list)
tab <- table(d.list[["cumulative"]])
barplot(tab)
pander(head(sort(tab,decreasing = TRUE)))

#' * Full length estimation
#' 
#' If we consider anything over 95% as full length based on the transloc
sel <- cov.idt[["mult"]][,1] >= 95

#' The full length multiple
sprintf("There are %s sequences that are full length",length(unique(sub("\\.path[0-9]+$","",ID.list[["mult"]][sel]))))

#' Most of the sequence are present twice. The one third that is present only once
#' is probably indicative of a full sequence, whose multiple hits are fragmented
plot.multidensity(list(cumulative=sapply(split(cov.idt[["mult"]][sel,1],
                                               sub("\\.path[0-9]+$","",ID.list[["mult"]][sel])),sum)),
                  xlab="coverage",
                  main="GMAP multiple coverage (FL only)")

#' The fragmented multiple
sprintf("There are %s sequences that are fragmented",length(unique(sub("\\.path[0-9]+$","",ID.list[["mult"]][!sel]))))

#' This is obviously messier
plot.multidensity(list(cumulative=sapply(split(cov.idt[["mult"]][!sel,1],
                                               sub("\\.path[0-9]+$","",ID.list[["mult"]][!sel])),sum)),
                  xlab="coverage",
                  main="GMAP multiple coverage (no FL)")

#' ## BEDtools results
#' This is the intersect of the GMAP results with the v1.0 spruce genome 
#' Eugene annotation (the gene only for both)
beds <- mclapply(dir("../bedtools",pattern="Pabies1.0-Trinity.2.*\\.tsv",full.names = TRUE),
                 read.delim,as.is=TRUE,header=FALSE,mc.cores=length(ID.list))

names(beds) <- gsub(".*\\.gff3\\.|-Eugene.*","",
                    dir("../bedtools",pattern="Pabies1.0-Trinity.2.*\\.tsv"))

#' Find the annotated gene IDs that are overlapped
gIDs <- lapply(lapply(beds,"[[",18),sub,pattern="ID=",replacement="")

#' * Barplot
#' 
#' The amount of genes overlaped by the GMAP alignments is surprisingly high for
#' for the multiples. We see the expected 30% fragmentation of gene loci (ratio
#' transloc over uniq)
barplot(elementNROWS(lapply(gIDs,unique)))

#' * Venn Diagram
#' 
#' Surprising to see so much overlap between the different alignment classes
grid.newpage()
grid.draw(venn.diagram(lapply(gIDs,unique),filename = NULL,fill=pal[1:3]))

#' How many genes
message(sprintf("There are %s overlapped genes",length(unique(unlist(gIDs)))))

#' * Confidence
#' 
#' Unsurprisingly, the uniq are mostly HC, whereas the mult and transloc are 
#' half-half LC and MC. The multiple also have the majority of LCs, which may
#' represent assembly artifacts (collapsed conserved regions of large gene families).
#' 
barplot(sapply(lapply(lapply(covs,"[[","ID"),unique),
               function(id){table(annot[annot$ID %in% id,"confidence"])}),
        beside = TRUE, col=pal[1:3],legend.text = c("HC","LC","MC"))

#' * Coverage
#' 
#' Let's assess those alignments that fully cover annotated genes
covs <- mclapply(beds,function(b){data.frame(gID=gsub("ID=|\\.path[0-9]+;.*","",b[,9]),
                                           gLen=b[,5] - b[,4] + 1,
                                           ID=sub("ID=","",b[,18]),
                                           Len=b[,14] - b[,13] + 1,
                                           stringsAsFactors = FALSE)},mc.cores=length(gIDs))

names(covs) <- names(gIDs)

#' Number of reads mapping a gene
#' 
#' Min overlap of 30%
nams <- unique(unlist(sapply(ID.list,sub,pattern="\\.path[0-9]+",replacement="")))
tab <- table(nams %in% unique(unlist(sapply(covs,"[","gID"))))
pander(tab)
sprintf("There are %s transcripts not mapping an annotated gene", round(tab["FALSE"]/sum(tab) * 100,digits=2))

#' # Expression
library(tximport)
tx <- tximport(list.files("../Salmon",pattern="*.sf",recursive=TRUE,full.names=TRUE),
         type="salmon",txOut=TRUE)
        
counts <- sapply(split.data.frame(
  t(as.data.frame(tx$counts)),
  f=factor(sub("_L00[1-4].*","",
               dirname(list.files("../Salmon",
                                  pattern="*.sf",
                                  recursive=TRUE))))),colSums)

nms <- nams[nams %in% rownames(counts)]

boxplot(log2(as.vector(counts[nms[!nms %in% unique(unlist(sapply(covs,"[","gID")))],])),
        log2(as.vector(counts[nms[nms %in% unique(unlist(sapply(covs,"[","gID")))],])))

source("~/Git/UPSCb/src/R/percentile.R")
percentile(log2(as.vector(counts[nms[!nms %in% unique(unlist(sapply(covs,"[","gID")))],])))

#' Number of transcripts not aligning
tab  <- table(names(trinity.seq) %in% unique(unlist(sapply(ID.list,sub,pattern="\\.path[0-9]+",replacement=""))))
pander(tab)
sprintf("There are %s transcripts not mapping in the genome", round(tab["FALSE"]/sum(tab) * 100,digits=2))


#' * GMAP hits coverage onto annotated genes
#' 
#' A majority of the unique alignments have a coverage slightly > 1, 
#' so we can expect to even have some UTRs
#' 
#' This is less pronounced for the multiple and transloc. Interestingly for both
#' they have a huge pick showing very little overlap with an existing gene. 
#' These could indicate an amount of assembly or annotation issues
#' 
boxplot(lapply(covs,function(cov){cov$gLen/cov$Len}),ylab="fraction covered",
        main="GMAP hit length / Annotation length")
abline(h=1,col="grey",lty=2)
plot.multidensity(lapply(covs,function(cov){cov$gLen/cov$Len}),
                  xlab="fraction covered",legend.lwd = 2,
                  main="GMAP hit length / Annotation length")
abline(v=1,col="grey",lty=2)

#' * Alignments covering at least 90% of the genes
FL.IDs <- mapply("[",gIDs,lapply(covs,function(cov){(cov$gLen/cov$Len)>0.9}))


#' It surprising that so many different alignment type overlap the same genes 
#' as "full-length". In particular for the "uniq", where this would rather indicate 
#' that these transcripts are only partial to either end of the split.
grid.newpage()
grid.draw(venn.diagram(FL.IDs,filename = NULL,fill=pal[1:3]))

#' * Full length coverage assessment
#' 
#' If we consider only the subset of genes unique to the uniq type, then about 
#' 80% are covered by full-length reads.
table(unique(FL.IDs$uniq) %in% setdiff(FL.IDs[["uniq"]],union(FL.IDs[["transloc"]],FL.IDs[["mult"]])))

#' Coming from about 2.5X more transcripts
table(FL.IDs$uniq %in% setdiff(FL.IDs[["uniq"]],union(FL.IDs[["transloc"]],FL.IDs[["mult"]])))

#' Which represent about 1/3 of all the uniq transcripts
sum(FL.IDs$uniq %in% setdiff(FL.IDs[["uniq"]],union(FL.IDs[["transloc"]],FL.IDs[["mult"]])))/length(gIDs$uniq)

#' write a DF
df <- cbind(data.frame(ID=names(trinity.seq),
                       length=width(trinity.seq),
                       GC=rowSums(alphabetFrequency(trinity.seq,as.prob = TRUE)[,c("C","G")]),
                 GMAP=names(trinity.seq) %in% unique(unlist(sapply(ID.list,sub,pattern="\\.path[0-9]+",replacement="")))),
      as.data.frame(do.call(cbind, mclapply(gff3s,function(g){
        match(names(trinity.seq),unique(sub("\\.path[0-9]+$","",
                                            getGffAttribute(g[g$type=="gene",],"ID"))))},
        mc.cores=3L))),
      data.frame(
        GeneID,
        GeneLength,
        OvlLength,
      ))
                 
#' Overlap gene
#' Confidence
)


test <- lapply(covs,function(cv){
  res <- t(sapply(split.data.frame(cv,cv$gID),function(x){apply(x,2,paste,collapse="-")}))
  res[,"gID"] <- sub("\\|.*","",res[,"gID"])
  res[match(names(trinity.seq),res[,"gID"]),c("ID","Len")]
  })

#' # Conclusion
#' 
#' 
#' ```{r empty, echo=FALSE,eval=FALSE}
#' ```
#' 
#' # Session Info
#' ```{r session info, echo=FALSE}`
#' sessionInfo()
#' ```
#'
