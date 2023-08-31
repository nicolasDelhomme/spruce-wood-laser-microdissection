#' ---
#' title: "Spruce LMD pathway analysis"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Environment
#' Set the working dir
setwd("/mnt/picea/projects/spruce/akaerkoenen/laser-capture-microdissection-olga/analysis")

#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/akaerkoenen/laser-capture-microdissection-olga/analysis")
#' ```

#' Libraries
suppressPackageStartupMessages(library(KEGGREST))
suppressPackageStartupMessages(library(pathview))
suppressPackageStartupMessages(library(S4Vectors))

#' Read the expression data
exp <- read.csv("HTSeq/library-size-normalized_variance-stabilized_merged-technical-replicates_noS1ray_data.csv",
                row.names=1)

#' Get the KEGG info
map <- read.delim("/mnt/picea/projects/gopher2/pabies/gene_to_kegg.tsv",header=FALSE,
                  stringsAsFactors = FALSE)
revmap <- strsplit(map[,2],"\\|")

revmap <- cbind(unlist(revmap),rep(map[,1],elementNROWS(revmap)))

#' # Data extraction
#' Get the KEGG PATHWAY ko00940 (Phenylpropanoid biosynthesis)
pathway <- keggGet("ko00940")

#' Extract the enzymes
enzymes <- strsplit(gsub(".*\\[|\\]","",pathway[[1]]$ORTHOLOGY)," ")
enzymes <- cbind(sub("EC:","",unlist(enzymes)),rep(names(enzymes),elementNROWS(enzymes)))

revmap[revmap[,1] %in% sub("\\.-","",enzymes[,1]),]

#' # Data plot

#' # Session Info
#' ```{r, session info, echo=FALSE}
#' sessionInfo()
#' ``` 