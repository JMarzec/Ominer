##########################################################################################
#
#File name:generate_venn.R
#
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#
#
#
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
##########################################################################################
###########################################################################################
#####Input arguments:
#####ebA - this is an object generated from limma/edgeR containing information regarding the number of genes that are up/down regulated in similar biological groups
#####project - this is the name given to the analysis
#####This code takes in the ebA object generated from limma/edgeR analysis 
#####Output generated is a Venn diagram showing overlaps between lists of differentially expressed genes between different biological groups
generate_venn <- function (ebA,project) 
{
	readdir <- paste("ominer_results",project,sep="/")
    ebA.decideTests <- read.table(paste(readdir,"/","DifferentialExpression/",ebA,sep=""),sep="\t", as.is = T, 
        header = T, strip.white = T)
    comps = colnames(ebA.decideTests)
    source("Venn.R")
    png(paste("ominer_results/",project,"/","cluster","/vennDiagram.png",sep=""))
    par(mar = c(5, 3, 1, 1))
    a <- vennCounts(ebA.decideTests[, 1:length(comps)], include = "both")
    vennDiagram(a, main = "Results for all probes")
    dev.off()
    png(paste("ominer_results/",project,"/","cluster","/vennDiagram_up.png",sep=""))
    par(mar = c(5, 3, 1, 1))
    a <- vennCounts(ebA.decideTests[, 1:length(comps)], include = "up")
    vennDiagram(a, main = "Results for up-regulated probes")
    dev.off()
    png(paste("ominer_results/",project,"/","cluster","/","/vennDiagram_down.png",sep=""))
    par(mar = c(5, 3, 1, 1))
    a <- vennCounts(ebA.decideTests[, 1:length(comps)], include = "down")
    vennDiagram(a, main = "Results for down-regulated probes")
    dev.off()
    pdf(paste("ominer_results/",project,"/","cluster","/vennDiagram.pdf",sep=""))
    a <- vennCounts(ebA.decideTests[, 1:length(comps)])
    vennDiagram(a, main = "Results for all probes")
    par(mar = c(5, 3, 1, 1))
    a <- vennCounts(ebA.decideTests[, 1:length(comps)], include = "up")
    vennDiagram(a, main = "Results for up-regulated probes")
    par(mar = c(5, 3, 1, 1))
    a <- vennCounts(ebA.decideTests[, 1:length(comps)], include = "down")
    vennDiagram(a, main = "Results for down-regulated probes")
    dev.off()
}
