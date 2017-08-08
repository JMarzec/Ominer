##########################################################################################
#
#File name: generate_ea_venn.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#
#
#
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
##########################################################################################
##########################################################################################
#Description: function to generate a venn diagram for up to four different comparisons from the results of differentially expressed genes from transcripts
#ebA.txt - is a file generated from running limma analysis giving the number(s) of genes up and down regulated both between and within the different biological groups
#project is the given project name
################################################################################################
######

generate_venn <- function (ebA,project) 
{
	readdir <- paste("ominer_results",project,sep="/")
    ebA.decideTests <- read.table(paste(readdir,"/","transcript/DifferentialExpression/",ebA,sep=""),sep="\t", as.is = T, 
        header = T, strip.white = T)
    comps = colnames(ebA.decideTests)
    source("Venn.R")
    png(paste("ominer_results/",project,"/","transcript/cluster","/vennDiagram.png",sep=""))
    par(mar = c(5, 3, 1, 1))
    a <- vennCounts(ebA.decideTests[, 1:length(comps)], include = "both")
    vennDiagram(a, main = "Results for all probes")
    dev.off()
    png(paste("ominer_results/",project,"/","transcript/cluster","/vennDiagram_up.png",sep=""))
    par(mar = c(5, 3, 1, 1))
    a <- vennCounts(ebA.decideTests[, 1:length(comps)], include = "up")
    vennDiagram(a, main = "Results for up-regulated probes")
    dev.off()
    png(paste("ominer_results/",project,"/","transcript/cluster","/","/vennDiagram_down.png",sep=""))
    par(mar = c(5, 3, 1, 1))
    a <- vennCounts(ebA.decideTests[, 1:length(comps)], include = "down")
    vennDiagram(a, main = "Results for down-regulated probes")
    dev.off()
    pdf(paste("ominer_results/",project,"/","transcript/cluster","/vennDiagram.pdf",sep=""))
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
