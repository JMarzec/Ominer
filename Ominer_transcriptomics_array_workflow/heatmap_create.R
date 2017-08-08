##########################################################################################
#
#File name: heatmap_create.R
#
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#
#
#
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
##########################################################################################
##########################################################################################
#Description: function to generate a heatmap from the list(s) of teh differentially expressed genes within a dataset: uses the filtered normalised expression matrix and list of differentially expressed genes as input 
#Arguments are:
#target_file = path to the target file 
#project is the given project name
#comp_file = path to the comparisons file 
################################################################################################KM plots - overall survival from a target file generates .txt file and KM plots for 5,10 and 15 years
######
heatmap_generate <- function (target, project,comp_file) 
{
    library("gplots")
    library("limma")
    library("hgu133plus2.db")
    readdir <- paste("ominer_results",project,sep="/")
    A.data <- read.table(paste(readdir,"/","norm","/filtered_data.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1)
    pd <- read.table(target, sep = "\t", as.is = TRUE, strip.white = TRUE, 
        header = T)
    colnames(A.data) = pd$Target
    
    comps <- read.table(comp_file, sep = "\t", as.is = T, header = F, 
        strip.white = T)
        for (i in 1:length(comps)) {
                print(comps[i])
        fg <- strsplit(comps[,i],"=")
        ci <- fg[[i]][1]
####need to do this for all comparisons
        #filteredx <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_annotated_filtered_DE.txt",sep=""),row.names = 1, header = T, sep = "\t") 
#filteredx <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_annotated_filtered_DE.txt",sep=""), header = T, sep = "\t",as.is=T) 
x <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_Filtered.txt",sep=""), header = T, sep = "\t",as.is=T) 

    
    #x <- read.table(paste(readdir,"/","DifferentialExpression/",de_genes,sep=""), blank.lines.skip = FALSE, sep = "\t", 
        #as.is = TRUE, header = TRUE)
    v = x$probeid
    esetSel = A.data[v, ]
    #symbol <- as.character(unlist(lapply(mget(v, env = get(paste(platform, 
     #   "SYMBOL", sep = ""))), function(symbol) {
      #  return(paste(symbol, collapse = "; "))
    #})))
    es = as.matrix(esetSel)
    #es = as.numeric(esetSel)
    colnames(es) = pd$Target
    png(paste(readdir, "/cluster/",ci,".heatmap.png",sep=""))
    res=900
    heatmap.2(es, col = redgreen(75), scale = "row", key = "TRUE", 
        symkey = "FALSE", density.info = "none", trace = "none", 
        cexRow = 0.4, cexCol = 0.5, labRow = v)
    dev.off()
}
}
