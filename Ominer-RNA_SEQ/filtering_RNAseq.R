##########################################################################################
#
#File name:filtering_RNAseq.R
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
#####filter - filtering method to use for the normalised expression matrix - for example 
#####filterval - value used to filter the data for example to keep the top 40% of most variable probes use filterval = 40
#####combat = 0/1 - o if the COMBAT algorithm is not to be used and i if the algorithm is to be run
#####project - this is the name given to the analysis
#####data - this refers to the filtered normalised expression matrix i.e. filtered_data.txt
#####This code filters the normalised expression matrix and keeps the top x% of most variable probes

library("affy")
library("affyPLM")
library("simpleaffy")

filtering <- function (filter, filterval,combat,project,data) {
if (combat=="1") {
    Normdata <- read.table(paste("ominer_results/",project,"/","norm","/","data_combat.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1)
        }
        else {
         Normdata <- read.table(paste("ominer_results/",project,"/","norm","/",data,sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1)
        }
    Norm <- new("ExpressionSet", exprs = data.matrix(Normdata))
    Norm.exps = exprs(Norm)
    dim(Norm.exps)
    Norm.exps = Norm.exps[substring(row.names(Norm.exps), 1, 
        4) != "AFFX", ]
    filterval = as.numeric(filterval)
    if (filter == "iqr") {
        library(genefilter)
        f2 <- function(x) (IQR(x) > filterval)
        ff <- filterfun(f2)
        selected <- genefilter(Norm.exps, ff)
        sum(selected)
        A.data <- Norm.exps[selected, ]
    }
    if (filter == "sd") {
        rsd <- apply((Norm.exps), 1, sd)
        sel <- order(rsd, decreasing = TRUE)[1:(nrow(Norm.exps) * 
            filterval/100)]
        A.data <- (Norm.exps[sel, ])
    }
    if (filter == "intensity") {
        library(genefilter)
        f2 <- pOverA(filterval, log2(100))
        ff <- filterfun(f2)
        selected <- genefilter(Norm.exps, ff)
        A.data <- Norm.exps[selected, ]
        sum(selected)
    }
    write.table(A.data, paste("ominer_results/",project,"/","norm","/","filtered_data.txt",sep=""), sep = "\t", quote = FALSE, 
        row.names = TRUE, col.names = NA)
    return(A.data)
}
