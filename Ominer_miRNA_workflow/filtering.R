##########################################################################################
#
#File name: filtering.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#
#
#
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
##########################################################################################
##########################################################################################
#Description: function to filter the normalised expression matrix using a filter method and a filtervalue
#filter - method used to filter the normalised expression matrix e.g. sd, igr etc
#filterval - value used to filter the normalised expression matrix: for exampe to retail the top 40% of most variable probes filterval = 40
#combat - 1 if COMBAT algorithm is to be applied and 0 if not
#project - name given to the analysis
################################################################################################KM plots - overall survival from a target file generates .txt file and KM plots for 5,10 and 15 years
######
library("affy")
library("affyPLM")
library("simpleaffy")

filtering <- function (filter, filterval,combat,project) {
if (combat=="1") {
    Normdata <- read.table(paste("ominer_results/",project,"/","norm","/","data_combat.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1)
        }
        else {
         Normdata <- read.table(paste("ominer_results/",project,"/","norm","/","normalised.txt",sep=""), sep = "\t", as.is = T, header = T, 
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
    print ("I AM HERE NOWWWWW BEFORE PRINTING THIS OUT")
    write.table(A.data, paste("ominer_results/",project,"/","norm","/","filtered_data.txt",sep=""), sep = "\t", quote = FALSE, 
        row.names = TRUE, col.names = NA)
    return(A.data)
    print ("I HAVE PRINTED OUT THIS DATA")
}
