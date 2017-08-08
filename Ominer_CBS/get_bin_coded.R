############################################################################################################
#
#File name: get_bin_coded.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: R script to recode the SNP positions into gains (1) and losses(0)
#####output is a binary coded matrix with samplenames as columns 
#####Arguments - project name is Cancer 
getBinaryCoded <- function(Cancer) {
	
path_to_change <- paste("ominer_results/",Cancer,"/output/",sep="")	
setwd(path_to_change)
    thresh <- read.table("threshold.txt", header = T, sep = "\t", 
        as.is = T)
    gain = thresh[, 1]
    loss = thresh[, 2]
    data <- read.table("R_input.txt", header = T, sep = "\t", as.is = T)
    first.data = data
    results = first.data
    a = colnames(first.data[, 4:ncol(first.data)])
    for (j in a) {
        results[(results[, j] > 0 & !is.na(results[, j]) & results[, 
            j] >= gain), j] = 1
        results[(results[, j] < 0 & !is.na(results[, j]) & results[, 
            j] <= loss), j] = -1
        results[(results[, j] > 0 & !is.na(results[, j]) & results[, 
            j] < gain), j] = 0
        results[(results[, j] < 0 & !is.na(results[, j]) & results[, 
            j] > loss), j] = 0
            results[is.na(results[,j]),j]=0; ###added

    }
    write.table(results, file = "results.txt", sep = "\t", quote = FALSE, 
        row.names = FALSE)
}
