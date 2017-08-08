#########################
#
#get_binary_filtered.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: R script to create regions by filtering the segmeneted data according to the minmum number of consecutive SNPs
#####Input argumnet - number to use for the minimum number of SNPs
#.txt file of binary coded regions - 0 for  loss and 1 for gain 
#################################################################################################################
filterBinary <- function (snp_number) 
{
  path_to_change <- paste("ominer_results/",Cancer,"/output/",sep="")	
setwd(path_to_change)
  
    results <- read.table("results.txt", header = T, sep = "\t", as.is = T)
    print("filter binary")
    print(snp_number)
    results.filtered = results
    a = colnames(results[, 4:ncol(results)])
    for (j in a) {
        x = results[, j]
        	#y=results[,j]
        z.rle <- rle(x)
        ends <- cumsum(z.rle$lengths)
        starts <- ends - z.rle$lengths + 1
        indexes <- with(z.rle, data.frame(starts, ends, lengths, 
            values))
        s = starts[z.rle$lengths < snp_number & z.rle$values != 0]
        e = ends[z.rle$lengths < snp_number & z.rle$values != 0]
        for (i in 1:length(s)) {
            x[s[i]:e[i]] = 0
            #y[s[i]:e[i]]=0
        }
        results.filtered[, j] = x
        	#first.data[,j]=y
    }
    write.table(results.filtered, "binary_coded_filtered.txt", 
        sep = "\t", row.names = FALSE)
}
