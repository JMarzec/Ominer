############################################################################################################
#
#File name threshold_calc_function.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: R script to calculate the threshold for which to call gains and losses
####Input argumnets are project
####thresholds for gain and loss are written to a .txt file
thresholdcalcfunction <- function(Cancer) {
	
path_to_change <- paste("ominer_results/",Cancer,"/cghweb/Matrix/",sep="")
	setwd(path_to_change)
	summary <- "R_input.txt"	
	
    first.data <- read.table(summary, header = T, sep = "\t", 
        as.is = T)
        if(ncol(first.data)>4) {
	   k = stack(first.data[,4:ncol(first.data)])
	   rowmed = k$values
	   }else {
	    rowmed=first.data[,4]
	   }

       pdf("density.pdf")
    plot(density(rowmed, na.rm = T, adjust = 0.5))
        dev.off()
    pdf("histogram.pdf")
    hist(rowmed, freq = F, density = 10, angle = 45, col = "blue")
   
    dev.off()
    gain = quantile(rowmed, probs = 0.90,na.rm=TRUE)
    loss = quantile(rowmed, probs = 0.10,na.rm=TRUE)
    amp = quantile(rowmed, probs = 0.95,na.rm=TRUE)
    del = quantile(rowmed, probs = 0.05,na.rm=TRUE)
   
   	
	write.table(data.frame(gain, loss), file = "threshold.txt", 
        sep = "\t", quote = FALSE, row.names = FALSE)
	
        file.copy("threshold.txt",paste("../../output/threshold.txt",sep=""))
    return(c(gain, loss))
}
