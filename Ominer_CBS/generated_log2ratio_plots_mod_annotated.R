############################################################################################################
#
#generated_log2ratio_plots_mod.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to generate frequency plots from log2ratio data
###Arguments - (1) data - log2ratio data, (2) target - .txt file containing the full path to target files for tumor and normal samples, (3) name to give the new plots, (4) 
####type of data e.g. CEL, log2raio, smoothed or segmented 
log2ratioPlots <-function(data,target,name,type) {
	  #name = "log2ratioPlot_"
	  targs <- read.table(target, sep = "\t", as.is = TRUE)
    targets = targs[1, ]
    print (targets)
    pd <- read.table(targets, sep = "\t", as.is = TRUE, strip.white = TRUE, 
        header = T)
        
        	if (type == "CEL") {
		if(platform=="500" || platform=="100") {
		patients=pd[,3]
		}
		else {
		patients <- pd[,2]
		}
        	}
        	else {
        		patients <- pd[,1]
        	}

                      
res <- read.table(data, header = T, sep = "\t", 
        as.is = T)	
     	
	for(x in 1:length(patients))
	{
# patient ID
		print(patients[x])
# res needs to be the R_input.txt file output not binary coded data
		
		res$Loss=rep(0,length(nrow(res)))
		res$Gain=rep(0,length(nrow(res)))
		
		data.plot=res[,c("ProbeID","Chromosome","Position",patients[x],"Gain","Loss")]
		data.plot.subset = data.plot[order(as.numeric(data.plot$Chromosome),as.numeric(data.plot$Position)),]
		
		data.plot.subset[,"Gain"] = data.plot.subset[,patients[x] ]
		data.plot.subset[,"Loss"] = data.plot.subset[,patients[x] ]
		data.plot.subset[,"Gain"] = replace(data.plot.subset[,"Gain"], data.plot.subset[,"Gain"] <0,0)
		data.plot.subset[,"Loss"] = replace(data.plot.subset[,"Loss"], data.plot.subset[,"Loss"] >0,0)
		
#label chromosomes
		labels_chr <- data.matrix(summary(as.factor(data.plot.subset$Chromosome)))
		
		test1<- data.frame(labels_chr,row.names(labels_chr) )
		test <- data.frame(unique(data.plot.subset$Chromosome))
		colnames(test) = c("Chromosome")
		colnames(test1) = c("count","Chromosome")
		F1 <- merge(test,test1, by.x="Chromosome", by.y="Chromosome", sort=FALSE)
		for(i in 2:length(row.names(F1)))
		{F1[i,2] = F1[i-1,2] + F1[i,2] ; }
		
		F1$label <- NULL ; F1[1,3] <- F1[1,2] / 2 ;
		for (i in 2:length(row.names(F1))){ F1[i,3] <- (F1[i,2]+F1[i-1,2])/2; }
		n=length(x)
		
		y1=data.plot.subset$Gain
		y2=data.plot.subset$Loss

		png(paste("ominer_results/",Cancer,"/FP/",name,patients[x],".png",sep=""),width=1024, height=300)
		plot(y1, type='h',  xaxt="n",  yaxt="n", col="green",   xlab='Chromosome Number',  xaxs = "i", yaxs = "i",ylim=range(min(y2-0.5,na.rm=TRUE), max(y1+0.5,na.rm=TRUE)), main = paste(name,patients[x],sep=" "),ylab='Log2Ratio',xaxs = "i", yaxs = "i")
#plot loss
		points(y2, type='h', col="red")
#label for chromsomes
		y = F1[,1]
		axis(1, at = c(F1[,2]), label =FALSE, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
		axis(1, at = c(F1[,3]), label =y, tick = FALSE );
		#axis(2, at = c(x), labels=x, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
		dev.off()
		
		pdf(paste("ominer_results/",Cancer,"/FP/",name,patients[x],".pdf",sep=""),width=1024, height=300)
		plot(y1, type='h',  xaxt="n",  yaxt="n", col="green",   xlab='Chromosome Number',  xaxs = "i", yaxs = "i",ylim=range(min(y2-0.5,na.rm=TRUE), max(y1+0.5,na.rm=TRUE)), main = paste(name,patients[x],sep=" "),ylab='Log2Ratio',xaxs = "i", yaxs = "i")
#plot loss
		points(y2, type='h', col="red")
#label for chromsomes
		y = F1[,1]
		axis(1, at = c(F1[,2]), label =FALSE, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
		axis(1, at = c(F1[,3]), label =y, tick = FALSE );
		#axis(2, at = c(x), labels=x, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
		dev.off()

	}
}
