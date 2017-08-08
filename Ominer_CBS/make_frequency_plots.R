############################################################################################################
#
#File make_frequency_plots.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to generate frequency plots from binary_code_data and binary coded filtered data
#Input a text file of log2ratios - columns are sample names, probeids, chromosome and position are also present as column headers.
#Sample input file R_input.txt
#inputs to function are: (1) results = .txt file of binary coded data, (2) target = .txt of target file, (3) Cancer - name of project, (4) type = datatype i.e. CEL and (5) name = name of plot







makeFreqPlot <-function(results,target,Cancer,type,name) {
 
   res <- read.table(results, header = T, sep = "\t", as.is = T)

  targs <- read.table(target, sep = "\t", as.is = TRUE)
    targets = targs[1, ]
    pd <- read.table(targets, sep = "\t", as.is = TRUE, strip.white = TRUE, 
        header = T)
	levels=unique(pd$Group)
	patientcount <- length(pd$Group)
    print(levels)
	par(mar=c(8,4,4,4))
	for (lev in 1:length(levels))
	{

		if(type == "CEL") {
		if(platform=="500" || platform=="100") {
				
			  x =(pd[pd$Name == levels[lev],3])
			  }
			  else {
			   #x = (pd[pd$Name == levels[lev], 2])
			    x = (pd[pd$Group == levels[lev], 2])
			   print(x)
			   
        }
	}
			
		res$Loss=rep(0,length(nrow(res)))
		res$Gain=rep(0,length(nrow(res)))
	
		data.plot=res[,c("ProbeID","Chromosome","Position",x,"Gain","Loss")]
		data.plot.subset = data.plot[order(as.numeric(data.plot$Chromosome),as.numeric(data.plot$Position)),]
	
		if(length(x)>1)
		{
			data.plot.subset[,"Gain"] = rowSums(data.plot.subset[,x ]==1)
			data.plot.subset[,"Loss"] = rowSums(data.plot.subset[,x ]==-1)
		}
	
		if(length(x)==1) {
			data.plot.subset[,"Gain"] = data.plot.subset[,x ]
			data.plot.subset[,"Loss"] = data.plot.subset[,x ]
			data.plot.subset[,"Gain"] = replace(data.plot.subset[,"Gain"], data.plot.subset[,"Gain"] ==-1,0)
			data.plot.subset[,"Loss"] = replace(data.plot.subset[,"Loss"], data.plot.subset[,"Loss"] ==1,0)
			data.plot.subset[,"Loss"] = replace(data.plot.subset[,"Loss"], data.plot.subset[,"Loss"] ==-1,1)
		}
	
	
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
	
		y1=100*(data.plot.subset$Gain/n)
		y2=paste("-",100*(data.plot.subset$Loss/n),sep="")
		
		#png(paste(name,levels[lev],".png",sep=""),width=1024, height=300)
		 png(paste("ominer_results/",Cancer, "/FP/",levels[lev],name,".png", sep = ""), width = 1024, 
            height = 300)

		#plot(y1, type='h',  xaxt="n",  yaxt="n", col="green", main = paste(Cancer,levels[lev],n,sep=" ") , ylim=range(-100, 100), xlab='Chromosome Number',  ylab='Fraction of Patient for Gain or Loss',xaxs = "i", yaxs = "i")
#plot loss

		plot(y1, type='h',  xaxt="n",  yaxt="n", col="green" , ylim=range(-100, 100), xlab='Chromosome Number',  ylab='Fraction of Patient for Gain or Loss',xaxs = "i", yaxs = "i")


		points(y2, type='h', col="red")
#label for chromsomes
		x= c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
		y = F1[,1]
		axis(1, at = c(F1[,2]), label =FALSE, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
		axis(1, at = c(F1[,3]), label =y, tick = FALSE );
		axis(2, at = c(x), labels=x, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
		print(patientcount)
		patientcount=as.numeric(patientcount)
		abline( h = patientcount, col = 'black', lwd = 2)
		patientcount<- patientcount *-1
		abline( h = patientcount, col = 'black', lwd = 2)
		dev.off()
                par(mar=c(8,4,4,4))
		pdf(paste(name,levels[lev],".pdf",sep=""))
		
		#plot(y1, type='h',  xaxt="n",  yaxt="n", col="green", main = paste(Cancer,levels[lev],n,sep=" ") , ylim=range(-100, 100), xlab='Chromosome Number',  ylab='Fraction of Patient for Gain or Loss',xaxs = "i", yaxs = "i")
				plot(y1, type='h',  xaxt="n",  yaxt="n", col="green", ylim=range(-100, 100), xlab='Chromosome Number',  ylab='Fraction of Patient for Gain or Loss',xaxs = "i", yaxs = "i")

#plot loss
		points(y2, type='h', col="red")
#label for chromsomes
		x= c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
		y = F1[,1]
		axis(1, at = c(F1[,2]), label =FALSE, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
		axis(1, at = c(F1[,3]), label =y, tick = FALSE );
		axis(2, at = c(x), labels=x, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
		abline( h = patientcount, col = 'black', lwd = 2)
		patientcount<- patientcount *-1
		abline( h = patientcount, col = 'black', lwd = 2)
		dev.off()
		
		
		 pdf(paste("ominer_results/",Cancer, "/FP/",levels[lev],name,".pdf", sep = ""), width = 1024, 
            height = 300)

		plot(y1, type='h',  xaxt="n",  yaxt="n", col="green", main = paste(Cancer,levels[lev],n,sep=" ") , ylim=range(-100, 100), xlab='Chromosome Number',  ylab='Fraction of Patient for Gain or Loss',xaxs = "i", yaxs = "i")
#plot loss
		points(y2, type='h', col="red")
#label for chromsomes
		x= c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
		y = F1[,1]
		axis(1, at = c(F1[,2]), label =FALSE, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
		axis(1, at = c(F1[,3]), label =y, tick = FALSE );
		axis(2, at = c(x), labels=x, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
		print(patientcount)
		patientcount=as.numeric(patientcount)
		abline( h = patientcount, col = 'black', lwd = 2)
		patientcount<- patientcount *-1
		abline( h = patientcount, col = 'black', lwd = 2)
		dev.off()
                par(mar=c(8,4,4,4))
		pdf(paste(name,levels[lev],".pdf",sep=""))
		
		plot(y1, type='h',  xaxt="n",  yaxt="n", col="green", main = paste(Cancer,levels[lev],n,sep=" ") , ylim=range(-100, 100), xlab='Chromosome Number',  ylab='Fraction of Patient for Gain or Loss',xaxs = "i", yaxs = "i")
#plot loss
		points(y2, type='h', col="red")
#label for chromsomes
		x= c(-100,-80,-60,-40,-20,0,20,40,60,80,100)
		y = F1[,1]
		axis(1, at = c(F1[,2]), label =FALSE, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
		axis(1, at = c(F1[,3]), label =y, tick = FALSE );
		axis(2, at = c(x), labels=x, tick = TRUE, las = 1, col = "black", lty = "dotted", tck = 1 );
		abline( h = patientcount, col = 'black', lwd = 2)
		patientcount<- patientcount *-1
		abline( h = patientcount, col = 'black', lwd = 2)
		dev.off()
		
	

		
	}
	
}
