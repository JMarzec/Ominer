#########################
#
#File mergearoma.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: R script to merge the log2ratios from two files to one for 500K and 100K chips
#(1) cghwebinput - full directory to where log2ratios are(2) targets .txt file with the full path to the target file for the normal samples and the tumor samples
#################################################################################################################

mergeAromaFiles <-function(cghwebinput,targets) {
	setwd(cghwebinput)
	for (i in 1:length(patients))
	{
		
		targs <- read.table(targets,sep="\t",as.is=TRUE)
	targets=targs[1,] 	
	pd <- read.table(targets, sep="\t", as.is=TRUE,strip.white=TRUE,header=T);
	patients <- pd[,2]
		chip1=pd[,1]
		chip2=pd[,2]
	
		chip2[i] = sub(".CEL","",chip2[i])
		chip1[i] = sub(".CEL","",chip1[i])
	
		hind=paste(chip2[i],"_",colnames(pd[2]),".txt",sep="")
		xba=paste(chip1[i],"_",colnames(pd[1]),".txt",sep="")
                print(hind)
 		print(xba)	
		H.data = read.table(hind, sep="\t", header=T,as.is=TRUE)
		X.data = read.table(xba, sep="\t", header=T, as.is=TRUE)
	
		HX.data <- rbind(H.data,X.data)
#remove the col names
	
		HX.data = subset(HX.data, ProbeID != 'ProbeID')	
		HX.data = subset(HX.data, Chromosome != 'Chromosome')	
		HX.data = subset(HX.data, Position != 'Position')	
		HX.data = subset(HX.data, LogRatio != 'LogRatio')	
#duplicated probes
		HX.data <- subset(HX.data, !duplicated(HX.data[,"ProbeID"]))
#duplicated genomic positions
		HX.data <- subset(HX.data, !duplicated(paste(HX.data[,"Chromosome"],HX.data[,"Position"])))
#run CGHweb
		HX.data$Position = as.numeric(HX.data$Position)
		HX.data$Chromosome = as.numeric(HX.data$Chromosome)
		HX.data$LogRatio = as.numeric(HX.data$LogRatio)
		HX.data$ProbeID = sub("SNP_A-","",HX.data$ProbeID)
                print(patients[i])
		write.table(HX.data,file=paste(patients[i],".txt",sep=""),sep="\t",row.names=FALSE)
#		unlink(hind)
#		unlink(xba)
	}
}	
	