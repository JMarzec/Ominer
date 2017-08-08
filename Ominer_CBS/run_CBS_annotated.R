
############################################################################################################
#
#File name runCBS.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: R script to segment the log2ratios - uses inly one segmentation model CBS
#####output is a matrix containing the segmented values for each of the samples 
runCGH <- function(outputdir,targets,Cancer,type,platform) {
library(CGHweb)
targs<-read.table(targets,sep="\t",as.is=TRUE)
 targets = targs[1, ]
 current_location <- getwd()
 
    pd <- read.table(targets, sep = "\t", as.is = TRUE, strip.white = TRUE, 
        header = T)
                if (platform == "100") {
        	patients <- pd[,3]
        }
        else {
              patients <- pd[,2]
        	}

cghwebinput = outputdir
setwd(cghwebinput)



for (i in 1:length(patients)) {
test_input=read.table(paste(patients[i],".txt",sep=""),sep="\t",header=T,as.is=T)
tested <- test_input[!is.na(test_input$Chromosome),]

x= runCGHAnalysis(tested, BioHMM = FALSE, UseCloneDists = FALSE, Lowess = FALSE,




					  Lwidth = 15, Wavelet = FALSE, Wlevels = 3, Runavg = FALSE,
					  Rwidth = 5, CBS = TRUE, alpha = 0.05, Picard = FALSE, Km = 20, 
					  S = -0.5, FusedLasso = FALSE, fluv = FALSE, FDR = 0.5,
					  rsm = FALSE, GLAD = FALSE, qlambda = 0.999, 
			  FASeg = FALSE, sig = 0.025, delta = 0.1, 
			  srange = 50, fineTune = FALSE, Quantreg = FALSE, lambda = 1,
			  minLR = -2, maxLR = 2, Threshold = 0.2, 
			  genomeType = "hg17", tempDir="../output",resultDir=paste(patients[i],"output",sep="_"))
			  }
			  
			  
			  
#setwd(paste(outputdir,Cancer,"/profiles/",sep=""))
setwd("../profiles/")
get_input_dir = "/var/www/cgi-bin/onlinetool/version_2/ominer_results/"

	library("R.utils")
	
		
		for (i in 1:length(patients)) {
		
		resultDir=paste(patients[i],"output", sep="_")
		cghweb_read=paste(get_input_dir,Cancer,"/cghweb/output/",resultDir,sep="");
		gunzip(paste(cghweb_read,"/Table_of_aCGH_smoothed_profiles.txt.gz",sep=""))
		file.copy(paste(cghweb_read,"/Table_of_aCGH_smoothed_profiles.txt",sep=""), paste(getwd(),"/",patients[i],".txt",sep=""))
	}
		
file.names=list.files(getwd())
	print(file.names)
#open first file and take the probeID, chr, position and summary
	first.data = read.table(file.names[1], sep="\t", header=T,as.is=TRUE)[,c("ProbeID","Chromosome","Position","Summary" )]
	colnames(first.data)=c("ProbeID","Chromosome","Position",sub(".txt","",file.names[1]))
	row.names(first.data)=first.data$ProbeID
	if(length(file.names)>1) {
		for (i in 2:length(file.names))
		{
			data=read.table(file.names[i], sep="\t", header=T,as.is=TRUE)[,c("ProbeID","Summary" )]	
			rownames(data) = data$ProbeID;	
			newcol=sub(".txt","",file.names[i]);	
			first.data[,newcol] = data[rownames(first.data),"Summary"]
		}
	}
	
	
	setwd("../Matrix")
	
	write.table(first.data,file="R_input.txt",sep="\t",quote=FALSE,row.names=FALSE)	
	file.copy("R_input.txt",paste("../../output/",sep=""))
	setwd("../output")
	system("rm -rf *_output")
		
	
}
			  