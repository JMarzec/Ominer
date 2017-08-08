##########################################################################################
#
#File name: estimate_ominer_function.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#
#
#
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
##########################################################################################
##########################################################################################
#Description: function to run the ESTIMATE algorithm to estimate tumor purity
#project - name given to the project
#platform - name of platform from which the array data to be analysed is from
#norm - file - "normalised.txt"
#
################################################################################################KM plots - overall survival from a target file generates .txt file and KM plots for 5,10 and 15 years
######
#######Function to run the estimate algorithm on a normalised gene expression matrix 
####Inputs to teh function are a :
####normalised gene expression matrix
####project == study name
#####platform = probe platform that has been used in the experiment
#####normalisation method that was used to normalise that data

#####Need to intergrate the traget file whre only the cancer samples are present in the input dataset as input to this function.
####This is to make sure that only the expression of teh cancer samples is taken from the normalised expression matices
#####estimate will be run as an option if users chosse affymetrix expression analysis - need to make sure that the reveant .db will be loaded
####epennding on the affy platform taht is being used.

######Ideas of hpw to implement this an extra tick box appears when the target file is visulaised on screen so user can mark which of the samples
#####are cancer samples OR just have another box appear if they tick estimate then they input in this box what they have called there cancer samples?
######When this has been decided unmark the hashed out code

#####Alos look at the estimate output again to see which ones are the best to visualise in the estimate paper they have boxplots of the 
####immune score and stroml score an tumor purity - prduce good quality graphical output for ominer 

run_estimate <- function (project,platform,norm,target_file) {
	#target.data <- as.matrix( read.table(paste("ominer_results/",project,"/QC/",target_file,sep=""),row.names = 1, header = T, sep = "\t") )
	#pd <- read.table(target_file, sep="\t", as.is=TRUE,strip.white=TRUE,header=T);

	#target.data <- target.data[target.data[,"Target"]=="cancer",]

	#ep.data = paste(norm,".exp",sep="")
exp.data <- as.matrix( read.table(paste("ominer_results/",project,"/norm/",norm,sep=""),row.names = 1, header = T, sep = "\t") )
  #  intersect.patients <- intersect( colnames(exp.data), rownames(target.data));
      # exp.data <- exp.data[ , intersect.patients ];
#target.data <- target.data[ intersect.patients , ];
source("bcctb.utils.R")
 
      #library("hgu133plus2.db")
        exp.data<-getEntrez(exp.data,platform,"SYMBOL")
      
write.table(exp.data, paste("ominer_results/",project,"/","estimate","/","exp_data.txt",sep=""), sep = "\t", quote = FALSE, 
        row.names = TRUE, col.names = NA)


library("estimate")
estfile = "est_file.out"
#####need the full path to exp_data.txt
filterCommonGenes(paste("ominer_results/",project,"/","estimate","/","exp_data.txt",sep=""),paste("ominer_results/",project,"/","estimate/",estfile,sep=""))
########if none are found in imm or str score then estimate score cannot be generated
scrfile ="scrfile.out.txt"
estimateScore(paste("ominer_results/",project,"/","estimate/",estfile,sep=""),paste("ominer_results/",project,"/","estimate/",scrfile,sep=""))

estimate.output =read.table(file=paste("ominer_results/",project,"/","estimate/",scrfile,sep=""),sep="\t",skip=2,header=T,row.names=1)
stroma.score=as.numeric(as.vector(estimate.output[1,2:length(estimate.output)]))
immune.score=as.numeric(as.vector(estimate.output[2,2:length(estimate.output)]))
estimate.score=as.numeric(as.vector(estimate.output[3,2:length(estimate.output)]))
tumour.purity=as.numeric(as.vector(estimate.output[4,2:length(estimate.output)]))

png(paste("ominer_results/",project,"/","estimate/","graphical_output.png",sep=""))
par(mfrow=c(3,1)) 
hist(stroma.score)
hist(immune.score)
hist(estimate.score)
dev.off()

png(paste("ominer_results/",project,"/","estimate/","boxplot_output.png",sep=""))
par(mfrow=c(3,1)) 
boxplot(stroma.score)
boxplot(immune.score)
boxplot(estimate.score)
dev.off()

####transpose estimate output for output purposes
transposed_output <- t(estimate.output)
write.table(format(transposed_output,digits=2),paste("ominer_results/",project,"/","estimate","/","tp_output.txt",sep=""),row.names = T,quote=FALSE, col.names = F, sep = "\t" );

#write out tumour purity to target file

#target.data <- cbind(target.data, "tumour.purity" = tumour.purity);
#write.table( target.data, file = target.file.out, row.names = T, col.names = NA, sep = "\t" );



}

