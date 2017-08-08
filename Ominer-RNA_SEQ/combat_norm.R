##########################################################################################
#
#File name:combat_norm.R
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
#####targets - this is a .txt file containing sample information, in an analysis using combat will contain four columns, filename, samplename and target (in the target column biological information regarding phenotype is found) and also a study column where samples 
#####from the same study have the same identifier
#####norm_matrix - this is the normalised expression matrix - "normalised.txt"
#####project - this is the name given to the analysis
#####This code  applies the combat algorithm to more than one dataset to combat batch effects 
run_combat <- function(targets,norm_matrix,project) {
	library(sva)
	targetFile <- read.table(targets,,header=T, row.name="Name", sep="\t")
	data <- read.table(paste("ominer_results/",project,"/","norm","/",norm_matrix,sep=""),sep="\t",as.is=TRUE,header=TRUE,row.names=1,check.names=FALSE)

datasets <- targetFile$Study

targets <- targetFile$Target

	rsd<-apply(data,1,sd)
data <- data[rsd>0,]

#####  Obtain known batch effects
batch = as.vector(datasets)

#####  Create model matrix for outcome of interest and other covariates besides batch
f.model <- factor(targets, levels = levels(factor(targets)))

mod <- model.matrix(~f.model)

print("Adjusting data for batch effects using ComBat...")
data_combat = ComBat(dat=data, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE, numCovs=NULL)
write.table(data_combat,file=paste("ominer_results/",project,"/","norm","/data_combat.txt",sep=""),sep="\t",row.names = TRUE, quote = FALSE);
	
	
	
	
	
}