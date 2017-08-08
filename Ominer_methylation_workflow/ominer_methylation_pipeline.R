############################################################################################################
#
#File name ominer_methylation_pipeline.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Master R script for the analysis of Raw idat files and unnormalised/normalised methylation data
#pre_requisites:
#chAMP
#limma
#Arguments needed 
#platform = 450k or 250k
#dataType = raw or normalised/unnormalised
#folder = pathe to the folder where files are uploaded from the user interface
#project = name given to analysis
#target = full path to target file
#csv = full path to the target file .csv format
#comp = full path to teh comparisons .comp file
#normalisation = normalisation method i.e. SWAN, BMIQ 
#filter = filter methos e.g SD, IQR
#filterval = i.e. 10,40 etc
#pvalue = pvalue to sort significant results
#foldchange - fc to sort sig results
#diffmethod = limma
#limmamethod = separate etc
#analysis = paired/unpaired
#replicates = yes/no
#combat = 0/1 - 1 it is applied 0 it is not
#champ - 1/0  - 1 it is applied and 0 it is not 

# Parse the command line options
for (e in commandArgs()) {
	ta = strsplit(e,"=",fixed=TRUE)
	ta
#	if(! is.na(ta[[1]][2])) {
		temp = ta[[1]][2]
		assign(ta[[1]][1],temp)
		cat("assigned ",ta[[1]][1]," the value of |",temp,"|\n")
#	} 
}

setwd("/var/www/cgi-bin/onlinetool/version_2/")
dir = paste("ominer_results",project,sep="/")
#setup directories to write output of results 
dir.create (dir,mode="775")   ###make a directory called project under OMINER-OUTPUT each time that this script is run
qcdir = paste(dir,"QC",sep="/")
dir.create(qcdir)
target_path = "../../../html/onlinetool/temp"
###for now during testing hash above out
#target_path = "ominer_results"
target_file = paste(target_path,project,"target.csv",sep="/")
limma_target_file =  paste(target_path,project,"target.txt",sep="/")
comp_path = "../../../html/onlinetool/temp"
#comp_path = "ominer_results"
comp_file =  paste(comp_path,project,"comp.txt",sep= "/")
normdir = paste(dir,"norm",sep="/")
dir.create(normdir)
file.copy(target_file,qcdir)
location <- getwd()
   

csv_file_path <- paste(target_path,project,sep="/")

#perform QC
library("ChAMP")
if (dataType == "raw") {
####Need to change line below to accomodate i.dat files uploaded from the web
print ("This is my target_file path")
print (csv_file_path)
champ_data <- champ.load(directory = csv_file_path, methValue = "B",resultsDir = qcdir,filterXY=TRUE,QCimages=TRUE,filterDetP = TRUE, detPcut = 0.01, removeDetP = 0, filterBeads=TRUE, beadCutoff=0.05, filterNoCG=FALSE)

#####write out unnormalised data

write.table(champ_data$beta,paste("ominer_results/",project,"/norm/normalised.txt",sep=""),row.names = TRUE, quote = FALSE,sep="\t")


######Normalize data 

upper_normalisation <- toupper(normalisation)


champ_normalized <- champ.norm(beta = champ_data$beta, rgSet = champ_data$rgSet, pd = champ_data$pd, mset = champ_data$mset,
sampleSheet = target_file, resultsDir = qcdir, methValue = "B", fromIDAT = TRUE, norm = upper_normalisation, fromFile = FALSE, betaFile,
filter = TRUE, filterXY = TRUE, QCimages = TRUE, plotBMIQ = TRUE)
}

upper_normalisation <- toupper(normalisation)
if (dataType == "normalised") {
print ("yes")
	####copy normalised data matrix to correct place
	norm_file_path = paste("../../../html/onlinetool/temp",project,sep="/")
        
normalisation =  paste(norm_file_path,"normalised.txt",sep= "/")
	        file.copy(normalisation,normdir)
		#####upload the relevant file via the web and use limma (convert bete vlues to m values????)
	

}

######write out all beta values to a nnormalized text file - ready for input to limma



	

	
	
	


DEdir = paste(dir,"DifferentialExpression",sep="/")
dir.create(DEdir)
source("run_limma_methy.R")

###Need to remove this once the input options have been changed from IDAT files and unnormalised to IDAT files and normalised
#####When it has been changed remove both if statements below and keep data as = "normalised.txt"
if (dataType == "raw") {
data = "normalised.txt"
}

if (dataType == "normalised") {
data = "normalised.txt"
}

run_limma(limma_target_file,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange,platform,project)
####run limma



######Clustering based on the normalised methylation values 

######Take this line out once the normalised file is saved as normalised.txt

#if (dataType == "unnormalised") {
    
#}

clusterdir = paste(dir,"cluster",sep="/")
dir.create(clusterdir)
source("cluster_code_methy.R")
 pd <- read.table(target_file, header = T, 
            sep = "\t")
             lev <- unique(pd$Target)
k <- length(lev)
print ("This is the number of clusters")
print(k)
cluster_all(limma_target_file,k,project)
testObject <- function(object)
	{
		exists(as.character(substitute(object)))
	}
        
if(testObject(venn)) {
    

source("generate_venn.R")

generate_venn("ebA_allresults.txt",project)
}



	if(testObject(gostats)) {
go_analysis = "0"
comps <- read.table(comp_file, sep = "\t", as.is = T, header = F, 
        strip.white = T)
        for (i in 1:length(comps)) {
                print(comps[i])
        fg <- strsplit(comps[,i],"=")
        ci <- fg[[i]][1]
####need to do this for all comparisons
        #filteredx <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"Filtered.txt",sep=""),row.names = 1, header = T, sep = "\t") 
#sigGene=rownames(filteredx)
filtered_x <- paste("ominer_results/",project,"/","DifferentialExpression/",ci,"_Filtered.txt",sep="")
if (file.exists(filtered_x)) {
	go_analysis = "1"
	godir = paste(dir,"GO",sep="/")
dir.create(godir)
}
if (go_analysis == "1") {
	source("methylation_GO_analysis.R")
gostats_analysis(platform,project,comp_file)
}
}
}

setwd(qcdir)
move_file_command <- paste("mv","target.csv","target.txt",sep=" ")
system(move_file_command)
setwd("../")
resultdir = "/var/www/cgi-bin/onlinetool/version_2/ominer_results/"
setwd(resultdir)
zipcommand=paste("zip -r ",project,".zip ",project,sep="")
print(zipcommand)
system(zipcommand)
file.rename(paste(resultdir,project,".zip",sep=""),paste(resultdir,project,"/",project,".zip",sep=""))

command="/var/www/cgi-bin/onlinetool/version_2/mail.pl"
url="Your results are now available at : http://o-miner.org/cgi-bin/onlinetool/version_2/methOutputView.pl?"
subject <- "O-miner Results"
mailcommand=paste(command," ",email," \"",subject,"\" \"",url,project,"\"",sep="")

#mailcommand=paste(command," ",email," ",project," \"",url,project,"\"",sep="")
system(mailcommand)
