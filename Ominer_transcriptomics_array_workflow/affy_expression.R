##########################################################################################
#
#File name: affy_expression.R
#
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#
#
#
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
##########################################################################################
##########################################################################################
#Description:Master R script for running analysis of array based data from Affymetrix platforms. This R script calls several different functions. In order to run this pipeline you need to have prepared two separate files:
#a target file and a comparisons file - examples of these files can be downloaded with this script. The script can be run from two different data entry points these are: (1) Raw CEL files and (2) Normalised expression
#matrix.A target file for straightforward analysis will have three columns (1) Name, (2) FileName and (3) Target. Additional columns are required in the target file if analysis is run as paired (target_paired.txt), if replicates (target_replicates.txt) are
#present and if combat is to be used in a meta-analysis (target_combat). Please note that the target file is always refered to as 'target.txt' regardless of whether analysis is run as paired or whether replicates are present or if combat is used. 
#Results are written out at each stage of the analysis to a subdirectory within the project folder that is created within the 'ominer_results' folder.
#pre-requisites: Users must have the following R packages installed on their system: (1) affy, (2) simpleaffy, (3) arrayMvout, (4) the annotation.db package for the array platform from which you are analysing data - for example hgu133plus2.db for the hgu133plus2
#array,(5) affyPLM, (6) arrayQualityMetrics, (7) annotate, (8) limma, (9) sva, (10) estimate, (11) gplots, (12) FC14.plotting.lib and  (13) survival. 
#Before running this pipeline you will need to have created a folder named 'ominer_results'.
#Code dependencies:
#Venn.R
#Code.R
#bcctb.utils.R
#The code dependencies above are distributed with this master script.
#The command line arguments needed to run this script are: 
#R--file-/scratch/ominer/transcriptomics/Affymetrix/affy_expression.R --args platform=hgu133plus2 dataType=CEL folder=/scratch/ominer/transcriptomics/Affymetrix/CEL project=a.sangaralingam-AFFY_TEST 
#target=/scratch/ominer/transcriptomics/Affymetrix/target.txt comp=/scratch/ominer/transcriptomics/Affymetrix/comp.txt normalization=rma filter=sd filterval=40 adjust=BH pvalue=0.05 foldchange=2 limmamethod=separate analysis=unpaired replicates=no combat=0
#survival=0 estimate=0 aqm=0 gostats=1
#Example of shell script containing command line arguments: mouse_exp_v2.sh
#Description of command line arguments:
#platform: Affymetrix platform from which data is being analysed
#dataType: Whether data analysed is raw data or a normalised expression matrix i.e. CEL or normalised
#folder: full path to where the input data is found
#project: name given to the analysis project
#target: full path to the target file
#normalization: name of normalisation method to use e.g. rma, gcrma 
#filter: method used to filter data e.g. sd, iqr
#filterval: % value to filter the normalised expression matrix i.e to reatin the topmost 40% filterval = 40
#adjust: method to adjust the false discovery rate (FDR) within limma e.g BH (Benjamini-Hochberg) or BY
#pvalue: adjusted p-value cutoff to report differentially expressed genes
#foldchange: log2 foldchange cutoff to report differentially expressed genes
#limmamethod:Method used within limma when comparing groups
#analysis: whether an unpaired or a paired analysis is executed
#replicates: yes/no indicating whether technical replicates are present within the dataset
#combat:0/1 if the combat algorithm is to be applied to the data
#aqm: 0/1 if arrayQualityMetrics 
#survival: 0/1 if survival analysis is to be executed
#estimate: 0/1 if the ESTIMATE algorithm is to be used
#gostats: 0/1 if gostats analysis is to be executed
#Description of output: Data is output to the ominer_results directory, a folder is automatically created within this directory with your project name. Subfolders are created within this directory such as (1) QC, (2) norm, (3) Differential Expression
############################################################################################################################################################
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
#####use this to check whether the directory already exists in the results folder
#####if dir.exists(paths)
dir.create (dir,mode="775")   ###make a directory called project under OMINER-OUTPUT each time that this script is run
qcdir = paste(dir,"QC",sep="/")
dir.create(qcdir)
target_path = "../../../html/onlinetool/temp"
target_file = paste(target_path,project,"target.txt",sep="/")
comp_path = "../../../html/onlinetool/temp"
comp_file =  paste(comp_path,project,"comp.txt",sep= "/")
normdir = paste(dir,"norm",sep="/")
dir.create(normdir)

errorFile = paste(dir,"/run_log.txt",sep="")
#fileerr = paaste(dir,"/error.txt",sep="")
file.create(errorFile)
#file.create(fileerr)

file.copy(target_file,qcdir)
target_path = "../../../html/onlinetool/temp"
folder= paste(target_path,project,sep="/")
DEdir = paste(dir,"DifferentialExpression",sep="/")
dir.create(DEdir)
mydir=getwd()

library("affy")
	library("simpleaffy")
	library("arrayMvout")
	library("affyPLM")
	
if (platform == "hugene1.1" || platform == "hugene2.0") {
#if (aqm == "1") {
library("oligo")
pd <- read.AnnotatedDataFrame(target_file,header=T, row.name="Name", sep="\t")
setwd(folder)
    lev <- unique(pd$Target)
x <- varMetadata(pd)
x <- data.frame(x, channel = "_ALL_")
varMetadata(pd) <- x
dat <- read.celfiles(filenames = pd$FileName, sampleNames = sampleNames(pd), phenoData=pd)
setwd(mydir)
#source("QC_affy_mod.R")
#runQC(dataType,dat,target_file,project)

##### Read in the probe annotation file
#probesAnnotFile = read.table("Affy_HuGene1ST_annot_collapsed.txt",sep="\t",as.is=T,header = T)
if (aqm== 1) {
aqmdir <- paste(qcdir,"AQM",sep="/")
library("arrayQualityMetrics")
preparedData = prepdata(expressionset = dat, intgroup = "Target", do.logtransform = TRUE)

QCboxplot <- aqm.boxplot(preparedData, subsample=20000, outlierMethod = "KS")
QCdensity <- aqm.density(preparedData)
QCheatmap <- aqm.heatmap(preparedData)
QCmaplot <- aqm.maplot(preparedData, subsample=20000, Dthresh=0.15, maxNumArrays=8, nrColumns=4)
QCmeansd <- aqm.meansd(preparedData)

qm = list("Boxplot"=QCboxplot, "Density"=QCdensity, "Heatmap"=QCheatmap, "MAplot"=QCmaplot, "MeanSD"=QCmeansd)

aqm.writereport(modules = qm, reporttitle = "arrayQualityMetrics report for HuGene_1.1_ST_v1 data", outdir = aqmdir, arrayTable = pData(dat))
#Data is written out to file index.html in the AQM directory
}
print ("YES, I got to Here")
source("hugene1.1_normalisation.R")
Normalisation(normalisation,target_file,dat,project)
source("hugene1.1_filtering.R")
filtering(filter,filterval,combat,project)
source("run_limma_affy_hugene.R")
data <- "filtered_data.txt"

run_limma(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange,platform,project)
clusterdir = paste(dir,"cluster",sep="/")
dir.create(clusterdir)
source("cluster_code_affy.R")

#####number of clusters k <- lev (this is the number of biological groups)
k <- length(lev)
print ("This is the number of clusters")
print(k)
cluster_all(target_file,k,project)

testObject <- function(object)
	{
		exists(as.character(substitute(object)))
	}

if(testObject(gostats)) {
    print ("user asked for gostats")

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
filtered_x <- paste("ominer_results/",project,"/DifferentialExpression/",ci,"_Filtered.txt",sep="")
if (file.exists(filtered_x)) {
	go_analysis = "1"
	print ("file filtered_x exists go_analysis =1")
	godir = paste(dir,"GO",sep="/")
dir.create(godir)
}
if (go_analysis == "1") {
	source("hugene1.0_GO_analysis.R")
hugene1.0_GO(project,comp)
}
}
}
resultdir = "/var/www/cgi-bin/onlinetool/version_2/ominer_results/"
setwd(resultdir)
zipcommand=paste("zip -r ",project,".zip ",project,sep="")
print(zipcommand)
system(zipcommand)
file.rename(paste(resultdir,project,".zip",sep=""),paste(resultdir,project,"/",project,".zip",sep=""))



command="/var/www/cgi-bin/onlinetool/version_2/mail.pl"
url="Your results are now available at : http://o-miner.org/cgi-bin/onlinetool/version_2/exprOutputView_dec_2013.pl?"

mailcommand=paste(command," ",email," ",project," \"",url,project,"\"",sep="")
#mailcommand=paste(command," ",email," ",url,project,"\"",sep="")
system(mailcommand)


}

 #outaqm = paste("ominer_results/",project,"/QC","/AQM",sep="")
                #aqm=try(arrayQualityMetrics(expressionset = dat, outdir = outaqm, force = TRUE, do.logtransform = FALSE))
        
# QC summary data frame
#qc.Adata <- qc(dat)
        
       # pdtable <- read.table(target_file,sep="\t",as.is=T,header=T,row.names="Name")
     #qcsummary<-data.frame(pdtable,avbg(qc.Adata),percent.present(qc.Adata),ratios(qc.Adata))
       #write.table(qcsummary, paste("ominer_results/",project,"/QC","/qcsummary.txt",sep=""), sep = "\t",quote=F,row.names=FALSE)


library(paste(platform,".db",sep=""),character.only = TRUE)
if (dataType == "CEL") {
 pd <- read.AnnotatedDataFrame(target_file, header = T, row.name = "Name", sep = "\t")
             lev <- unique(pd$Target)
            ####get the number of clusters here for the cluster plot K = number of biological groups....
            setwd(folder)
	    #cat("Reading CEL files",file=errorFile,sep="\n",append=TRUE)

        dat <- ReadAffy(filenames = paste(pd$FileName, sep = ""), 

            sampleNames = sampleNames(pd), phenoData = (pd))
setwd(mydir)

#if (combat == "1") {
####### Run QC for each of the datasets separately
#source("combat_qc_extra_v1.R")
#runQC(dataType,platform,analysis,project)


#new_pdtable <- read.table(target_file,header=T, sep="\t")
#lev <- unique(new_pdtable$Target)
#no_of_files <- length(lev)
#revised = "0"

#for (i in 1:length(lev)) {

#out_file <- paste("outliers_",i,".merge.txt",sep="")
#out <- read.table(paste("ominer_results/",project,"/QC/",out_file,sep=""),sep = "\t", as.is = T, header = T, strip.white = TRUE)
#if (nrow(out) > 1) {

#revised = "1"
#samples_to_remove <- out$samp
#list_of_samples_to_remove <- setdiff(new_pdtable$FileName,samples_to_remove)
#a <- new_pdtable$FileName %in% list_of_samples_to_remove
#y <- subset(new_pdtable,a)
#t <- y[1:4]
#new_pd <- new_pdtable[new_pdtable$FileName %in% list_of_samples_to_remove,]

#}
#write.table(t,paste("ominer_results/",project,"/QC","/target_qc.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
#}

#if (revised == 1) {
#target_qc_file = paste("ominer_results/",project,"/QC","/target.txt",sep="")
#pd <- read.AnnotatedDataFrame(target_qc_file, header = T, row.name = "Name", 
            #sep = "\t")
            #setwd(folder)
        #dat <- ReadAffy(filenames = paste(pd$FileName, sep = ""), 
            #sampleNames = sampleNames(pd), phenoData = (pd))
            #setwd(mydir)
#}
#setwd(qcdir)
#system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=qc.pdf qc_*.pdf")
#system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=PLM.pdf PLM_*.pdf")
#system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=deg.pdf deg_*.pdf")

#system("cat NUSE_*.txt > NUSE.txt")
#system("cat RLE_*.txt > RLE.txt")
#setwd("../../../")

#}
	




	
    
    
    
    
    
    
    


###if (combat == "0") {  Hashed out 22.01.15 not sure why QC would only be run if combat was not used 




source("QC_affy_mod.R")
#echeck <- paste("/",errorFile,sep="")

cat ("Running QC steps",file=errorFile,sep="\n",append=TRUE)

runQC(dataType,dat,target_file,project)   ####all output goes to the QC directory
#}




if (aqm =="1") {

	source("aqm.R")

	#run_aqm(target_file,dat,analysis,project)

	

}
revised = "0"
####If some of the samples did not pass QC then remove these from the analysis
if (combat == "0") {
out <- read.table(paste("ominer_results/",project,"/QC","/outliers.txt",sep=""), sep = "\t", as.is = T, 
        header = T)
    colnames(out) = c("Sample", "averageBG", "ScaleFactor", "Present", 
        "HSAC07", "GAPDH", "NUSE", "RLE", "RLE_IQR", "RNAslope")
    if (nrow(out) > 0) {
    	revised = "1"
    	}
    	 if (nrow(out) == 0){
    		revised = "0"
    	    		file.copy(target_file,qcdir)
    	}
    	
    	if (revised == 1) {

    	target_qc_file = paste("ominer_results/",project,"/QC","/target_qc.txt",sep="")
        	 pd <- read.AnnotatedDataFrame(target_qc_file, header = T, row.name = "Name", 
            sep = "\t")
            setwd(folder)
        dat <- ReadAffy(filenames = paste(pd$FileName, sep = ""), 
            sampleNames = sampleNames(pd), phenoData = (pd))
            setwd(mydir)
    	}    	
    
    	print ("This is revised")
	print (revised)
	}
	
    	source("affy_normalisation_mod.R")
cat ("Normalising data",file=errorFile,sep="\n",append=TRUE)

if (revised == 0) {
Normalisation(normalisation,target_file,dat,project)
}
if (revised == 1) {
Normalisation(normalisation,target_qc_file,dat,project)
}
}

if (dataType == "normalised") {
	####copy normalised data matrix to correct place
	norm_path = "../../../html/onlinetool/temp"
normalisation =  paste(norm_path,project,"normalised.txt",sep= "/")
	        file.copy(normalisation,normdir)
		dataset = paste("../../../html/onlinetool/temp",project,"normalised.txt",sep="/")
######if Users upload entire normalised matrix and then delete some entries then the columns they deleted need to be reoved from the normalised matrix
t_file <- read.table(target_file, header = T, sep = "\t",as.is=TRUE,row.names=NULL)

n_matrix <- read.table(dataset,header=T,sep="\t",as.is=TRUE,row.names=1)
lev <- unique(t_file$Target)


	changed_file = "0"

total_length <- length(colnames(n_matrix))
target_length <- nrow(t_file)

if (total_length > target_length) {
	changed_file = "1"
}

if (changed_file == "1") {
detection_names = NULL
n <- length(t_file$Name)
for (i in seq(from=1,to=n,by=1)) {
	
	#wanted <- paste("\\",t_file$Name[i],"\\b",sep="")
	wanted <- paste("^",t_file$Name[i],"$",sep="")
	print ("This is the name of the file")
print (wanted)
####EXACT MATCH
#ed <- paste("^",wanted,"$",sep="")
detection_names <- c(detection_names,grep(wanted,colnames(n_matrix)))
print ("THIS is detection names")
print (detection_names)
#detection_names <- c(detection_names,grep(ed,colnames(n_matrix)))
}


 new_dat <- subset(n_matrix[detection_names])
 #new_dat <- subset(n_matrix,select=wanted)
print ("THIS IS THE FILTERED DATA FILE")
print (new_dat)
#TargetID <- rownames(n_matrix)
#new_matrix <- cbind(TargetID,new_dat)

write.table(new_dat,paste("ominer_results/",project,"/norm","/normalised.txt",sep=""), sep="\t", quote=FALSE, row.names=T)

}
	}

#norm_matrix = paste(normalisation,".exp",sep="")
norm_matrix = "normalised.txt"
if (estimate == "1") {
	estimatedir = paste(dir,"estimate",sep="/")
dir.create(estimatedir)
	source("estimate_ominer_function.R")
	cat ("Running ESTIMATE for tumor purity",file=errorFile,sep="\n",append=TRUE)
	run_estimate(project,platform,norm_matrix,target_file)
	}
	
if (survival == "1") {
	survivaldir = paste(dir,"survival",sep="/")
dir.create(survivaldir)
cat ("Running survival analysis",file=errorFile,sep="\n",append=TRUE)
source("KM_surv_plots.R")
survival_os_plot(target_file,project)
	
}	
  if (combat == "1") {

	source("combat_norm.R")
	if (revised =="0") {
	cat ("Correcting for batch effect using COMBAT",file=errorFile,sep="\n",append=TRUE)
		run_combat(target_file,norm_matrix,project)
	}

	if (revised == "1") {
		target_file = target_qc_file
			cat ("Correcting for batch effect using COMBAT",file=errorFile,sep="\n",append=TRUE)
		run_combat(target_file,norm_matrix,project)
	}

}


print ("I AM NOW DOING THIS PART")

source("filtering.R")
cat ("filtering normalised expression matrix",file=errorFile,sep="\n",append=TRUE)

filtering(filter,filterval,combat,project)

if (dataType == "filtered") {
	####copy filtered matrix to correct place
	}



source("run_limma_affy.R")

cat ("Running limma for differential expression analysis",file=errorFile,sep="\n",append=TRUE)

data = "filtered_data.txt"

if (dataType == "CEL") {
if (revised == "0") {
print ("I am doing unrevised")
run_limma(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange,platform,project)

}
if (revised == "1") {
print ("I am doing revised")
	target_file = target_qc_file
	run_limma(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange,platform,project)
}
}

if (dataType == "normalised") {
	run_limma(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange,platform,project)
}
clusterdir = paste(dir,"cluster",sep="/")

dir.create(clusterdir)

source("cluster_code_affy.R")
cat ("Clustering data",file=errorFile,sep="\n",append=TRUE)
#####number of clusters k <- lev (this is the number of biological groups)
k <- length(lev)
print ("This is the number of clusters")
print(k)
cluster_all(target_file,k,project)

comps <- read.table(comp_file, sep = "\t", as.is = T, header = F, 
        strip.white = T)
        for (i in 1:length(comps)) {
                print(comps[i])
        fg <- strsplit(comps[,i],"=")
        ci <- fg[[i]][1]
####need to do this for all comparisons
        #filteredx <- paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"Filtered.txt",sep="")
#sigGene=rownames(filteredx)
if (platform == "st1") {
filtered_x <- paste("ominer_results/",project,"/","transcript/DifferentialExpression/",ci,"_Filtered.txt",sep="")
}
else {
filtered_x <- paste("ominer_results/",project,"/","DifferentialExpression/",ci,"_Filtered.txt",sep="")
}

if (file.exists(filtered_x)) {
filtered_x <- read.table(paste("ominer_results/",project,"/","DifferentialExpression/",ci,"_Filtered.txt",sep=""),header = T, sep = "\t",as.is=T)
length_AGI <- length(filtered_x$AGI)
if (length_AGI > 2) {
source("heatmap_create_TRIMMED.R")
cat ("Generating heatmap",file=errorFile,sep="\n",append=TRUE)
if (dataType == "CEL") {
if (revised == "0") {
heatmap_generate(target_file,project,comp_file)
}
if (revised == "1") {
target_file = target_qc_file
heatmap_generate(target_file,project,comp_file)
}
}
if (dataType == "normalised") {
heatmap_generate(target_file,project,comp_file)
}
}
}
}

testObject <- function(object)
	{
		exists(as.character(substitute(object)))
	}
        
if(testObject(venn)) {
    


source("generate_venn.R")



generate_venn("ebA_allresults.txt",project)
cat ("Generating venn diagram",file=errorFile,sep="\n",append=TRUE)

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
if (platform == "st1") {
filtered_x <- paste("ominer_results/",project,"/","transcript/DifferentialExpression/",ci,"_Filtered.txt",sep="")
}
else {
filtered_x <- paste("ominer_results/",project,"/","DifferentialExpression/",ci,"_Filtered.txt",sep="")

}

if (file.exists(filtered_x)) {
	go_analysis = "1"
	godir = paste(dir,"GO",sep="/")
dir.create(godir)
}
if (go_analysis == "1") {
	source("GOStats_function.R")
	cat ("Running Gene Ontology analysis",file=errorFile,sep="\n",append=TRUE)
gostats_analysis(platform,norm_matrix,project,comp_file)
}
}
}

resultdir = "/var/www/cgi-bin/onlinetool/version_2/ominer_results/"
setwd(resultdir)
zipcommand=paste("zip -r ",project,".zip ",project,sep="")
print(zipcommand)
system(zipcommand)
file.rename(paste(resultdir,project,".zip",sep=""),paste(resultdir,project,"/",project,".zip",sep=""))



#uncomment out the four lines below.....
command="/var/www/cgi-bin/onlinetool/version_2/mail.pl"
url="Your results are now available at : http://o-miner.org/cgi-bin/onlinetool/version_2/exprOutputView_dec_2013.pl?"
subject <- "O-miner Results"
mailcommand=paste(command," ",email," \"",subject,"\" \"",url,project,"\"",sep="")

system(mailcommand)


#source("heatmap_create.R")



#DE_list = "TNPvstumFiltered.txt"  ###if user picks comparison DE_list = paste(comp[i],"_all_DE.txt") - comp[i] realted to file comp.txt containing comparisons

#heatmap_generate(target,DE_list,project)



	