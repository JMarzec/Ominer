##########################################################################################
#
#File name: illumina_exp.R
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
#Description:Master R script for running analysis of array based data from Illumina platforms. This R script calls several different functions. In order to run this pipeline you need to have prepared two separate files:
#a target file and a comparisons file - examples of these files can be downloaded with this script. The script can be run from two different data entry points these are: (1) Raw CEL files and (2) Normalised expression
#matrix.A target file for straightforward analysis will have three columns (1) Name, (2) FileName and (3) Target. Additional columns are required in the target file if analysis is run as paired (target_paired.txt), if replicates (target_replicates.txt) are
#present and if combat is to be used in a meta-analysis (target_combat). Please note that the target file is always refered to as 'target.txt' regardless of whether analysis is run as paired or whether replicates are present or if combat is used. 
#Results are written out at each stage of the analysis to a subdirectory within the project folder that is created within the 'ominer_results' folder.
#pre-requisites: Users must have the following R packages installed on their system: (1) arrayQualityMetrics , (3) arrayMvout, (4) biomaRt
#array,(5) affyPLM, (6) limma, (7) sva,  (8) gplots, (9) FC14.plotting.lib and  (10) survival, (11) gostats, (12) lumi 
#Before running this pipeline you will need to have created a folder named 'ominer_results'.
#Code dependencies:
#Venn.R
#Code.R
#bcctb.utils.R
#The code dependencies above are distributed with this master script.
#The command line arguments needed to run this script are: 
#R --file=/data/BCI-BioInformatics/ominer/illumina_exp.R  --args   platform=ht12v3 dataType=normalised folder=/var/www/html/onlinetool/temp/a.sangaralingam-illumina_expression_array project=a.sangaralingam-illumina_expression_array
#target=/var/www/html/onlinetool/temp/a.sangaralingam-illumina_expression_array/target.txt comp=/var/www/html/onlinetool/temp/a.sangaralingam-illumina_expression_array/comp.txt normalisation=rsn filter=sd filterval=40 adjust=BH pvalue=0.05
#foldchange=0.5 limmamethod=separate analysis=unpaired replicates=no combat=0 aqm=0 survival=0 gostats=0
#Example of shell script containing command line arguments: a.sangaralingam-illumina_expression_array.sh
#Description of command line arguments:
#platform: Illumina platform from which data is being analysed i.e. array platform illumina-humanht-12-v3 = ht12v3
#dataType: Whether data analysed is normalised or unnormalised
#folder: full path to where the input data is found
#project: name given to the analysis project
#target: full path to the target file
#comp: full path to comparisons file
#normalisation: name of normalisation method to use e.g. rma, gcrma 
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
mydir = getwd()

dir = paste("ominer_results",project,sep="/")
dir.create (dir)   ###mkae a directory called project under OMINER-OUTPUT each time taht this script is run
qcdir = paste(dir,"QC",sep="/")
dir.create(qcdir)
mydir=getwd()
normdir = paste(dir,"norm",sep="/")
dir.create(normdir)
dedir = paste(dir,"DifferentialExpression",sep="/")
dir.create(dedir)
clusterdir = paste(dir,"cluster",sep="/")
dir.create(clusterdir)
target_path = "../../../html/onlinetool/temp"
target_file = paste(target_path,project,"target.txt",sep="/")
comp_path = "../../../html/onlinetool/temp"
comp_file =  paste(comp_path,project,"comp.txt",sep= "/")
errorFile = paste(dir,"/run_log.txt",sep="")
file.create(errorFile)

#setwd(folder)


#norm_file <- read.table(dataset,sep="\t",as.is=T,header=F,strip.white=T)
dataset = paste("../../../html/onlinetool/temp",project,"normalised.txt",sep="/")
######if Users upload entire normalised matrix and then delete some entries then the columns they deleted need to be reoved from the normalised matrix
t_file <- read.table(target_file, header = T, sep = "\t",as.is=TRUE,row.names=NULL)
n_matrix <- read.table(dataset,header=T,sep="\t",as.is=TRUE)
lev <- unique(t_file$Target)


	changed_file = "0"

total_length <- length(colnames(n_matrix))
target_length <- nrow(t_file)

if (total_length > (3*(target_length)+1)) {
	changed_file = "1"
}

if (changed_file == "1") {
detection_names = NULL
n <- length(t_file$Name)
for (i in seq(from=1,to=n,by=1)) {
	
	#wanted <- paste(t_file$Name[i],sep="")
	wanted <- paste("\\",t_file$Name[i],"\\b",sep="")
	print ("THESE ARE WANTED")
print (wanted)
detection_names <- c(detection_names,grep(wanted,colnames(n_matrix)))
}


 new_dat <- subset(n_matrix[detection_names])

TargetID <- n_matrix$TargetID
new_matrix <- cbind(TargetID,new_dat)

write.table(new_matrix,paste("ominer_results/",project,"/norm","/refined.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)

}


		
#setwd(mydir)


if (changed_file == "0") {
		dataset = paste("../../../html/onlinetool/temp",project,"normalised.txt",sep="/")
		}
		if (changed_file == "1") {
			dataset = paste("ominer_results/",project,"/norm","/refined.txt",sep="")
		}

##########If samples are found to fail quality control then they are automatically taken out of the analysis. However, in this version of the script don't need to rerun the qualiy control analysis

if (dataType == "unnormalised") {
	#dataset = paste("../../../html/onlinetool/temp",project,"normalised.txt",sep="/")
	#dataset = paste("ominer_results/",project,"/norm","/normalised.txt",sep="")
	#dataset = second_matrix
	print (dataset)
source("illumina_exp_qc.code_mod.R")
illumina_qc(dataset,project,target_file,analysis,normalisation)
cat ("Running quality control steps",file=errorFile,sep="\n",append=TRUE)

#out <- read.table(paste("ominer_results/",project,"/QC","/outliers.txt",sep=""), sep = "\t", as.is = T, 
 #       header = T)
  #  colnames(out) = c("Sample", "mean", "standard deviation", "detection_rate", 
   #     "distance_to_sample_mean")
    #if (nrow(out) > 0) {
    #	revised = "1"
    #		#file.copy(target_file,qcdir)
    #	}
    #	 if (nrow(out) == 0){
    	#	revised = "0"
    	 #   		file.copy(target_file,qcdir)
    	 
    	#}
    	
    	
    	revised = "0"
    	
    	    	outlier_file <- paste("ominer_results/",project,"/QC","/outliers.txt",sep="")
		if (file.exists(outlier_file)) {
    	revised = "1"
    	
    	}
    		if (revised == "1") {
    			

    	target_qc_file = paste("ominer_results/",project,"/QC","/target_qc.txt",sep="")
    	normalised_mod_file = paste("ominer_results/",project,"/QC","/new_dataset",sep="")
}
  if (revised == "0") {
	
	file.copy(target_file,qcdir)
}          	   
}





if (aqm == "1") {
	
	source("illum_aqm.R")
	#dataset = paste("../../../html/onlinetool/temp",project,"normalised.txt",sep="/")
	#dataset = paste("ominer_results/",project,"/QC","/modified_dataset.txt",sep="")
	illumina_aqm(dataset,target_file,project)
}

####if data is already normalised need to rename it as normalised_exp and copy it to the ominer-results/project/norm directory
if (dataType == "normalised") {
	normalised_path = paste("../../../html/onlinetool/temp/",project,"normalised.txt",sep="/")
	file.copy(normalised_path,normdir)
}

if (survival == "1") {
	survivaldir = paste(dir,"survival",sep="/")
dir.create(survivaldir)
source("KM_surv_plots.R")
survival_os_plot(target_file,project)
cat ("Running survival analysis",file=errorFile,sep="\n",append=TRUE)
	}	






norm_matrix = "normalised.txt"
if (combat == "1") {
	source("combat_norm.R")
	run_combat(target_file,norm_matrix,project)
	cat ("Running batch effect correction using COMBAT",file=errorFile,sep="\n",append=TRUE)
}
source("filtering_illum.R")
filtering(filter,filterval,combat,project)
cat ("Filtereing normalised expression matrix",file=errorFile,sep="\n",append=TRUE)

source("run_limma_illumina.R")
data = "filtered_data.txt"

if (dataType == "unnormalised") {
print ("I AM UNNORMALIZED")
if (revised == "0") {
run_limma(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange,platform,project)
cat ("Running differential expression analysis",file=errorFile,sep="\n",append=TRUE)
}
if (revised == "1") {
	run_limma(target_qc_file,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange,platform,project)
	cat ("Running differential expression analysis",file=errorFile,sep="\n",append=TRUE)
}
}


if (dataType == "normalised"){
    run_limma(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange,platform,project)
    cat ("Running differential expression analysis",file=errorFile,sep="\n",append=TRUE)
    
}

testObject <- function(object)
	{
		exists(as.character(substitute(object)))
	}
        
if(testObject(venn)) {
    

source("generate_venn.R")

generate_venn("ebA_allresults.txt",project)
cat ("Generating venn diagrams",file=errorFile,sep="\n",append=TRUE)
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
filtered_x <- paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_Filtered.txt",sep="")
if (file.exists(filtered_x)) {
	go_analysis = "1"
	godir = paste(dir,"GO",sep="/")
dir.create(godir)
}

if (go_analysis == "1") {
	
	source("illumina_GO_analysis.R")
	illum_GO(project,comp)
	cat ("Running Gene Ontology analysis",file=errorFile,sep="\n",append=TRUE)
	}
	
}
}



source("cluster_illum_code.R")
if (revised == "0") {
#####number of clusters k <- lev (this is the number of biological groups)
k <- length(lev)
print ("This is the number of clusters")
print(k)	
cluster_all(target_file,k,project)
cat ("Running clustering analysis",file=errorFile,sep="\n",append=TRUE)
}
if (revised == "1") {
	k <- length(lev)
print ("This is the number of clusters")
print(k)
	cluster_all(target_qc_file,k,project)
	}
	
comps <- read.table(comp_file, sep = "\t", as.is = T, header = F, 
        strip.white = T)
        for (i in 1:length(comps)) {
                print(comps[i])
        fg <- strsplit(comps[,i],"=")
        ci <- fg[[i]][1]
####need to do this for all comparisons
        #filteredx <- paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"Filtered.txt",sep="")
#sigGene=rownames(filteredx)

comps <- read.table(comp_file, sep = "\t", as.is = T, header = F, 
        strip.white = T)
        for (i in 1:length(comps)) {
                print(comps[i])
        fg <- strsplit(comps[,i],"=")
        ci <- fg[[i]][1]
####need to do this for all comparisons
        #filteredx <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"Filtered.txt",sep=""),row.names = 1, header = T, sep = "\t") 
#sigGene=rownames(filteredx)
filtered_x <- paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_Filtered.txt",sep="")
if (file.exists(filtered_x)) {
filtered_x <- read.table(paste("ominer_results/",project,"/","DifferentialExpression/",ci,"_Filtered.txt",sep=""),header = T, sep = "\t",as.is=T)
length_AGI <- length(filtered_x$AGI)
if (length_AGI > 2) {
source("heatmap_create_TRIMMED_ILLUMINA.R")
cat ("Generating heatmap",file=errorFile,sep="\n",append=TRUE)

if (revised == "0") {
heatmap_generate(target_file,project,comp_file)
}
if (revised == "1") {
target_file = target_qc_file
heatmap_generate(target_file,project,comp_file)
}
}
}
}
}	
	
	
	
resultdir = "/var/www/cgi-bin/onlinetool/version_2/ominer_results/"
setwd(resultdir)
zipcommand=paste("zip -r ",project,".zip ",project,sep="")
print(zipcommand)
system(zipcommand)
file.rename(paste(resultdir,project,".zip",sep=""),paste(resultdir,project,"/",project,".zip",sep=""))


command="/var/www/cgi-bin/onlinetool/mail.pl"
url="Your results are now available at : http://o-miner.org/cgi-bin/onlinetool/version_2/exprOutputView_dec_2013.pl?"
mailcommand=paste(command," ",email," ",project," \"",url,project,"\"",sep="")
#mailcommand=paste(url,project,sep="")
system(mailcommand)
#source("heatmap_create.R")

#DE_list = "NormvstumFiltered.txt"  ###if user picks comparison DE_list = paste(comp[i],"_all_DE.txt") - comp[i] realted to file comp.txt containing comparisons
#heatmap_generate(target,DE_list,project)

