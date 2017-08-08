##########################################################################################
#
#File name: miRNA_commands.R
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
#Description:Master R script for running analysis of array based data from Affymetrix miRNA platforms miRNA 1.0 and miRNA 3.0. This R script calls several different functions. In order to run this pipeline you need to have prepared two separate files:
#a target file and a comparisons file - examples of these files can be downloaded with this script. The script can be run from two different data entry points these are: (1) Raw CEL files and (2) Normalised expression
#matrix.A target file for straightforward analysis will have three columns (1) Name, (2) FileName and (3) Target. Additional columns are required in the target file if analysis is run as paired (target_paired.txt), if replicates (target_replicates.txt) are
#present and if combat is to be used in a meta-analysis (target_combat). Please note that the target file is always refered to as 'target.txt' regardless of whether analysis is run as paired or whether replicates are present or if combat is used. 
#Results are written out at each stage of the analysis to a subdirectory within the project folder that is created within the 'ominer_results' folder.
#pre-requisites: Users must have the following R packages installed on their system: (1) oligo, (2) biomaRt, (3) limma, (4) survival, (5) sva, (6) gostats, (7) ggplots (8) FC14.plotting.lib. 
#Before running this pipeline you will need to have created a folder named 'ominer_results'.
#Code dependencies:
#Venn.R
#Code.R
#bcctb.utils.R
#The code dependencies above are distributed with this master script.
#The command line arguments needed to run this script are:
#R --file=/scratch/ominer/transcriptomics/Affymetrix/miRNA_commands.R --args  dataType=CEL folder=/var/www/html/onlinetool/temp/a.sangaralingam-miRNA_test_submission
#project=a.sangaralingam-miRNA_test_submission target=/var/www/html/onlinetool/temp/a.sangaralingam-miRNA_test_submission/target.txt comp=/var/www/html/onlinetool/temp/a.sangaralingam-miRNA_test_submission/comp.txt normalisation=RMA
#filter=sd filterval=40 adjust=BH pvalue=0.05 foldchange=2.0 limmamethod=separate analysis=unpaired replicates=no combat=0 aqm=0 gostats=1 survival = 0
#Example of shell script containing command line arguments: a.sangaralingam-miRNA_test_submission.sh
#Description of command line arguments:
#platform: Affymetrix miRNA platform from which data is being analysed i.e. for Affymetrix Genechip miRNA 1.0 use "mirna2" and for Affymetrix Genechip miRNA 3.0 use "mirna2"
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
######################################################################################################################################################################################################################################################################

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

#setwd("/scratch/ominer/transcriptomics/miRNA/")
setwd("/var/www/cgi-bin/onlinetool/version_2/")
dir = paste("ominer_results",project,sep="/")
dir.create (dir)
normdir = paste(dir,"norm",sep="/")
dir.create(normdir)

errorFile = paste(dir,"/run_log.txt",sep="")
file.create(errorFile)

target_path = "../../../html/onlinetool/temp"
target_file = paste(target_path,project,"target.txt",sep="/")
comp_path = "../../../html/onlinetool/temp"
comp_file =  paste(comp_path,project,"comp.txt",sep= "/")
mydir=getwd()

library(oligo);
if (platform =="mirna2") {
library(pd.mirna.3.0);
}
if (platform == "mirna1") {
library(pd.mirna.1.0);
}
library(biomaRt);
library(limma);
qcdir = paste(dir,"QC",sep="/")
dir.create(qcdir)
file.copy(target_file,qcdir)

# read in annot.data
pd <- read.table(target_file, sep="\t", header=T);
lev <- unique(pd$Target)
if (dataType == "CEL") {
row.names(pd) <- pd[["FileName"]];
setwd(folder)
# read in .cel files (from oligo) -nb CDF not available for 3.0)
affybatch <- read.celfiles( row.names(pd) );

setwd(mydir)
cat ("Read in CEL files",file=errorFile,sep="\n",append=TRUE)
source("miRNA_ominer.R")
#norm = "gcrma"
run_miRNA(target_file,project)

#if (survival == "yes") {
#source("KM_surv_plots.R")
#survival_os_plot("target_file,project")
	
#}	
}
if (dataType == "normalised") {
	####copy normalised data matrix to correct place
	norm_path = "../../../html/onlinetool/temp"
normalisation =  paste(norm_path,project,"normalised.txt",sep= "/")
	        file.copy(normalisation,normdir)
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
	cat ("Running batch effect correction COMBAT",file=errorFile,sep="\n",append=TRUE)
}


source("filtering.R")
filtering(filter,filterval,combat,project)
cat ("Filtering normalised expression matrix",file=errorFile,sep="\n",append=TRUE)
if (dataType == "filtered") {
	####copy filtered matrix to correct place
	}



DEdir = paste(dir,"DifferentialExpression",sep="/")
dir.create(DEdir)
source("miRNA_limma_annotation.R")
data = "filtered_data.txt"
run_limma(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange,project)
cat ("Running differential expression analysis",file=errorFile,sep="\n",append=TRUE)
clusterdir = paste(dir,"cluster",sep="/")
dir.create(clusterdir)


testObject <- function(object)
	{
		exists(as.character(substitute(object)))
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
	source("miRNA_GO.R")
	cat ("Running Gene Ontology analysis",file=errorFile,sep="\n",append=TRUE)
	miRNA_GO(project,comp_file)

}
}
}

#comps <- read.table(comp_file, sep = "\t", as.is = T, header = F, 
        #strip.white = T)
        #for (i in 1:length(comps)) {
                #print(comps[i])
        #fg <- strsplit(comps[,i],"=")
        #ci <- fg[[i]][1]
####need to do this for all comparisons
        #filteredx <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"Filtered.txt",sep=""),row.names = 1, header = T, sep = "\t") 
#sigGene=rownames(filteredx)
#if (platform == "st1") {
#filtered_x <- paste("ominer_results/",project,"/","transcript/DifferentialExpression/",ci,"_Filtered.txt",sep="")
#}
#else {
#filtered_x <- paste("ominer_results/",project,"/","DifferentialExpression/",ci,"_Filtered.txt",sep="")
#}
#if (file.exists(filtered_x)) {
#source("heatmap_create.R")
#cat ("Generating heatmap",file=errorFile,sep="\n",append=TRUE)


#heatmap_generate(target_file,project,comp_file)
	

#}
#}



        
if(testObject(venn)) {
    

source("generate_venn.R")
generate_venn("ebA_allresults.txt",project)
cat ("Generating Venn diagrams",file=errorFile,sep="\n",append=TRUE)
}

source("cluster_code.R")
k <- length(lev)
cluster_all(target_file,k,project)
cat ("Clustering data",file=errorFile,sep="\n",append=TRUE)
#cluster_all(target_file,"2",project)


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
length_AGI <- length(filtered_x$ID)
if (length_AGI > 2) {
source("heatmap_create_TRIMMED_MIRNA.R")
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



command="/var/www/cgi-bin/onlinetool/version_2/mail.pl"
url="Your results are now available at : http://o-miner.org/cgi-bin/onlinetool/version_2/exprOutputView_dec_2013.pl?"

mailcommand=paste(command," ",email," ",project," \"",url,project,"\"",sep="")
#mailcommand=paste(command," ",email," ",url,project,"\"",sep="")
system(mailcommand)





#source("heatmap_create.R")

#de_genes = "TumourvsNormal_annotated_all_DE.txt"  ###if user picks comparison DE_list = paste(comp[i],"_all_DE.txt") - comp[i] realted to file comp.txt containing comparisons
#heatmap_generate(target_file,de_genes,project)
