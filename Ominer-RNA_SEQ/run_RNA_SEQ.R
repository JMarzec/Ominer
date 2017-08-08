##########################################################################################
#
#File name: run_RNA_SEQ.R
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
#Description:Master R script for the post-processing of RNA-Seq data. This R script calls several different functions. In order to run this pipeline you need to have prepared two separate files:
#a target file and a comparisons file - examples of these files can be downloaded with this script. The script can be run from two different data entry points these are: (1) Raw CEL files and (2) Normalised expression
#matrix.A target file for straightforward analysis will have three columns (1) Name, (2) FileName and (3) Target. Additional columns are required in the target file if analysis is run as paired (target_paired.txt), if replicates (target_replicates.txt) are
#present and if combat is to be used in a meta-analysis (target_combat). Please note that the target file is always refered to as 'target.txt' regardless of whether analysis is run as paired or whether replicates are present or if combat is used. 
#Results are written out at each stage of the analysis to a subdirectory within the project folder that is created within the 'ominer_results' folder.
#pre-requisites: Users must have the following R packages installed on their system: (1) sva, (2) edgeR, (3) limma, (4) goseq
#(5) affyPLM, (6) biomart, (7) annotate,  (8) gplots, (12) FC14.plotting.lib and  (13) survival. 
#Before running this pipeline you will need to have created a folder named 'ominer_results'.
#Code dependencies:
#Venn.R
#Code.R
#bcctb.utils.R
#The code dependencies above are distributed with this master script.
#The command line arguments needed to run this script are: 
#R--file-/scratch/ominer/transcriptomics/Affymetrix/run_RNA_SEQ.R --args dataType=unnormalised folder=/var/www/html/onlinetool/temp/a.sangaralingam-RNA_SEQ_O
MINER_TEST-8089 project=a.sangaralingam-RNA_SEQ_OMINER_TEST-8089 target=/var/www/html/onlinetool/temp/a.sangaralingam-RNA_SEQ_OMINER_TEST-8089/target.txt comp=/var/www/html/onlinetool/temp/a.sangaralingam-RNA
_SEQ_OMINER_TEST-8089/comp.txt normalisation=0 filter=sd filterval=40 adjust=BH pvalue=0.05 foldchange=2.0 diffmethod=edge limmamethod= analysis=unpaired replicates=no combat=0 survival=0 gostats=1
#Example of shell script containing command line arguments: a.sangaralingam-RNA_SEQ_OMINER_TEST-8089.sh
#Description of command line arguments:
#dataType: Whether data analysed is normalised or unnormalised i.e. normalised unnormalised
#folder: full path to where the input data is found
#project: name given to the analysis project
#target: full path to the target file
#comp: full path to the comparison file
#filter: method used to filter data e.g. sd, iqr
#filterval: % value to filter the normalised expression matrix i.e to reatin the topmost 40% filterval = 40
#adjust: method to adjust the false discovery rate (FDR) within limma e.g BH (Benjamini-Hochberg) or BY
#pvalue: adjusted p-value cutoff to report differentially expressed genes
#foldchange: log2 foldchange cutoff to report differentially expressed genes
#diffmethod: method to use for differential expression analysis - there is a choice of two different methods: LIMMA & edgeR - arguments are: limma or edge
#limmamethod:Method used within limma when comparing groups, e.g separate
#analysis: whether an unpaired or a paired analysis is executed
#replicates: yes/no indicating whether technical replicates are present within the dataset
#combat:0/1 if the combat algorithm is to be applied to the data
#survival: 0/1 if survival analysis is to be executed
#gostats: 0/1 if gostats analysis is to be executed
#Description of output: Data is output to the ominer_results directory, a folder is automatically created within this directory with your project name. Subfolders are created within this directory such as (1) QC, (2) norm, (3) Differential Expression
############################################################################################################################################################# Parse the command line options
for (e in commandArgs()) {
	ta = strsplit(e,"=",fixed=TRUE)
	ta
#	if(! is.na(ta[[1]][2])) {
		temp = ta[[1]][2]
		assign(ta[[1]][1],temp)
		cat("assigned ",ta[[1]][1]," the value of |",temp,"|\n")
#	} 
}

#setwd("/scratch/ominer/transcriptomics/RNA_seq/")
setwd("/var/www/cgi-bin/onlinetool/version_2/")

dir = paste("ominer_results",project,sep="/")

dir.create (dir)   ###make a directory called project under OMINER-OUTPUT each time that this script is run
qcdir = paste(dir,"QC",sep="/")
dir.create(qcdir)
normdir = paste(dir,"norm",sep="/")
dir.create(normdir)
godir = paste(dir,"GO",sep="/")
dir.create(godir)
target_path = "../../../html/onlinetool/temp"
#target_path = "/var/www/cgi-bin/onlinetool/version_2"
#target_path = "/scratch/ominer/transcriptomics/RNA_seq/"
target_file = paste(target_path,project,"target.txt",sep="/")
#target_file = paste(target_path,target,sep="/")
#comp_path = "/var/www/cgi-bin/onlinetool/version_2/"
comp_path = "../../../html/onlinetool/temp"

errorFile = paste(dir,"/run_log.txt",sep="")
file.create(errorFile)

#comp_path = "/scratch/ominer/transcriptomics/RNA_seq/"
comp_file =  paste(comp_path,project,"comp.txt",sep= "/")
counts_path = 
counts_path = "../../../html/onlinetool/temp"
#counts_path = "/scratch/ominer/transcriptomics/RNA_seq/"
counts_file = paste(counts_path,project,"normalised.txt",sep="/")   ####will probably want to change this to target_path/counts_matrix.txt?

file.copy(counts_file,normdir)
file.copy(target_file,qcdir)
data = "normalised.txt"
mydir=getwd()
changed_file = "0"

t_file <- read.table(target_file, header = T, sep = "\t",as.is=TRUE,row.names=NULL)

             lev <- unique(t_file$Target)
n_matrix <- read.table(counts_file,header=T,sep="\t",as.is=TRUE)

total_length <- ncol(n_matrix)
target_length <- nrow(t_file)

print ("This is the length of the target file")
print (target_length)
print ("This are the dimensions of the counts data")
print (n_matrix)
if (total_length > target_length) {
	changed_file = "1"
}

if (changed_file == "1") {
detection_names = NULL
n <- length(t_file$Name)
for (i in seq(from=1,to=n,by=1)) {
	
	#wanted <- paste("\\",t_file$Name[i],"\\b",sep="")
wanted <- paste("^",t_file$Name[i],"$",sep="")
detection_names <- c(detection_names,grep(wanted,colnames(n_matrix)))
}


 new_dat <- subset(n_matrix[detection_names])

gene_id <- n_matrix[,1]    ##### Indexing like this allows counts matrix even with ambigously named first column(s) to work
new_matrix <- cbind(gene_id,new_dat)

new_size <- ncol(new_matrix)
print ("This is the size of the new matrix")
print (new_size)
write.table(new_matrix,paste("ominer_results/",project,"/norm","/refined.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
}

if (changed_file == "0") {
		dataset = paste("ominer_results/",project,"/norm","/normalised.txt",sep="")
		data <- "normalised.txt"
		}
		if (changed_file == "1") {
			dataset = paste("ominer_results/",project,"/norm","/refined.txt",sep="")
			data <- "refined.txt"
		}


if (combat == "1") {
	source("combat_norm.R")
	norm_matrix <- data
	run_combat(target_file,norm_matrix,project)
	cat ("Running batch effect correction using COMBAT",file=errorFile,sep="\n",append=TRUE)
}

source("filtering_RNAseq.R")
cat ("filtering normalised expression matrix",file=errorFile,sep="\n",append=TRUE)

filtering(filter,filterval,combat,project,data)
data <- "filtered_data.txt"

DEdir = paste(dir,"DifferentialExpression",sep="/")
dir.create(DEdir)
if (diffmethod == "edge") {
	
	source("edgeR_function.R")
	run_edgeR (data,target,analysis,comp_file,replicates,adjust,pvalue,project) 
	
	source("annotate_edgeR.R")
	annotated_edgeR(comp_file,project)
	cat ("Running differential expression analysis using edgeR",file=errorFile,sep="\n",append=TRUE)
}

if (diffmethod == "limma") {
	source("run_limma_affy_RNAseq_modS.R")
	run_limma(target,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange,platform,project,dataType)
	cat ("Running differential expression analysis using limma",file=errorFile,sep="\n",append=TRUE)
}

if (diffmethod == "DESeq") {
	#source("DESEQ_function.R")
}

if (diffmethod == "DEXSeq") {
	#source("DEXseq.R")
	
}
testObject <- function(object)
	{
		exists(as.character(substitute(object)))
	}
	
		if(testObject(gostats)) {


	if (diffmethod == "edge") {
	source("goseq.R")
	run_rna_seq_GO(project,comp_file,diffmethod,data)
	cat ("Running Gene Ontology analysis",file=errorFile,sep="\n",append=TRUE)
########need to find the GO-seq function that I wrote for this 	
}
if (diffmethod == "limma") {
	source("goseq.R")
	run_rna_seq_GO(project,comp_file,diffmethod,data)
	cat ("Running Gene ontology analysis",file=errorFile,sep="\n",append=TRUE)
	
}



}
if(testObject(venn)) {
    


source("generate_venn.R")



generate_venn("ebA_allresults.txt",project)
cat ("Generating venn diagram",file=errorFile,sep="\n",append=TRUE)

}


clusterdir = paste(dir,"cluster",sep="/")
dir.create(clusterdir)
source("cluster_code_mod_RNAseq.R")
cat ("Clustering data",file=errorFile,sep="\n",append=TRUE)
#####number of clusters k <- lev (this is the number of biological groups)
k <- length(lev)
print ("This is the number of clusters")
print(k)
cluster_all(target_file,k,project)



#source("cluster_RNA-seq.R")

#cluster_all(target,"2",project,counts_file)
cat ("Running clustering analysis",file=errorFile,sep="\n",append=TRUE)
######Do heatmaps - from the normalized count matrix 
#if (limmamethod == "limma") {
#	if (venn == 1)
	#source("generate_venn.R")
	#generate_venn("decideTestsSummary.txt",project)
#}

if (survival == "1") {
	survivaldir = paste(dir,"survival",sep="/")
dir.create(survivaldir)
source("KM_surv_plots.R")
survival_os_plot(target_file,project)
cat ("Running survival analysis",file=errorFile,sep="\n",append=TRUE)
	
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
length_AGI <- length(filtered_x$ensembl_gene_id)
if (length_AGI > 2) {
source("heatmap_create_TRIMMED_RNA_SEQ.R")
cat ("Generating heatmap",file=errorFile,sep="\n",append=TRUE)


heatmap_generate(target_file,project,comp_file)


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




#command="/var/www/cgi-bin/onlinetool/mail.pl"
#url="Your results are now available at : http://o-miner.org/cgi-bin/onlinetool/exprOutputView_dec_2013.pl?"

mailcommand=paste(command," ",email," ",project," \"",url,project,"\"",sep="")
##mailcommand=paste(command," ",email," ",url,project,"\"",sep="")
system(mailcommand)

