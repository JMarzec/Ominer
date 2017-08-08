##########################################################################################
#
#File name: EA_run.R
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
#a target file and a comparisons file - examples of these files can be downloaded with this script. The script can be run from Raw CEL files only.
#A target file for straightforward analysis will have three columns (1) Name, (2) FileName and (3) Target. Additional columns are required in the target file if analysis is run as paired (target_paired.txt), if replicates (target_replicates.txt) are
#present and if combat is to be used in a meta-analysis (target_combat). Please note that the target file is always refered to as 'target.txt' regardless of whether analysis is run as paired or whether replicates are present or if combat is used. 
#Results are written out at each stage of the analysis to a subdirectory within the project folder that is created within the 'ominer_results' folder.
#pre-requisites: Users must have the following R packages installed on their system: (1) Aroma.affymetrix, (2) affy, (3) limma, (4) sva, (5) survival and (6) Gostats 
#Before running this pipeline you will need to have created a folder named 'ominer_results'.
#Code dependencies:
#Venn.R
#Code.R
#The code dependencies above are distributed with this master script.
#The command line arguments needed to run this script are:
#R--file-/scratch/ominer/transcriptomics/Affymetrix/EA_run.R platform=st1 dataType=CEL folder=/var/www/html/onlinetool/temp/a.sangaralingam-TEST_EXON_ARRAY project=a.sangaralingam-TEST_EXON_ARRAY target=/var/www/html/onlinetool/temp/a.sangaralingam-TEST_EXON_ARRAY/target.txt
#comp=/var/www/html/onlinetool/temp/a.sangaralingam-TEST_EXON_ARRAY/comp.txt normalisation=0 filter=sd filterval=40 adjust=BH pvalue=0.05 foldchange=0.5 limmamethod=separate analysis=unpaired replicates=no combat=0 aqm=0 
#survival=0 estimate=0 aqm=0 gostats=1
#Example of shell script containing command line arguments: EA_example_run.sh
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
#gostats: 0/1 if gostats analysis is to be executed
#Description of output: Data is output to the ominer_results directory, a folder is automatically created within this directory with your project name. Subfolders are created within this directory such as (1) exon, (2) splicing, (3) transcript and (4) survival
#within the exon,splicing and transcript folders, subfolders (A) cluster, (B) DifferentialExpression and (C) norm are found.
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

#setwd("/scratch/ominer/transcriptomics/exon_array/")
setwd("/var/www/cgi-bin/onlinetool/version_2/")
mydir = getwd()

##aromadir is the path to the directory where aroma annotation directory (annotationData) and the raw data directory (rawData) is located

aromadir = "/var/onlinetool/version_2/Aroma/"


##create all output directories
dir = paste("ominer_results","/",project,sep="")
dir.create (dir)   ###make a directory called project under OMINER-OUTPUT each time that this script is run
trdir = paste(dir,"transcript",sep="/")
dir.create(trdir)
exdir = paste(dir,"exon",sep="/")
dir.create(exdir)
spdir = paste(dir,"splicing",sep="/")
dir.create(spdir)
qcdir = paste(trdir,"QC",sep="/")
dir.create(qcdir)
trnormdir = paste(trdir,"norm",sep="/")
dir.create(trnormdir)
ex_normdir = paste(exdir,"norm",sep="/")
dir.create(ex_normdir)
sp_normdir = paste(spdir,"norm",sep="/")
dir.create(sp_normdir)
sp_de = paste(spdir,"DifferentialExpression",sep="/")
dir.create(sp_de)
trdir_de = paste(trdir,"DifferentialExpression",sep="/")
dir.create(trdir_de)
exdedir = paste(exdir,"DifferentialExpression",sep="/")
dir.create(exdedir)
#trdedir = paste(dedir,"transcripts",sep="/")
#spdir = paste(dedir,"splicing",sep="/")
#dir.create(spdir)
#dir.create(trdedir)
#dir.create(exdedir)
tr_clusterdir = paste(trdir,"cluster",sep="/")
dir.create(tr_clusterdir)
ex_exclusterdir = paste(exdir,"cluster",sep="/")
dir.create(ex_exclusterdir)
tr_godir = paste(trdir, "GO",sep="/")
dir.create(tr_godir)
#dir.create(trclusterdir)
target_path = "../../../html/onlinetool/temp"
target_file = paste(target_path,project,"target.txt",sep="/")
comp_path = "../../../html/onlinetool/temp"
comp_file =  paste(comp_path,project,"comp.txt",sep= "/")
file.copy(target_file,qcdir)




platform = "HuEx-1_0-st-v2"
chipType = platform
		
# Set up directory structure for pipeline based on input type
	if(dataType == "CEL") {
				dir.create(paste(aromadir,"rawData/",project,"/",sep=""))
		dir.create(paste(aromadir,"rawData/",project,"/",platform,sep=""))
		
		print(paste("Reading CEL files from :",folder))
		
#read the chip names and patient names from the target file & copy over to correct Aroma rawData directory
		
			
			
			 pd <- read.table(target_file, header = T, row.name = "Name", 
        sep = "\t")
					
			name=pd$FileName
			print(name)
#Copy file over from user directory to aroma directory
			for (i in 1:length(name))
			{
				print(paste(folder,"/",name[i],sep=""))
				print(paste(aromadir,"rawData/",project,"/",platform,"/",name[i],".CEL",sep=""))
				
				file.copy(paste(folder,"/",name[i],sep=""), paste(aromadir,"rawData/",project,"/",platform,"/",name[i],sep=""))
				#print(chip[i])
			}
		
	
		print("..Done reading CEL files from folder & copying over to aroma structure")	
		#file.copy(paste(folder,"/chip.txt",sep=""),paste(resultdir,project,"/chip.txt",sep=""))
            }    


source("exon_qc.R")
qc(chipType,project,aromadir,mydir)
errorFile = paste(dir,"/run_log.txt",sep="")
file.create(errorFile)
cat ("Running quality control steps",file=errorFile,sep="\n",append=TRUE)
if (survival == "1") {
	survivaldir = paste(dir,"survival",sep="/")
dir.create(survivaldir)
source("KM_surv_plots.R")
survival_os_plot(target_file,project)
cat ("Running survival analysis",file=errorFile,sep="\n",append=TRUE)
	}	
if (combat == "yes") {
	source("combat_norm.R")
	run_combat(target_file,norm_matrix,project)
	cat ("Running batch effect correction using COMBAT",file=errorFile,sep="\n",append=TRUE)
}

source("run_limma_exon.R")

data = "normalised.txt"
run_limma_exon(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange)
cat ("Running differential expression analysis for exons",file=errorFile,sep="\n",append=TRUE)

source("run_limma_transcript.R")
#data = "trfit.txt"
run_limma_transcript(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange)
cat ("Running differential expression analysis for transcripts",file=errorFile,sep="\n",append=TRUE)

source("run_limma_firma.R")
#data = "FmFit.txt"
run_limma_firma(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange)
cat ("Running differential splicing analysis",file=errorFile,sep="\n",append=TRUE)





######run differential splicing analysis
#######modify the normalised exon and transcript matrices so apling index can be calculated
#readdir <- paste("ominer_results",project,sep="/")
#data <- "ExFittest.txt"
#ex_Normdata <- read.table(paste(readdir,"/","norm/",data,sep=""), sep = "\t", as.is = T, header = T, strip.white = TRUE)
#ncolumns <- ncol(ex_Normdata)
#newdata <- ex_Normdata[6:ncolumns]
#ll_newdata <- cbind(ex_Normdata$unitName,newdata)
#write.table(all_newdata, paste(readdir,"/norm","/newexFit.txt", sep=""),
 #       sep = "\t", row.names = FALSE, quote = FALSE)
#tr_data <- "trfittest.txt"

#tr_Normdata <- read.table(paste(readdir,"/","norm/",tr_data,sep=""), sep = "\t", as.is = T, header = T, strip.white = TRUE)
#ncolumns <- ncol(tr_Normdata)
#tr_newdata <- tr_Normdata[6:ncolumns]
#alltr_newdata <- cbind(tr_Normdata$unitName,tr_newdata)
#write.table(alltr_newdata, paste(readdir,"/norm","/newetrFit.txt", sep=""),
 #       sep = "\t", row.names = FALSE, quote = FALSE)


#cmd <- paste("perl NI_calc_mod.pl",project);
#system(cmd)

#ni_file <- paste("ominer_results",project,"norm","ni.txt",sep="/")
#col_names_needed <- colnames(tr_Normdata)
#cl_names <- col_names_needed[6:ncolumns]
#ni_data <- read.table(ni_file,sep = "\t", as.is = T, header = F, strip.white = TRUE)
#new_ni <- ni_data[1:10]
#colnames(ni_data) <- c("psid","trid",cl_names)
#new_ni <- ni_data[1:10]
######calculate splicing index
#	pd<-read.table(target_file,header=T, row.name="Name",sep="\t") # read in targets file
#	ncols_ni <- ncol(new_ni)
#	names <- new_ni[1:2]
#	colnames(names) <- c("psid","trid")
#	data_ni <- new_ni[3:ncols_ni]
#	colnames(data_ni) <-pd$Target
#	all <- cbind(names,data_ni)
#	comps<-read.table(comp_file,sep="\t",as.is=T,header=F,strip.white=T)
#for (i in 1:length(comps)) {
 #               print(comps[i])
  #      fg <- strsplit(comps[,i],"=")
   #     ci <- fg[[i]][1]
#print(ci)

#split_ci <- strsplit(ci,"vs")
#type_1 <- split_ci[[1]] [1]
#type_2 <- split_ci[[1]] [2]
#detection_names_type1 <- NULL
#detection_names_type1 <- c(detection_names_type1,grep(type_1,colnames(all)))
#new_dat_type1 <- subset(all[detection_names_type1])
#detection_names_type2 <- NULL
#detection_names_type2 <- c(detection_names_type1,grep(type_2,colnames(all)))
#new_dat_type2 <- subset(all[detection_names_type2])






#gene_id <- all[1]  ###transcript_cluster_ID
#exon_id <- all[2]  ###probeset_id
#####for each row in above matrices calculate median value
#control_group_mat <- as.matrix(new_dat_type1)
#tumor_group <- as.matrix(new_dat_type2)
#mean_control <- apply(control_group_mat,1,median)
#mean_PDAC <- apply(tumor_group,1,median)
#mc <- as.numeric(mean_control)
#mpd <- as.numeric(mean_PDAC)
#SI <- log2(t(t(mc)/mpd)) ###need to subtract for my data as already log2 transformed 
#SI <- log2(mc-mpd)
#all_data <- cbind(gene_id,exon_id,SI)
#colnames(all_data) <- c("transcript_cluster_id","exon_id","SI")

#
#write.table(all_data, paste(readdir,"/DifferentialExpression/",comps[i],"_SI.txt", 
 #           sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
#}
#write.table(all_data, "all_SI.txt",sep="\t",quote=F,row.names=T)






#source("get_annotation_exon.R")
#get_annotation(project,comp)

if (gostats == "1") {
	
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
filtered_x <- paste("ominer_results/",project,"/","transcript/DifferentialExpression/",ci,"_Filtered.txt",sep="")
if (file.exists(filtered_x)) {
	go_analysis = "1"
	godir = paste(trdir,"GO",sep="/")
dir.create(godir)
}

if (go_analysis == "1") {
	
	source("Exarray_GO.R")
	Exon_GO(project,comp)
	cat ("Running Gene Ontology analysis",file=errorFile,sep="\n",append=TRUE)
	}
	
}
}


source("cluster_code_exon.R")
cluster_all_Exon(target_file,"2",project)
cat ("Running clustering for exons",file=errorFile,sep="\n",append=TRUE)
source("cluster_code_transcript.R")
cluster_all_tr(target_file,"2",project)
cat ("Running clustering for transcripts",file=errorFile,sep="\n",append=TRUE)



source("heatmap_create.R")
#de_genes = "PDACVSCONTROL_TopTable.txt"
#heatmap_generate(target,de_genes,data,project)
#source("get_annotation_exon.R")
#file = "PDACVSCONTROLTopTable_all_annotated.txt"#
#get_annotation(file)

source("generate_ea_venn.R")
generate_venn("ebA_allresults.txt",project)
cat ("Generating Venn diagrams",file=errorFile,sep="\n",append=TRUE)

resultdir = "/var/www/cgi-bin/onlinetool/version_2/ominer_results/"
setwd(resultdir)
zipcommand=paste("zip -r ",project,".zip ",project,sep="")
print(zipcommand)
system(zipcommand)
file.rename(paste(resultdir,project,".zip",sep=""),paste(resultdir,project,"/",project,".zip",sep=""))

command="/var/www/cgi-bin/onlinetool/version_2/mail.pl"
url="Your results are now available at : http://o-miner.org/cgi-bin/onlinetool/version_2/exprOutputView_dec_2013.pl?"

mailcommand=paste(command," ",email," ",project," \"",url,project,"\"",sep="")

system(mailcommand)



