############################################################################################################
#
#File name CNV_shell.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Master R script for the analysis of raw copy number data from Affymetrix platforms and processed copy number data from any platform. This R script calls several functions:
#Pre-requisites
#A number of R packages need to be pre-installed prior to using this pipeline - these are:
#aroma.affymetrix
#DNAcopy
#ASCAT - if ASCAT is to be used the the ASCAT code needs to be downloaded 


#In order to run this script a number of arguments are required these are:
#Normaldir = This is the name of the directory containing "Normal" samples i.e. these are samples taken from blood etc.
#threshold = threshold used during segmentation to call a gain or a loss
#project = this is the name given to the analysis project
#analysisMethod = This can be either ASCAT or CBS
#folder = this is the path to where the input data is found
#baseline = user (this is when a user defined baseline is used) - or hapmap can be used as this argument 
#email = this argument is not relevant if being run from the command line
#cdfarray1 = name of the affymetrix latform id raw data is used
#patient = number of samples found in dataset
#platform = the Affymetrix platform that is being used (1) argument - six if Affy GenomeWideSNP6 is used, (2) 100 for 100K platform, (3) 500 for the 500K platform, (4) 250sty for the 250STY chip, (5) 250nsp for the 250NSP chip
#cdfarray2 = this argument is used if the user is analysing a platform that is made up of two chips e.g. 500k which is made up of 250k nsp and 250k sty. If user is analysing data from a platform that is not made up of two chips
#then this argument is left empty.
#snp_number = this is to set the minimum number of snps the user would like to view in each region.
#analysisType = this is set to either paired/unpaired
#hap = this argument is only used if the user chooses to compare the samples to a hapmap dataset - the id of the HapMap dataset is then used as the argument
#targetFile = name of the .txt file this takes the same format as the target file for the expression pipeline - two columns are present if a platform containg two chips are being analysed.
#pair_array1 = name of .txt file and full directory path to this file containing the pair information for the dataset to be analysed.
#type = this argument refers to the entry point of the pipeline (1) if raw data files are being processed this is argument should be CEL, (2) log2ratios are used argument should be log2ratios, (3) segmented values then argumnet = segmented
#(4) smoothed values - argument = smoothed.
#mcr = type of Minimal Common Region (MCR) method to be used on the data e.g. CGHregions
#Normal = this gives the full path to the .txt file containing information regarding the normal samples.
#Please note before running this pipeline you will need to have created a folder named "ominer_results"




#analysisType_method in webpage refers the whether user wanted to use ASCAT or CBS method == analysisType as arguments here equal to paired/unpaired


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




Cancer=project

######From previous version, set up files   ##############
rootdir="/var/www/cgi-bin/onlinetool/version_2/"
setwd(rootdir)
aromadir = "/var/onlinetool/version_2/Aroma/"

resultdir = paste("ominer_results",project,sep="/")

dir.create(resultdir)
errorfile=paste(resultdir,"/",Cancer,"/run_log.txt",sep="")
####Set the minimum number of SNPs to call a region
SNP_nb=as.numeric(snp_number)
#####Set the threshold for which gains and losses are called
threshold=as.numeric(threshold)



#####Set each of the directories for each of the subfolders of the analysis 
dir.create(paste(resultdir,"/Regions/",sep=""))# when automating
dir.create(paste(resultdir,"/Regions/","segments",sep=""))# when automating
dir.create(paste(resultdir,"/Regions/","loh",sep=""))# when automating
dir.create(paste(resultdir,"/cghweb/",sep=""))
cghwebdir = paste(resultdir,"/cghweb/",sep="")
	dir.create(paste(resultdir,"/Cluster/",sep=""))
	dir.create(paste(resultdir,"/Heatmaps",sep=""))
	dir.create(paste(resultdir,"/FP",sep=""))
	dir.create(paste(resultdir,"/bed",sep=""))
	dir.create(paste(resultdir,"/MCR",sep=""))
	dir.create(paste(resultdir,"/threshold/",sep=""))# when automating
	dir.create(paste(resultdir,"/output/",sep=""))# when automating
	qcdir = paste(resultdir,"/QC",sep="")
		dir.create(qcdir)
	print (qcdir)
	target_path = "../../../html/onlinetool/temp"
target = paste(target_path,project,"targets.txt",sep="/")
file.copy(target,qcdir)
file.copy(target,resultdir)
annotation_file = paste(target_path,project,"annotations.out",sep="/")
file.copy(annotation_file,resultdir)
object_out_file = paste(target_path,project,"obj.out",sep="/")
file.copy(object_out_file,resultdir)
chip_file = paste(target_path,project,"chip.txt",sep="/")
file.copy(chip_file,resultdir)
threshold_file = paste(target_path,project,"threshold.txt",sep="/")
threshold_dir = paste(resultdir,"/threshold",sep="")
file.copy(threshold_file,threshold_dir)


if (analysisType == "paired") {
	if (analysisMethod == "ASCAT") {
		#####set path for the text file containing normal samples
		normal_target_path = "../../../html/onlinetool/temp"
		normal_target = paste(normal_target_path,project,"normal.txt",sep="/")
		file.copy(normal_target,qcdir)
	}
}

source("mysetupdir.R")
mysetupdir(platform,locations,targets,type.analysisType,folder)
######### Functions ################

#This function is to prepare the input to ASCAT  
#This function collates the log2ratio - LRR files generated from Aroma for each of the samples in the dataset and also collates all the BAF files
# for each of the samples - my_collect_ascat is for the processing of tumor samples

my_collect_ascat <-function (dir, colname) 
{
    setwd(dir)
    file_list <- list.files(".",paste(colname,"*",sep=""))
    first.data <- read.table(file_list[1], header = TRUE, sep = "\t")
    if (length(file_list) > 1) {
        for (i in 2:length(file_list)) {
            data <- read.table(file_list[i], header = TRUE, sep = "\t")
            newcol = sub(".txt", "", file_list[i])
            if (colname == "LRR") {
                latestcol = sub("LRR_","",newcol)
            }
            if (colname == "BAF") {
                latestcol = sub("BAF_","",newcol)
            }
            first.data[, latestcol] = data[rownames(first.data), 
                colname]
                length <- colnames(first.data)
                if (length[3] == "LRR") {
                	new_name <- file_list[1]
                	latest <- sub("LRR_","",new_name)
                	length[3] <- sub(".txt","",latest)
                	print(length[3])
                	colnames(first.data) <- length
                }
                if (length[3] == "BAF") {
                	new_name <- file_list[1]
                	latest <- sub("BAF_","",new_name)
                	length[3] <- sub(".txt","",latest)
                	print(length[3])
                	colnames(first.data) <- length
                }
                
                
                print ("I am here NOW")
        }
    }


    Name <- seq(0, nrow(first.data)-1, by = 1)
    new_first_data <- cbind(Name, first.data)
    write.table(first.data, paste(colname, "R_input.txt", 
        sep = ""), sep = "\t", quote = FALSE, row.names = TRUE)
        print ("I GOT NAME HERE")
}
#This function is to prepare the input to ASCAT
# This function collates the log2ratio - LRR files generated from Aroma for each of the samples in the dataset and also collates all the BAF files
# for each of the samples - my_collect_ascat is for the processing of normal samples

my_collect_ascat_normal <-function (dir, colname) 
{
    setwd(dir)
    file_list <- list.files(".",paste(colname,"*",sep=""))
    first.data <- read.table(file_list[1], header = TRUE, sep = "\t")
    if (length(file_list) > 1) {
        for (i in 2:length(file_list)) {
            data <- read.table(file_list[i], header = TRUE, sep = "\t")
            newcol = sub(".txt", "", file_list[i])
            if (colname == "LRR") {
                latestcol = sub("LRR_","",newcol)
            }
             if (colname == "BAF") {
                latestcol = sub("BAF_","",newcol)
            }
            first.data[, latestcol] = data[rownames(first.data), 
                colname]
                length <- colnames(first.data)
                if (length[3] == "LRR") {
                	new_name <- file_list[1]
                	latest <- sub("LRR_","",new_name)
                	length[3] <- sub(".txt","",latest)
                	print(length[3])
                	colnames(first.data) <- length
                }
                if (length[3] == "BAF") {
                	new_name <- file_list[1]
                	latest <- sub("BAF_","",new_name)
                	length[3] <- sub(".txt","",latest)
                	print(length[3])
                	colnames(first.data) <- length
                }
                print ("I am here NOW")
        }
    }

   
    Name <- seq(0, nrow(first.data)-1, by = 1)
    new_first_data <- cbind(Name, first.data)
    write.table(first.data, paste("normal_",colname, "R_input.txt", 
        sep = ""), sep = "\t", quote = FALSE, row.names = TRUE)
        print ("I GOT NAME HERE")
}





##########################################################END OF  function ##########################################################



#####################################Code below writes locations.txt out for the dataset projects analysed #######################################



targs<-read.table(target,sep="\t",as.is=TRUE)
 targets = targs[1, ]
 

    pd <- read.table(targets, sep = "\t", as.is = TRUE, strip.white = TRUE, 
        header = T)

fact = pd
targets=paste(qcdir,"/targets.txt",sep="")
if(type == "CEL") {
	hapmap=0
# Parameters to allow for 1 or 2 arrays
	if(is.na(cdf_array2) || (cdf_array2 == "NA")) {
		cdf_arrays=c(cdf_array1)
	} else {
		cdf_arrays=c(cdf_array1,cdf_array2)
	}
	
	if(analysisType=="paired") {
		if (analysisMethod == "CBS") {
					
		write.table(as.matrix(c(Cancer),nrow=1,ncol=1),file=paste(qcdir,"/locations",sep=""),row.names=F,col.names=F,quote=F)
		locations=paste(qcdir,"/locations",sep="")
		}
		
		if (analysisMethod == "ASCAT") {
		write.table(as.matrix(c(Cancer,Normaldir),nrow=2,ncol=1),file=paste(qcdir,"/locations",sep=""),row.names=F,col.names=F,quote=F)
			locations=paste(qcdir,"/locations",sep="")	
		}
		
		}
		
		
		
		targets=paste(qcdir,"/targets.txt",sep="")
		
	}
        if (analysisMethod == "CBS") {
	if (analysisType == "unpaired") {
		
			write.table(as.matrix(c(Cancer,Normaldir),nrow=2,ncol=1),file=paste(qcdir,"/locations",sep=""),row.names=F,col.names=F,quote=F)
			locations=paste(qcdir,"/locations",sep="")
			}
				if(baseline =="hapmap") {
					hapmap=1
				}
                                }
	
	patients=pd[,1]
	rownames(fact)=pd[,1]




#####if the thershold for calling gains/losses has been left empty by the user - the default is set to 0.2
if(is.na(threshold)) {
	cgh_threshold=0.2
}


###### End of From previous version, set up files   ##############


#library("ominer")
#setwd("/var/onlinetool")           ### Fix setwd
#setupDir(platform,locations,targets,type,analysisTypeType,folder,Normal,paste("/var/onlinetool/Results/",Cancer,sep=""))   ###folder == directory where uploads are loaded to from interface, Normal name of normals directory, Cancer name of project
mysetupDir()
#####to upload a dataset from GEO

######for this bit - it is okay for an unpaired analysis with ASCAT BUT for a paired analysis this needs to work differently as all four files need to be in one txt file i.e. 1 LRR file each for paired and unpaired
#######and also 1 BAF file each one for paired and one for unpaired
LLR_output=paste(resultdir,Cancer,sep="")






if (type == "CEL") {
	if (analysisMethod == "ASCAT") {    ######run Aroma to prepare the files for ASCAT for Cancer samples 
#	source("/var/onlinetool/predictGG.R")
#	source("/var/onlinetool/ascat.R")
#	source("/var/onlinetool/aspcf.R")
	setwd(aromadir)
	project <- Cancer
	dir.create(paste("output/",project,sep=""))
        outf_dir <- paste("output/",project,sep="")
	#dataSet=paste(project,",ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2",sep="") ###will need the "ra" for this for affySNP6
	if (platform == "six") {
	dataSet=paste(project,",ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2",sep="")
	}
	else {
	dataSet=paste(project,",ACC,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2",sep="")
}
	locations = paste("../../../www/cgi-bin/onlinetool/version_2/ominer_results/",Cancer,"/QC/","locations",sep="")
		targets = paste("../../../www/cgi-bin/onlinetool/version_2/",targets,sep="")
		ascat_dir <- getwd()
		print ("I am here")
		print (ascat_dir)
source("../../../www/cgi-bin/onlinetool/version_2/run_ASCAT_input.R")

print (platform)
	run_ASCAT_input (locations,targets,platform,dataSet) 
	# Move files from the top level Aroma directory to an output dir
	output_files = list.files(".","*.txt")
        
    for (outputfile in output_files)
    {
       	file.copy(paste(aromadir,outputfile,sep="/"),paste(aromadir,"output",project,outputfile,sep="/"))
       	file.remove(paste(aromadir,outputfile,sep="/"))  ####Unblock this!!!!!!
    }
    
    if (platform == "100" || platform == "500") {
       two_chip_dir <- outf_dir
    
    
    target <- paste ("../../../www/html/onlinetool/temp/",project,"/","targets.txt",sep="")
    targs<-read.table(target,sep="\t",as.is=TRUE)
 targets = targs[1, ]
 

    pd <- read.table(targets, sep = "\t", as.is = TRUE, strip.white = TRUE, 
        header = T)
    setwd(two_chip_dir)
    file.names=list.files(getwd())
    
    chip1 <- colnames(pd[1])  ##xba chip
    chip2 <- colnames(pd[2])  ####hind chip
    filenames <- pd[3]
    
    patients <- length(pd$Name) 
        
    for (i in 1:length(pd$Name)) {

        BAF_file_xba  <- paste("BAF_",pd$Name[i], "_",chip2,".txt",sep="") 
        BAF_file_hind <- paste("BAF_",pd$Name[i], "_",chip1,".txt",sep="")
        bfx <- read.table(BAF_file_xba)
        bfh <- read.table(BAF_file_hind)
        BAF_file <- rbind(bfx,bfh)
        write.table(BAF_file,file=paste("BAF_",pd$Name[i],".txt",sep=""),sep="\t",row.names=FALSE)
        
        
        LRR_file_xba  <- paste("LRR_",pd$Name[i], "_",chip2,".txt",sep="") 
        LRR_file_hind <- paste("LRR_",pd$Name[i], "_",chip1,".txt",sep="")
        lfx <- read.table(LRR_file_xba)
        lfh <- read.table(LRR_file_hind)
        LRR_file <- rbind(lfx,lfh)
        
    write.table(LRR_file,file=paste("LRR_",pd$Name[i],".txt",sep=""),sep="\t",row.names=FALSE)
         
    }
    
    system ("rm BAF_P*_Mapping*.txt")  ####remove intermediary files so there is no confusion when creating BAF and LRR files 
    system ("rm LRR_P*_Mapping*.txt")
    setwd("../../")
    }
#	#my_collect_ascat(paste("output/",project,sep=""),"LRR",1781656) 
###dir is the name of the directory containing the LRR and BAF files that are to be put together in one matrix, colname == BAF or LRR np is the number of probes - 1 in your array
	
        my_collect_ascat(paste("output/",project,sep=""),"LRR") 
        ###dir is the name of the directory containing the LRR and BAF files that are to be put together in one matrix, colname == BAF or LRR np is the number of probes - 1 in your array
	setwd(aromadir)
	
	my_collect_ascat(paste("output/",project,sep=""),"BAF") ###dir is the name of teh directory conatining teh LRR and BAF files that are to be put together in one matrix, colname == BAF or LRR np is the number of probes - 1 in your array
	setwd(aromadir)
    setwd(paste("output",project,sep="/"))

## Create input.txt - lists the two files that are needed
    file1=paste(aromadir,"output/",project,"/LRRR_input.txt",sep="")
    file2=paste(aromadir,"output/",project,"/BAFR_input.txt",sep="")

file1="LRRR_input.txt"
file2="BAFR_input.txt"

if (analysisType == "unpaired") {
	write.table(as.matrix(c(file1,file2),nrow=2,ncol=1),file=paste(folder,"/input_ascat.txt",sep=""),row.names=F,col.names=F,quote=F)
    input_file=paste(folder,"/input_ascat.txt",sep="")
    print(analysisType)
    print(input_file)
    }
    ### Platform - need to standardise these to always be what is required for aroma, or add another variable
input_file = paste("../../../../../www/html/onlinetool/temp/",Cancer,"/input_ascat.txt",sep="")


if (analysisType == "paired") {
	setwd(aromadir)
	if (analysisMethod == "ASCAT"){
		dir.create(paste("output/",Normaldir,sep=""))
		if (platform == "six") {
	dataSet=paste(Normaldir,",ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2",sep="")
	}
	else{
		dataSet=paste(Normaldir,",ACC,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2",sep="")
		}
		
locations = Normaldir
	locations = paste("../../../www/cgi-bin/onlinetool/version_2/ominer_results/",Cancer,"/QC/","locations",sep="")
	
		targets = paste("../../../www/cgi-bin/onlinetool/version_2/ominer_results/",Cancer,"/QC","/normal.txt",sep="")
		ascat_dir <- getwd()
		print ("I am here")
		print (ascat_dir)
source("../../../www/cgi-bin/onlinetool/version_2/run_ASCAT_input_normals.R")
	run_ASCAT_input (locations,targets,platform,dataSet)
	 ###locations and targets are the same as for previous run_aroma function dataSet == projectName and LRR_output 

# Move files from the top level Aroma directory to an output dir #####unhash from here to the next bracket line 932
	output_files = list.files(".","*.txt")
    for (outputfile in output_files)
   {
       	file.copy(paste(aromadir,outputfile,sep="/"),paste(aromadir,"output",Normaldir,outputfile,sep="/"))
      	file.remove(paste(aromadir,outputfile,sep="/"))
    }

	my_collect_ascat_normal(paste("output/",Normaldir,sep=""),"LRR") ###dir is the name of teh directory conatining teh LRR and BAF files that are to be put together in one matrix, colname == BAF or LRR np is the number of probes - 1 in your array
	setwd(aromadir)
	
	my_collect_ascat_normal(paste("output/",Normaldir,sep=""),"BAF") ###dir is the name of teh directory conatining teh LRR and BAF files that are to be put together in one matrix, colname == BAF or LRR np is the number of probes - 1 in your array
	setwd(aromadir)
    setwd(paste("output",Normaldir,sep="/"))

## Create input.txt - lists the two files that are needed
    file1=paste(aromadir,"output/",project,"/LRRR_input.txt",sep="")
    file2=paste(aromadir,"output/",project,"/BAFR_input.txt",sep="")

file3="normal_LRRR_input.txt"
file4="normal_BAFR_input.txt"
	write.table(as.matrix(c(file1,file2,file3,file4),nrow=2,ncol=1),file=paste(folder,"/input_ascat.txt",sep=""),row.names=F,col.names=F,quote=F)
    input_file=paste(folder,"/input_ascat.txt",sep="")
		
	}
        
        
}

	if (analysisType == "paired") {   ####unhash from here to line 949
        
if(analysisMethod == "ASCAT") {
	write.table(as.matrix(c(file1,file2,file3,file4),nrow=2,ncol=1),file=paste(folder,"/input_ascat.txt",sep=""),row.names=F,col.names=F,quote=F)
   file.copy("normal_LRRR_input.txt",paste("../",project,sep=""))
    file.copy("normal_BAFR_input.txt",paste("../",project,sep=""))
    print(analysisType)
    print(input_file)
    
    ### Platform - need to standardise these to always be what is required for aroma, or add another variable
input_file = paste("../../../../../www/html/onlinetool/temp/",Cancer,"/input_ascat.txt",sep="")
change_back_dir <- paste("../",project,sep="")
setwd(change_back_dir)
	}

	}
        
        
        
        
        
	
	
		
	my_run_paired_unpaired_ascat(input_file,analysisType,platform,type) ###input = name of a txt file containing the names of teh BAF/LRR matrice(s), analysisType = paired/unpaired, platform == platform chosen
	setwd("../../../../../www/cgi-bin/onlinetool/version_2/")
        output_dir_remove <- paste("../../../../onlinetool/version_2/Aroma/output","/",project,"/",sep="")
        system("rm -rf output_dir_remove")
        
	source("ASCAT_output_parser.R")
		
	
	ASCAT_segment_parser(analysisType,project)
	source("ASCAT_segments_plot.R")
	print ("I found this file")
	mydirg = getwd()
	print ("this is my directory")
	print (mydirg)
	ASCAT_segments_plot(analysisType,project)
	print ("I have run ASCAT_segments_plot.R")
	
	
	
}
	}
	
if (type=="log2ratio") {
	if (analysisMethod == "ASCAT") {
		if (analysisType == "unpaired") {
		LRRR_input_path = "../../../html/onlinetool/temp"
		LRRR_file = paste(LRRR_input_path,"/",project,"normalised_raw.txt",sep="/")    #####change this but interface keep saving LRRR.txt as "normalised.txt"
		file.copy(LRRR_file,qcdir)
		BAFR_input_path = "../../../html/onlinetool/temp"
		BAFR_file = paste(BAFR_input_path,"/",project,"BAFR_input.txt",sep="/")
		file.copy(BAFR_file,qcdir)
		file1="normalised_raw.txt"
		file2="BAFR_input.txt"
                file1_path = paste("ominer_results/",project,"/QC/normalised_raw.txt",sep="")
                file2_path = paste("ominer_results/",project,"/QC/BAFR_input.txt",sep="")
                
	write.table(as.matrix(c(file1_path,file2_path),nrow=2,ncol=1),file=paste(qcdir,"/input_ascat.txt",sep=""),row.names=F,col.names=F,quote=F)
                input_file <- paste(qcdir,"/input_ascat.txt",sep="")
		my_run_paired_unpaired_ascat(input_file,analysisType,platform,type) ######will need to change this os that it looks for input files in the qcdir 
		#######then do output parsre freq_plots etc
		
		
		######format of input is a text file containing the names of the BAF and LRRR matrices for tumour and normal samples.
		######format of file input for an unpaired analysisType is name of LRR tumor matrix file 
											#####				name of BAF tumor matrix file
											
		#######format of input file for a paired analysisType is 
												###name of tumor LogR matrix file
												###name of tumor BAf matrix file
												###name of normal LRR matrix file
												###name of normal BAF matrix file									
		
		
		
		
		
		
	
	}
	if (analysisType == "paired") {
				LRRR_input_path = "../../../html/onlinetool/temp"
		LRRR_file = paste(LRRR_input_path,"/",project,"LRRR_input.txt",sep="/")
		file.copy(LRRR_file,qcdir)
		BAFR_input_path = "../../../html/onlinetool/temp"
		BAFR_file = paste(BAFR_input_path,"/",project,"BAFR_input.txt",sep="/")
		file.copy(BAFR_file,qcdir)
		normal_LRRR_file = paste(LRRR_input_path,"/",project,"normal_LRRR_input.txt",sep="/")
		file.copy(normal_LRRR_file,qcdir)
		normal_BAFR_file <- paste(BAFR_input_path,"/",project,"normal_BAFR_input.txt",sep="/")
		file.copy(normal_BAFR_file,qcdir)
		file1="LRRR_input.txt"
		file2="BAFR_input.txt"
		file3="LRRR_normal_input.txt"
		file4="BAFR_normal_input.txt"
	write.table(as.matrix(c(file1,file2,file3,file4),nrow=2,ncol=1),file=paste(qcdir,"/input_ascat.txt",sep=""),row.names=F,col.names=F,quote=F)
        input_file <- "qcdir/input_ascat.txt"
         if(is.na(threshold)) {
 		
                }
                else {
                    print ("Using user defined thresholds")
                    threshold=as.numeric(threshold)
                    gain= threshold
                    loss=threshold*-1
                    }
	my_run_paired_unpaired_ascat(input_file,analysisType,platform) #####will need to change this so it looks for input files in qcdir
####then do output parser etc and freq plots etc



		}
}
}


if (analysisMethod == "paired_CBS") {
	if (type == "ratios") {
		
	}
}

if (analysisMethod == "CBS") {
		if (type == "CEL") {
	setwd(aromadir)

		dfg <- getwd()
		print(dfg)
		source("../../../www/cgi-bin/onlinetool/version_2/runAroma_logratios.R")

		locations = paste("../../../www/cgi-bin/onlinetool/version_2/",locations,sep="")
		targets = paste("../../../www/cgi-bin/onlinetool/version_2/ominer_results/",Cancer,"/targets.txt",sep="")
		print (targets)
		qcdir = paste("../../../www/cgi-bin/onlinetool/version_2/",qcdir,sep="")
		inputdir <- paste(cghwebdir,"/input/",sep="")
		outputdir = paste("../../../www/cgi-bin/onlinetool/version_2/",inputdir,sep="")
		#cghwebdir = paste(resultdir,"/cghweb/",sep="")
		run_aroma(platform,analysisType,targets,locations,outputdir,hapmap,qcdir)
		setwd("../../../www/cgi-bin/onlinetool/version_2/")
if (platform == "100") {
	if (platform == "500") {
		
				source("merge_aroma.R")
		cghwebinput = outputdir
				
	}
	
}
}
}

current_dir <- getwd()
print ("This is my current directory")
print (current_dir)
if (analysisMethod == "CBS") {
if (type == "CEL"|| type == "log2ratio") {
	#if (type == "log2ratio") {
		#locations = paste("../../../www/cgi-bin/onlinetool/version_2/ominer_results/",Cancer,locations,sep="")
		targets = paste("ominer_results/",Cancer,"/QC/targets.txt",sep="")
		#outputdir = paste("../../../www/cgi-bin/onlinetool/version_2/",cghwebdir,"/input/",sep="")
outputdir = paste(cghwebdir,"/input/",sep="")

#cghwebinput = outputdir
		this_my_dir <- getwd()
		print ("this is my CBS location")
		print (this_my_dir)
	print ("I am trying to run CBS")
#source("../../../www/cgi-bin/onlinetool/version_2/run_CBS.R")
source("run_CBS.R")
 	runCGH(outputdir,targets,Cancer,type,platform)
 	setwd("../../../../")
 	
 	 	
 	}
 	
 	
 	if (type == "CEL"||type == "log2ratio"||type == "segmented") {
 		#setwd("../../../www/cgi-bin/onlinetool/version_2/")
                
                if(is.na(threshold)) {
 		source("threshold_calc_function.R")
 		thresholdcalcfunction(Cancer)
 		print ("i calculated the thresholds")
                }
                else {
                    print ("Using user defined thresholds")
                    threshold=as.numeric(threshold)
                    gain= threshold
                    loss=threshold*-1
                    amp=1
                    deletion=-1
		    LOCV <- getwd()
		    print ("THIS IS MY THRESHOLDS LOCATION")
		    print (LOCV)
                     write.table(data.frame(gain, loss, amp, deletion), paste("www/cgi-bin/onlinetool/version_2/ominer_results/",Cancer,"/cghweb/Matrix/threshold.txt",sep=""), 
        sep = "\t", quote = FALSE, row.names = FALSE)
	LOCATION_1 <- getwd()
	print ("THIS IS WHERE I AM COPYING THE THRESHOLDS FILE FROM")
	print (LOCATION_1)
                      file.copy("threshold.txt",paste("www/cgi-bin/onlinetool/version_2/ominer_results/",Cancer,"cghweb/Matrix/threshold.txt",sep=""))
                }
 	setwd("www/cgi-bin/onlinetool/version_2")
	LOCATION <- getwd()
	print ("THIS IS WHERE I AM NOW")
	print (LOCATION)
 	source("get_bin_coded.R")
 	getBinaryCoded(Cancer)	
 	print ("I got binary coded")
 	setwd("../../../")
 	source("filter_res.R")
	filterBinary(snp_number,Cancer)
	setwd("../../../")
 	}
 	
 	if (type == "CEL"||type == "log2ratio"||type == "segmented"||type =="smoothed") {
 		source("get_binary_filtered.R")
 	filterBinary(snp_number)
 	print ("I filtered binary coded")
 	setwd("../../../")
 	source("getregions_mod.R")
 	getRegions(snp_number,type)
 	print (" I got regions")
 	setwd("../../../")
 	source("make_frequency_plots.R")
 	target = paste("ominer_results/",Cancer,"/targets.txt",sep="")
 	results = paste("ominer_results/",Cancer,"/output/results.txt",sep="")
 	makeFreqPlot(results,target,Cancer,type,name="frequency_plot")
 	#print ("I generated frequency_plots")
 	#####generate frequency plots that are filtered 
 	results = paste("ominer_results/",Cancer,"/output/binary_coded_filtered.txt",sep="")
 	makeFreqPlot(results,target,Cancer,type,name="FrequencyPlotsFiltered")
 	source("freqplot_cgh.R")
 	####Need to add a name argument to this function so can generate both filtered & unfiltered
 	freq_plot_cgh(target,results,type, Cancer)

 	####currently theses are unfiltered freq plots 
 	#####Neeed to do one for filtered data for "group view"
 	####e.g res.unfiltered <- results = paste("ominer_results/",Cancer,"/output/binary_coded.txt",sep="")
 		makeFreqPlot(results,target,Cancer,type,name="FrequencyPlots")
 		#####also do unfiltered for each csome
 		freq_plot_cgh(target,results,type, Cancer)
 		
 		results = paste("ominer_results/",Cancer,"/output/binary_coded_filtered.txt",sep="")
makeFreqPlot(results,target,Cancer,type,name="FrequencyPlotsFiltered")
 		
 		source("generate_log2ratio_plots_mod.R")
	#####Make log2ratio plots using binary coded unfiltered text file as input
	#if (type)
	
	if(type=="CEL" || type=="log2ratio" || type=="segmented") {
	data <- paste(resultdir,"/output/R_input.txt",sep="")
	#log2ratioPlots(data)
	log2ratioPlots(data,target,name="log2ratioPlot_",type)	
	print ("I generated log2ratio plots")
	######make filtered log2ratio plots using the binary coded filtered file as input
	data <- paste(resultdir,"/output/R_input_filtered.txt",sep="")
	log2ratioPlots(data,target,name="log2ratio_filtered",type)
	}
	else {
		target = paste("ominer_results/",Cancer,"/QC/targets.txt",sep="")
		data <- paste(resultdir,"/output/results.txt",sep="")
		#print ("This is what I am inputting")
		#print (data)
		log2ratioPlots(data,target,name="log2ratioPlot_",type)	
	data <- paste(resultdir,"/output/binary_coded_filtered.txt",sep="")
log2ratioPlots(data,target,name="log2ratio_filtered",type)
	}
	
	
	
	#setwd("../../../")
	results = paste("ominer_results/",Cancer,"/output/binary_coded_filtered.txt",sep="")
	source("generate_CNA_cluster.R")
	getCluster(Cancer,target,results)
	mydir <- getwd()
	print ("I am here NOW")
	print(mydir)	
		print ("I generated a cluster plot")
		
		#system ("./R --no-save --args project < ../../generate_plots.R")
	fd <- "../output/R_input.txt"
	target <- "../QC/targets.txt"
	if (mcr == "msa") {
	source("runMSA_method.R")
	runMSA(fd,target,type,platform,Cancer)
	}
	if (mcr == "CGHRegions") {
	source("runCGHRegions.R")
	runCGHRegions(type,Cancer,platform)
	}
	
	}
}
#}

if (analysisType == "CBS") {
	if (type == "ratios") {
		
	}
}

FPdir <- paste(resultdir,"/FP",sep="")
setwd(FPdir)
FPzip <- "zip -r FLFrequencyPlots.zip *FrequencyPlotschr*.png"
system(FPzip)
setwd("../../../")
#command="/home/cgh/cgh-bin/mail.pl"
#url="Your results are now available at : http://o-miner.org/cgi-bin/onlinetool/version_2/cghOutputView.pl?"


resultdir = "/var/www/cgi-bin/onlinetool/version_2/ominer_results/"
setwd(resultdir)
zipcommand=paste("zip -r ",Cancer,".zip ",Cancer,sep="")
print(zipcommand)
system(zipcommand)
file.rename(paste(resultdir,Cancer,".zip",sep=""),paste(resultdir,Cancer,"/output/",Cancer,".zip",sep=""))

command="/var/www/cgi-bin/onlinetool/version_2/mail.pl"
url="Your results are now available at : http://o-miner.org/cgi-bin/onlinetool/version_2/cghOutputView.pl?"
subject <- "O-miner Results"
mailcommand=paste(command," ",email," \"",subject,"\" \"",url,project,"\"",sep="")



system(mailcommand)
