#########################
#
#File runAroma_logratios.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to calculate the log2ratios from raw data files
#(1) Platform - name of Affymetrix array platform used, (2) analysisType - paired/unpaired,(3) targets - this is a .txt file with the full path of the location of the target file for 
###tumor samples and another full path for the .target.txt file for normal samples,(4) locations ia .txt file containing the full path to normal samples and another for the full path
####to tumor samples(5) full pathof the directory to write output of aroma to (5) hapmap - 1 if Hapmap is to be used and 0 if user chooses the baseline, (6) is the fullpath to the QC ####directory.2
#output is a matrix of log2ratios
#################################################################################################################

run_aroma <- function(platform,analysisType,targets,locations,outputdir,hapmap,qcdir)  { 
	cghwebinput=outputdir
	pathdir=outputdir	
	library("aroma.affymetrix")
	library("sfit")	

	verbose <- Arguments$getVerbose(-10, timestamp=TRUE);
	setOption(aromaSettings, "memory/ram", 200.0);
	locs<-read.table(locations,sep="\t",as.is=T)		
	Cancer=locs[1,] ###location of tumour files
	Normals=locs[2,] ###location of normal samples #here i am assuming that locs[2,] could be hapmaptype"JPT"?
	#if (analysisType =="paired"){Normals = Cancer;}	

	targs <- read.table(targets,sep="\t",as.is=TRUE)
	targets=targs[1,] 	
	pd <- read.table(targets, sep="\t", as.is=TRUE,strip.white=TRUE,header=T);

	if (platform=="100") { 
		cdf_arrays=c("Mapping50K_Hind240","Mapping50K_Xba240")
		plm="RmaCnPlm";
		if(analysisType=="paired") {
			chip_pairs=targs[2,] #2 files for 100k and 500k # asuming p1 is hind and p2 is xba are pairs files same order as cdf_arrays
		}
	}
	if (platform=="50xba") { 
		cdf_arrays="Mapping50K_Xba240"
		plm="RmaCnPlm";
		if(analysisType=="paired") {
			chip_pairs=targs[2,] #2 files for 100k and 500k # asuming p1 is hind and p2 is xba are pairs files same order as cdf_arrays
		}
	}
	if (platform=="50hind") { 
		cdf_arrays="Mapping50K_Hind240"
		plm="RmaCnPlm";
		if(analysisType=="paired") {
			chip_pairs=targs[2,] #2 files for 100k and 500k # asuming p1 is hind and p2 is xba are pairs files same order as cdf_arrays
		}
	}
	
	if (platform=="500") { 
		cdf_arrays=c("Mapping250K_Sty","Mapping250K_Nsp")
		plm="RmaCnPlm"	
		if(analysisType=="paired") {
			chip_pairs=targs[2,] #2 files for 100k and 500k # asuming p1 and p2 are pairs files same order as cdf_arrays
		}
	}

	if (platform=="250sty") { 
		cdf_arrays="Mapping250K_Sty"
		plm="RmaCnPlm"			
		if(analysisType=="paired") {
			chip_pairs=targs[2,] #2 files for 100k and 500k # asuming p1 and p2 are pairs files same order as cdf_arrays
		}
	}

	if (platform=="250nsp") { 
		cdf_arrays="Mapping250K_Nsp"
		plm="RmaCnPlm"			
		if(analysisType=="paired") {
			chip_pairs=targs[2,] #2 files for 100k and 500k # asuming p1 and p2 are pairs files same order as cdf_arrays
		}
	}

	
	if (platform=="six") { 
		cdf_arrays="GenomeWideSNP_6,Full"
		plm="AvgCnPlm"
		if(analysisType=="paired") {
			chip_pairs=targs[2,] #1 file
		}
	}
	
	if (platform=="five") { 
		cdf_arrays="GenomeWideSNP_5,Full,r2"
		plm="AvgCnPlm"
		if(analysisType=="paired") {
			chip_pairs=targs[2,] #1 file
		}
	}
	
	if (platform=="131") { 
		cdf_arrays="Mapping10K_Xba131"	
		plm="RmaCnPlm"
		if(analysisType=="paired") {
			chip_pairs=targs[2,] #1 file
		}
	}
	
	if (platform=="142") { 
		cdf_arrays="Mapping10K_Xba142"	
		plm="RmaCnPlm"
		if(analysisType=="paired") {
			chip_pairs=targs[2,] #1 file
		}
	}					
						
#get the names of the hapmap right for paired arrays	
# Fullnames translator
	fnt <- function(names, ...) {
		names <- gsub("_(HIND|XBA|STY|NSP)", "", names);
		names;	
	} 
	
	
# Preprocessing using CRMAv2 
	dsC=list();
	dsN=list();
	dsN_average=list();
	for (j in 1:length(cdf_arrays))
	{
		chipType=cdf_arrays[j];		
		cdf <- AffymetrixCdfFile$byChipType(chipType); 
		csR <- AffymetrixCelSet$byName(Cancer,  cdf=cdf)
cs <- csR
output_qc_plot <- paste("../../../www/cgi-bin/onlinetool/version_2/ominer_results/",Cancer,"/QC/",chipType,".pdf",sep="")
print ("this is the location of my QC density plots")
print (output_qc_plot)
pdf (output_qc_plot)
par(mar=c(4,4,1,1)+0.1)
plotDensity(cs,lwd=2,ylim=c(0,0.40))
stext(side=3,pos=0,getFullName(cs))
dev.off()


		dsC[[chipType]] = doCRMAv2(Cancer, cdf=cdf, combineAlleles = TRUE, plm=plm,verbose=verbose); 		

		if (hapmap=="1" && analysisType=="unpaired") 
			{
				print ("I am trying to find my hapmap here and am failing")
				
				dsN[[chipType]] <-  AromaUnitTotalCnBinarySet$byName(Normals,tags="ACC,-XY,BPN,-XY,RMA,A+B,FLN,-XY", chipType=chipType);
				print ("this is the file that I am looking for NOW/HAPMAP")
				print(dsN[[chipType]])
				#dsN[[chipType]] = doCRMAv2(Normals, cdf=cdf, combineAlleles = TRUE, plm=plm,verbose=verbose); 
				setFullNamesTranslator(dsN[[chipType]], fnt);
				print ("I am failing here")
			}
		#if (hapmap=="0" && analysisType =="paired")	
			#{
				#dsN[[chipType]] = doCRMAv2(Normals, cdf=cdf, combineAlleles = TRUE, plm=plm,verbose=verbose); 
			#}
	if (analysisType == "unpaired" && hapmap == "0") {
		   print ("THE ERROR IS HERE I THINK THAT I AM AN UNPAIRED ANALYSIS")
		dsN[[chipType]] = doCRMAv2(Normals, cdf=cdf, combineAlleles = TRUE, plm=plm,verbose=verbose); 
             
		
		
	}
	
	
	
		if(analysisType =="unpaired")
		{
			dsN_average[[chipType]] = getAverageFile(dsN[[chipType]])
		}			
		if(analysisType =="paired")
		{
#make sure it is in the right matching pair # Splitting the dataset in normal-tumor pairs

			infoSamples <- read.table(chip_pairs, sep="\t", header=TRUE, as.is=TRUE,strip.white=TRUE);
#Rows in this file are in the same order as arrays in the aroma dataset.
			idxC <-  sub(".CEL","",infoSamples[,1]);
			idxN <-  sub(".CEL","",infoSamples[,2]);
#			stopifnot(all(idxC %in% getNames(dsC[[chipType]])));
#			stopifnot(all(idxN %in% getNames(dsN[[chipType]])));							
			dsN[[chipType]] <- extract(dsC[[chipType]], idxN);
			dsC[[chipType]] <- extract(dsC[[chipType]], idxC);
#Circular binary segmentation T vs N: Paired model (first N sample pairs with first T sample)					  	
		}
	}
	print("cbsmodel")
	if(analysisType =="paired")		{	cbs <- CbsModel(dsC,dsN);}
	if(analysisType =="unpaired")	{	cbs <- CbsModel(dsC,dsN_average);}		
	
#both cdfs in case 100 or 500k
	cdfs <- lapply(cdf_arrays, FUN=function(chipType) {AffymetrixCdfFile$byChipType(chipType)})
	gis <- lapply(cdfs, getGenomeInformation)
	dataf <- lapply(gis, FUN=function(chipType) {readDataFrame(chipType)})
	sis <- lapply(cdfs, getUnitNames)

	data=do.call("rbind", dataf)
	data$unitName=unlist(sis)	
	
	arrays = getNames(cbs)

	
# CAUTION - parallel libraries may interfere with Aroma Exit function.  	
	#library(doMC)
	#registerDoMC(15)
	#foreach (i = 1:length(arrays))  %dopar%	
	
	here <- getwd()
	print ("I am in this directory NOW\n")
		print (here)
		
		
for (i in 1:length(arrays) ) #for each patient - do this in parallel
	{
			for (j in 1:22 ) #for each chromosome 
		{
				rawCNs <- extractRawCopyNumbers(cbs, array=i, chromosome=j)#both cdfs
				rawCNs <- as.data.frame(rawCNs)	#both cdfs			
				all = subset(data,chromosome==j)
				test=merge(all, rawCNs, by.x="physicalPosition", by.y="x", all =T)
				#write.table(test,paste(cghwebinput,arrays[i],".test_txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE,append=TRUE,col.names=TRUE)
				
				new_matrix= cbind(test$physicalPosition,test$chromosome.x,test$unitName,test$cn)
				
				colnames(new_matrix) = c("Position","Chromosome","ProbeID","LogRatio")
				#write.table(all,paste(cghwebinput,arrays[i],".all_txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE,append=TRUE,col.names=TRUE)
				#write.table(rawCNs,paste(cghwebinput,arrays[i],".raw_CNs_txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE,append=TRUE,col.names=TRUE)

				write.table(new_matrix[,c("ProbeID","Chromosome","Position","LogRatio")],paste(cghwebinput,arrays[i],".txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE,append=TRUE,col.names=FALSE)
			}}
				
#prepare CGHinput
	patients=pd$Name ##### make sure all cel files have been renamed to pd$Name.
	date()
	for (i in 1:length(patients))  #- parallel write out file
	#foreach (i = 1:length(patients))  %dopar%
		{
			HX.data = read.table(paste(cghwebinput,"/",patients[i],".txt",sep=""), sep="\t", header=FALSE,as.is=TRUE)		
			colnames(HX.data) = c("ProbeID","Chromosome","Position","LogRatio")
			HX.data <- subset(HX.data, !duplicated(HX.data[,"ProbeID"]))
			HX.data <- subset(HX.data, !duplicated(paste(HX.data[,"Chromosome"],HX.data[,"Position"])))
			HX.data$Position = as.numeric(HX.data$Position)
			HX.data$Chromosome = as.numeric(HX.data$Chromosome)
			HX.data$LogRatio = as.numeric(HX.data$LogRatio)
			HX.data$ProbeID = sub("SNP_A-","",HX.data$ProbeID)
			HX.data$ProbeID = sub("CN_","",HX.data$ProbeID)
			HX.data$ProbeID = sub("AFFX-SNP_","",HX.data$ProbeID)
			write.table(HX.data,paste(pathdir,patients[i],".txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
	}
		
}
