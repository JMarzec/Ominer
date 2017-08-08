############################################################################################################
#
#File name mysetupdir.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: R script to set up required directory structure for analysis using the O-miner pipeline. This is required due to the rigid
####structure demanded by aroma.affymetrix
#######Set up directory structure for either the ASCAT method or the CBS method based on which the user has chosen
#######This is very important as the rest of the pipeline relies on finding the files in these directories - also aroma which is used for the #######processing of raw files from both ASCAT and CBS relies on files being found in a specific directory structure.
#######Also, copies over user supplied data for any stage of analysis from the CBS pipeline i.e. log2ratios,smoothed and segmented
#######Arguments: (1)platform (2)locations(3)targets(4)type(5)analysisType,(6) folder - this is the folder that raw/processed data is uploaded to from the interface
mysetupDir <-function () {



	
# Set up directory structure for pipeline based on input type
	if(type == "CEL") {
        if(analysisMethod == "ASCAT") {
			#if (analysisType == "paired") {
				print ("I am here")
			dir.create(paste(aromadir,"rawData/",Normaldir,sep=""))
			print("created dir for normals")
			print(Normaldir)
			print ("I created the directory")
				
			}
		dir.create(paste(aromadir,"rawData/",Cancer,sep="")) 
		print("check normal dir")
               
		if(analysisMethod=="CBS" && analysisType == "unpaired" && baseline== "user") {
			dir.create(paste(aromadir,"rawData/",Normaldir,sep=""))
			print("created dir for normals")
			print(Normaldir)
		}
	    
		
		print(paste("Reading CEL files from :",folder))
		
#read the chip names and patient names from the target file & copy over to correct Aroma rawData directory
		for (ar in 1:length(cdf_arrays)) 
		{
			cdf_array = cdf_arrays[ar];
			dir.create(paste(aromadir,"rawData/",Cancer,"/",cdf_array,sep=""))
                       
			if(analysisMethod=="CBS" && analysisType == "unpaired" && baseline=="user") {
			  dir.create(paste(aromadir,"rawData/",Normaldir,"/",cdf_array,sep=""))
			}
                        
			if(analysisType == "paired") {
				if (analysisMethod == "ASCAT") {
			  dir.create(paste(aromadir,"rawData/",Normaldir,"/",cdf_array,sep=""))
			  
			}
			}
			
			
			print("created directory for array")
			print(cdf_array)
			chip=pd[,ar]
			chip= sub(".CEL","",chip)
			chip=paste(chip,".CEL",sep="");
			name=pd$Name
			print ("this is my name")
			print(chip)
#Copy file over from user directory to aroma directory
			for (i in 1:length(patients))
			{
				print(paste(folder,"/",chip[i],sep=""))
				print(paste(aromadir,"rawData/",Cancer,"/",colnames(pd[ar]),"/",name[i],".CEL",sep=""))
				
				
				origin_file <- paste(folder,"/",chip[i],sep="")
				print ("This is my orginal path")
				print (origin_file)
				print 
				
				
				
				file.copy(paste(folder,"/",chip[i],sep=""), paste(aromadir,"rawData/",Cancer,"/",colnames(pd[ar]),"/",name[i],".CEL",sep=""))
				print(chip[i])
				
				
				
			}
		}
		
#get the normal baseline copied over as well
		
			
			if(analysisType == "paired") {
				copyDir=Cancer
			} else {
				copyDir=Normaldir
			}
	
			
			
			if (analysisMethod == "ASCAT") {   ####if ASCAT paired need to have two separate directories one for cancer samples and another for normals
				if(analysisType == "paired") {
				#copyDir=Cancer
                                copyDir=Normaldir
			}
                        else {
                            copyDir=Cancer   #####id unpaired all go in one directory
                        }
			}
			if(analysisMethod == "CBS" && baseline=="user"||analysisMethod== "ASCAT" && analysisType == "paired")
		{
#specified by user
                        
			norm<-read.table(Normal,header=T,sep="\t",as.is=T)
			print(copyDir)
			if(platform=="500" || platform=="100") {
			chip1=norm[,1]
			chip2=norm[,2]
			chip2= sub(".CEL","",chip2)
			chip2=paste(chip2,".CEL",sep="");
			chips=c(chip1,chip2)
			} else {
				chip1=norm[,1]
				chips=c(chip1)
			}
			chip1= sub(".CEL","",chip1)
			chip1=paste(chip1,".CEL",sep="");
			
			
			name=norm$Name
			for (i in 1:length(chip1))
			{
				if(platform=="500" || platform=="100") {
				   file.copy(paste(folder,"/",chip2[i],sep=""), paste(aromadir,"rawData/",copyDir,"/",colnames(pd[2]),"/",name[i],".CEL",sep=""))
				}	
				file.copy(paste(folder,"/",chip1[i],sep=""), paste(aromadir,"rawData/",copyDir,"/",colnames(pd[1]),"/",name[i],".CEL",sep=""))
                               
                                
			}
		}else {

			for (ar in 1:length(cdf_arrays)) 
			{
				cdf_array = cdf_arrays[ar];
				dir.create(paste(aromadir,"probeData/",Cancer,",ACC,-XY/",sep=""))
				dir.create(paste(aromadir,"probeData/",Cancer,",ACC,-XY/",cdf_array,sep=""))
				filelist=list.files(paste(aromadir,"probeData/",Normaldir,",ACC,-XY/",cdf_array,sep=""))
				for(i in 1:length(filelist)) {

					file.copy(paste(aromadir,"probeData/",Normaldir,",ACC,-XY/",cdf_array,"/",filelist[i],sep=""),paste(aromadir,"probeData/",Cancer,",ACC,-XY/",cdf_array,"/",filelist[i],sep=""))
					}
				dir.create(paste(aromadir,"probeData/",Cancer,",ACC,-XY,BPN,-XY",sep=""))
				dir.create(paste(aromadir,"probeData/",Cancer,",ACC,-XY,BPN,-XY/",cdf_array,sep=""))
				filelist=list.files(paste(aromadir,"probeData/",Normaldir,",ACC,-XY,BPN,-XY/",cdf_array,sep=""))
				for(i in 1:length(filelist)) {

					file.copy(paste(aromadir,"probeData/",Normaldir,",ACC,-XY,BPN,-XY/",cdf_array,"/",filelist[i],sep=""),paste(aromadir,"probeData/",Cancer,",ACC,-XY,BPN,-XY/",cdf_array,"/",filelist[i],sep=""))
				}
				
				dir.create(paste(aromadir,"plmData/",Cancer,",ACC,-XY,BPN,-XY,RMA,A+B",sep=""))
				dir.create(paste(aromadir,"plmData/",Cancer,",ACC,-XY,BPN,-XY,RMA,A+B/",cdf_array,sep=""))
				filelist=list.files(paste(aromadir,"plmData/",Normaldir,",ACC,-XY,BPN,-XY,RMA,A+B/",cdf_array,sep=""))
				for(i in 1:length(filelist)) {
	
					file.copy(paste(aromadir,"plmData/",Normaldir,",ACC,-XY,BPN,-XY,RMA,A+B/",cdf_array,"/",filelist[i],sep=""),paste(aromadir,"plmData/",Cancer,",ACC,-XY,BPN,-XY,RMA,A+B/",cdf_array,"/",filelist[i],sep=""))
				}
				
				dir.create(paste(aromadir,"plmData/",Cancer,",ACC,-XY,BPN,-XY,RMA,A+B,FLN,-XY",sep=""))
				dir.create(paste(aromadir,"plmData/",Cancer,",ACC,-XY,BPN,-XY,RMA,A+B,FLN,-XY/",cdf_array,sep=""))
				filelist=list.files(paste(aromadir,"plmData/",Normaldir,",ACC,-XY,BPN,-XY,RMA,A+B,FLN,-XY/",cdf_array,sep=""))
				for(i in 1:length(filelist)) {
		
					file.copy(paste(aromadir,"plmData/",Normaldir,",ACC,-XY,BPN,-XY,RMA,A+B,FLN,-XY/",cdf_array,"/",filelist[i],sep=""),paste(aromadir,"plmData/",Cancer,",ACC,-XY,BPN,-XY,RMA,A+B,FLN,-XY/",cdf_array,"/",filelist[i],sep=""))
				}
				}	
		}		
				
				
				
		print("..Done reading CEL files from folder & copying over to aroma structure")	
		
	}
	
	dir.create(paste(cghwebdir,"/input/",sep=""))
	dir.create(paste(cghwebdir,"/output/",sep=""))
	dir.create(paste(cghwebdir,"/profiles/",sep=""))
	dir.create(paste(cghwebdir,"/Matrix/",sep=""))# when automating
	
# For log2ratios copy the relevant files over to the cghweb input directory	
	if(type=="log2ratio")  {
		allNormalisedData<-read.delim(paste(folder,"/normalised_raw.txt",sep=""),sep="\t",as.is=T,header=T,na.strings = "NA")
		for (i in 1:length(patients))
		{
					patientdata<-data.frame(allNormalisedData[,c("ProbeID","Chromosome","Position",patients[i] )])
			colnames(patientdata) = c("ProbeID","Chromosome","Position","LogRatio")
			#outputdir = paste("../../../www/cgi-bin/onlinetool/version_2/",cghwebdir,"/input/",sep="")

			write.table(patientdata,file=paste(cghwebdir,"/input/",patients[i],".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
		}		

				}
	
	if(type=="segmented") {
		file.copy(paste(folder,"/R_input.txt",sep=""),paste(cghwebdir,"/Matrix/R_input.txt",sep=""))
		file.copy(paste(folder,"/R_input.txt",sep=""),paste(resultdir,"/output/R_input.txt",sep=""))
    }
	if(type=="smoothed") {
		file.copy(paste(folder,"/binary_coded.txt",sep=""),paste(resultdir,"/output/results.txt",sep=""))
	}	
	
	unlink(paste(resultdir,Cancer,sep=""))
	dir.create(paste(resultdir,Cancer,sep=""))
	dir.create(paste(resultdir,Cancer,"/Regions/",sep=""))# when automating
	dir.create(paste(resultdir,Cancer,"/Cluster/",sep=""))
	dir.create(paste(resultdir,Cancer,"/Heatmaps",sep=""))
	dir.create(paste(resultdir,Cancer,"/FP",sep=""))
	dir.create(paste(resultdir,Cancer,"/threshold/",sep=""))# when automating
	dir.create(paste(resultdir,Cancer,"/output/",sep=""))# when automating
	
		file.copy(paste(folder,"/obj.out",sep=""),paste(resultdir,Cancer,"/obj.out",sep=""))
	file.copy(paste(folder,"/cgh.out",sep=""),paste(resultdir,Cancer,"/cgh.out",sep=""))
	file.copy(paste(folder,"/annotations.out",sep=""),paste(resultdir,Cancer,"/annotations.out",sep=""))
	
	
	if(analysisType=="paired") {
		if (analysisMethod == "CBS") {
			if (analysisMethod == "PSCBS") {
						write.table(as.matrix(c(Cancer),nrow=1,ncol=1),file=paste(qcdir,"/locations",sep=""),row.names=F,col.names=F,quote=F)
		locations=paste(qcdir,"/locations",sep="")
		}
		}
		if (analysisMethod == "ASCAT") {
		write.table(as.matrix(c(Cancer,Normaldir),nrow=2,ncol=1),file=paste(qcdir,"/locations",sep=""),row.names=F,col.names=F,quote=F)
			locations=paste(qcdir,"/locations",sep="")	
		}
		
		}
if(analysisType=="unpaired") {
	if (analysisMethod == "ASCAT") {
        write.table(as.matrix(c(Cancer,Normaldir),nrow=2,ncol=1),file=paste(qcdir,"/locations",sep=""),row.names=F,col.names=F,quote=F)
			locations=paste(qcdir,"/locations",sep="")
                        }
                        }
	
	
	
}

