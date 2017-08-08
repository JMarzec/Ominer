#######script to process the ASCAT_output and create files that can be used as input to the annotation system for the genomics_pipeline
#######first write out all inormation for each sample to a separate .txt file

######this will fail if rerunning the same project because *_segments.txt already exists there

############################################################################################################
#
#File name run_ASCAT_output_parser.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: R script to parse the output from ASCAT and get it ready into a format to generate plots
####Arguments - (1) analysisType this is either paired/unpaired to inform which .txt file should be analysed and (2) project - this is the unique name assigned to the project


ASCAT_segment_parser <- function(analysisType,project) {

if (analysisType == "paired") {
ASCAT_segments_info <- read.table(paste("ominer_results/",project,"/Regions/segments_final_paired.out.txt",sep=""))
}
if (analysisType == "unpaired")  {
ASCAT_segments_info <- read.table(paste("ominer_results/",project,"/Regions/segments_final_unpaired.out.txt",sep=""))	
}


sample_id <- ASCAT_segments_info$SampleID


sample_ID_list <- NULL
sample_ID_list <- unique(sample_id)



for (i in 1:length(sample_ID_list)) {
	file_list_ID <- NULL
		file_list_ID <- subset(ASCAT_segments_info, SampleID == sample_ID_list[i] ,)
		
		write.table(file_list_ID,paste("ominer_results/",project,"/Regions/segments/",sample_ID_list[i],sep=""),sep="\t",quote=FALSE,row.names=TRUE)
		print (i)
	}
	
#######parse each file separately

segments_dir <- paste("ominer_results/",project,"/Regions/segments",sep="")


setwd(segments_dir)

system("rm *_gains.txt")
system ("rm *_losses.txt")
system ("rm *_CN_LOH.txt")

seg_file_list <- NULL
seg_file_list = list.files()
for (i in 1:length(seg_file_list)) {
	filename <- seg_file_list[i]

	
seg_info <- read.table(filename, sep= "\t",header=T)
gaindata <- subset(seg_info,CNEventTypePloidyCorrected == "Gain",)





proc_gain_data <- cbind(gaindata[1],gaindata[6],gaindata[2],gaindata[3],gaindata[4],gaindata[11],gaindata[13],gaindata[14])
colnames(proc_gain_data) <- c("patient","nbre","Chromosomes","pos_start","pos_end","correctedtcnploidy","CNA","NOTE")
write.table(proc_gain_data,paste("../","gains",filename,".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)


lossdata <- subset(seg_info,CNEventTypePloidyCorrected == "Loss",)
proc_loss_data <- cbind(lossdata[1],lossdata[6],lossdata[2],lossdata[3],lossdata[4],lossdata[11],lossdata[13],lossdata[14])
colnames(proc_loss_data) <- c("patient","nbre","Chromosomes","pos_start","pos_end","correctedtcnploidy","CNA","NOTE")
write.table(proc_loss_data,paste("../","losses",filename,".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)



CN_LOH <- subset(seg_info,CNEventTypePloidyCorrected == "CN_LOH",)
proc_CLOH_data <- cbind(CN_LOH[1],CN_LOH[6],CN_LOH[2],CN_LOH[3],CN_LOH[4],CN_LOH[11],CN_LOH[14],CN_LOH[13])
colnames(proc_CLOH_data) <- c("patient","nbre","Chromosomes","pos_start","pos_end","correctedtcnploidy","CNA","NOTE")
write.table(proc_CLOH_data,paste("../","LOH",filename,".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)

	}
	setwd("../../../../")
}