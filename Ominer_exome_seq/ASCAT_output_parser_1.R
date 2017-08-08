############################################################################################################
#
#File: ASCAT_output_parser_1.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to manipulate the output of ASCAT and get it ready for annotation and plotting in frequency plots
#####Arguments - (1) main sample name - sample i.e if filenames are CSCC_1 then sample is equal to "CSCC", (2) project - this is the unique name given to each project.
######this will fail if rerunning the same project because *_segments.txt already exists there
ASCAT_segment_parser <- function(sample,project) {

segment_file_name <- paste("segments_",sample,".txt",sep="")


filename <- (paste("ominer_results/",project,"/Regions/",segment_file_name,sep=""))



seg_info <- read.table(filename, sep= "\t",header=T)
gaindata <- subset(seg_info,CNEventTypePloidyCorrected == "Gain",)

outputdir <- paste("ominer_results/",project,"/Regions/",sep="")
proc_gain_data <- cbind(gaindata[1],gaindata[6],gaindata[2],gaindata[3],gaindata[4],gaindata[11],gaindata[13],gaindata[14])
colnames(proc_gain_data) <- c("patient","nbre","Chromosomes","pos_start","pos_end","correctedtcnploidy","CNA","NOTE")
write.table(proc_gain_data,paste(outputdir,"gains",sample,".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)


lossdata <- subset(seg_info,CNEventTypePloidyCorrected == "Loss",)
proc_loss_data <- cbind(lossdata[1],lossdata[6],lossdata[2],lossdata[3],lossdata[4],lossdata[11],lossdata[13],lossdata[14])
colnames(proc_loss_data) <- c("patient","nbre","Chromosomes","pos_start","pos_end","correctedtcnploidy","CNA","NOTE")
write.table(proc_loss_data,paste(outputdir,"losses",sample,".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)



CN_LOH <- subset(seg_info,CNEventTypePloidyCorrected == "CN_LOH",)
proc_CLOH_data <- cbind(CN_LOH[1],CN_LOH[6],CN_LOH[2],CN_LOH[3],CN_LOH[4],CN_LOH[11],CN_LOH[14],CN_LOH[13])
colnames(proc_CLOH_data) <- c("patient","nbre","Chromosomes","pos_start","pos_end","correctedtcnploidy","CNA","NOTE")
write.table(proc_CLOH_data,paste(outputdir,"LOH",sample,".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)

	
	
}