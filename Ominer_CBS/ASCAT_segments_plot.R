############################################################################################################
#
#File name ASCAT_segments_plot_annotated.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: R script to format the ASCAT output and get it ready for input to the copynumber package
####Arguments - (1) analysisType this is either paired/unpaired to inform which .txt file should be analysed and (2) project - this is the unique name assigned to the project

ASCAT_segments_plot <- function(analysisType,project) {
library(biomaRt)
up <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

if (analysisType == "paired") {
fseg <- read.table(paste("ominer_results/",project,"/Regions/segments_final_paired.out.txt",sep=""))
}
if (analysisType == "unpaired")  {
fseg <- read.table(paste("ominer_results/",project,"/Regions/segments_final_unpaired.out.txt",sep=""))	
}

fseg.res <- vector()
#new_line <- NULL
#cat_file <- NULL
#new_file <- NULL
theFilters = c("chromosome_name","start","end")   ######or illumina_humanht_v3
theAttributes = c("chromosome_name","band","start_position","end_position")

new_test1 <- vector()

test1 <- getBM(attributes=theAttributes,filters=theFilters,values=list(fseg$chr,fseg$start,fseg$end),mart=up)
farm <- substr(test1$band, 1, 1)
number <- nrow(fseg)
wanted_band <- farm[1:number]
tmp.res <- cbind(as.character(fseg$SampleID),fseg$chr,wanted_band,fseg$start,fseg$end,as.numeric(fseg$nSNP),as.numeric(fseg$correctedTcnForPloidy))
 colnames(tmp.res) <- c("sampleID","chrom","arm","start.pos","end.pos","n.probes","mean")
new_seg_info <- as.data.frame(tmp.res)
new_seg_info$sampleID <- as.character(tmp.res[,1])
new_seg_info$chrom <- as.integer(tmp.res[,2])
new_seg_info$arm <- as.character(tmp.res[,3])
new_seg_info$start.pos <- as.integer(tmp.res[,4])
new_seg_info$end.pos <-  as.integer(tmp.res[,5])
new_seg_info$n.probes <- as.numeric(tmp.res[,6])
new_seg_info$mean <- as.numeric(tmp.res[,7])
#write.table(new_seg_info,file="ominer_results/a.sangaralingam-ASCAT_paired_250K/Regions/segments_processed_out.txt")
write.table(tmp.res, paste("ominer_results/",project,"/Regions/segments_processed_out.txt",sep=""))


}







