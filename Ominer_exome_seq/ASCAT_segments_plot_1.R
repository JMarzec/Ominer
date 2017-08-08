#########script to parse the output of ASCAT in orser to format a text file for input to the copynumber package

ASCAT_segments_plot <- function(sample,project) {
library(biomaRt)
#up <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
segment_file_name <- paste("segments_",sample,".txt",sep="")

filename <- (paste("ominer_results/",project,"/Regions/",segment_file_name,sep=""))

fseg  <- read.table(filename, sep= "\t",header=T)
print ("THIS IS THE FILENAME")
print (fseg)
outputdir <- paste("ominer_results/",project,"/FP/",sep="")


fseg.res <- vector()
#new_line <- NULL
#cat_file <- NULL
#new_file <- NULL
theFilters = c("chromosome_name","start","end")   ######or illumina_humanht_v3
theAttributes = c("chromosome_name","band","start_position","end_position")

new_test1 <- vector()

test1 <- getBM(attributes=theAttributes,filters=theFilters,values=list(fseg$chr,fseg$start,fseg$end),mart=ensembl)
print ("THIS IS TEST1")
print (test1)

farm <- substr(test1$band, 1, 1)
print ("THIS IS FARM")
print (farm)
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
 
outputdir <- paste("ominer_results/",project,"/Regions/",sep="")

write.table(tmp.res, paste(outputdir,sample,"_processed_segments.txt",sep=""))


}







