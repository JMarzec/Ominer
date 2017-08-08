############################################################################################################
#
#File: ASCAT_segments_plot_rev.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to run an annotate the segments generated from ASCAT with chromosomenumber, end pos, start pos and band using biomaRt ready for input to copynumber package
###arguments - project name 

ASCAT_segments_plot <- function(project) {
library(biomaRt)
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
#up <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#ensembl = useDataset("hsapiens_gene_ensembl",mart=mart)
	#theAttributes = c(platform_annot,"refseq_mrna_predicted","chromosome_name","band","hgnc_symbol", "description")
        


fseg <- read.table(paste("ominer_results/",project,"/Regions/segments_final.out.txt",sep=""),row.names=NULL)	


fseg.res <- vector()
#new_line <- NULL
#cat_file <- NULL
#new_file <- NULL
theFilters = c("chromosome_name","start","end")   ######or illumina_humanht_v3
theAttributes = c("chromosome_name","band","start_position","end_position")

new_test1 <- vector()

test1 <- getBM(attributes=theAttributes,filters=theFilters,values=list(fseg$chr,fseg$start,fseg$end),mart=mart)
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







