############################################################################################################
#
#File: readvcf_modified.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to run an automated processed output from varscan for each sample and prepare the data for processing with ASCAT. i.e. generate the input files for ASCAT
####these are the tumor BAF values, tumor LRR values and normal BAF values.
######Arguments include (1) folder this is the file that contains the data that has been uploaded on the server, (2) sample - this is the name given to each sample e.g. if filename
######is CSCC_1 then sample will equal "CSCC", (3) outfile_name this is the name given to the output file,(4) project this is the unique name given to the project, (5) is the ####directory to which files are output.

readvcf <- function(folder,sample,outfile_name,project, outputdir) {




inputfile <- paste(folder, sample, sep="/")
tvcf <- read.table(inputfile, sep="\t")

wchr <- which( tvcf[,2]=="1" | tvcf[,2]=="2" | tvcf[,2]=="3" | tvcf[,2]=="4" | tvcf[,2]=="5" | tvcf[,2]=="6" | tvcf[,2]=="7" | tvcf[,2]=="8" | tvcf[,2]=="9" | tvcf[,2]=="10" | tvcf[,2]=="11" | tvcf[,2]=="12" | tvcf[,2]=="13" | tvcf[,2]=="14" | tvcf[,2]=="15" | tvcf[,2]=="16" | tvcf[,2]=="17" | tvcf[,2]=="18" | tvcf[,2]=="19" | tvcf[,2]=="20" | tvcf[,2]=="21" | tvcf[,2]=="22" )
tvcf <- tvcf[wchr, ]

chr.name <- tvcf[,2]
tvcf[,2] <- as.character(chr.name)
tvcf[,1] <- c(1:nrow(tvcf))
tvcf[,3] <- as.character(tvcf[,3])
tvcf[,7] <- as.character(tvcf[,7])
tvcf[,11] <- as.character(tvcf[,11])

fvcf <- as.data.frame(tvcf)

med.can <- median(fvcf$V10)
med.nor <- median(fvcf$V6)

lrr <- log10((fvcf$V10/med.can)/(fvcf$V6/med.nor))

mlrr <- cbind(name=fvcf$V1, chr=fvcf$V2, pos=fvcf$V3, test=lrr)
nbaf <- cbind(name=fvcf$V1, chr=fvcf$V2, pos=fvcf$V3, test=fvcf$V7)
tbaf <- cbind(name=fvcf$V1, chr=fvcf$V2, pos=fvcf$V3, test=fvcf$V11)
colnames(mlrr) <- c("name", "chr", "pos", outfile_name)
colnames(nbaf) <- c("name", "chr", "pos", outfile_name)
colnames(tbaf) <- c("name", "chr", "pos", outfile_name)

filelogr <- paste(outputdir, "/logr_", outfile_name, sep="")
filenbaf <- paste(outputdir, "/nbaf_", outfile_name, sep="")
filetbaf <- paste(outputdir, "/tbaf_", outfile_name, sep="")

cat (filelogr, "\n")
write.table(mlrr, filelogr,  sep="\t", quote=F, row.names=F)
write.table(nbaf, filenbaf,  sep="\t", quote=F, row.names=F)
write.table(tbaf, filetbaf,  sep="\t", quote=F, row.names=F)
}