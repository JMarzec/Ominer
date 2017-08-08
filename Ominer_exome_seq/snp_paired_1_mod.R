############################################################################################################
#
#File: snp_paired_1_mod.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to run ASCAT for the post-processing of exome sequencing data
#Arguments: (1) sample - this is the main name given to each of the samples i.e. if the samples all have a names of CSCC_1, CSCC_2, CSCC_3 etc then 'sample' refers to "CSCC", (2) project is the name given to the results folder and
#is unique to each dataset that is analysed, (3) output dir - this is the full path to where the results folder should be written out to and (4) workdir - this is where the input files are found - input files in this case are 
#the processed output from varscan analysis
snp_paired <- function(sample, project, outputdir,workdir) {

# fthrPloidy : threshold to define gain/loss on total copy number after the ploidy correction (correctedTcnForPloidy). 
#              Gain: correctedTcnForPloidy > fthrPloidy , Loss: correctedTcnForPloidy < fthrPloidy 
# fthr       : threshold to define gain/loss on total copy number (tcn) without the ploidy correction. 
#              Gain: tcn-2 > fthr , Loss: tcn-2 < fthr 

thrPloidy <- 0.6
thr <- 0.8
fthrPloidy=as.numeric(thrPloidy)
fthr=as.numeric(thr)
outputdir <- paste("ominer_results", project, sep="/")
inputdir <- paste("ominer_results", project, sep="/")
cat(paste("ASCATcode/ascat_v2.R", sep=""))
source(paste("ASCATcode/process.ASCAT.func.1.R", sep=""));
source(paste("ASCATcode/ascat_v2.R", sep=""));

Description: Set up teh relevant directory structure
filelogr <- paste(inputdir, "/logr_", sample , sep="")
cat(filelogr)
filenbaf <- paste(inputdir, "/nbaf_", sample, sep="")
filetbaf <- paste(inputdir, "/tbaf_", sample, sep="")
filesegs <- paste(outputdir, "/Regions","/segments_", sample, ".txt",sep="")
fileploi <- paste(outputdir, "/ploidy_", sample, ".txt", sep="")
filecont <- paste(outputdir, "/content_", sample, ".txt", sep="")
ftumorfiles <- paste(outputdir, "/FP","/Tumor",sep="")
fgermlinefiles <- paste(outputdir,"/FP", "/Germline_",sep="")

fsunrisefiles <- paste(outputdir, "/FP","/sunrise", sep="")
fprofilefiles <- paste(outputdir, "/FP","/ASCATprofile",sep="")
frawprofilefiles <- paste(outputdir, "/FP","/rawprofile",sep="")
faspcffiles <- paste(outputdir, "/FP","/ASPCF", sep="")

#Description: Set up parameters to call ASCAT
ascat.bc = ascat.loadData(filelogr, filetbaf, filelogr, filenbaf);
ascat.plotRawData(ascat.bc, ftumorfiles, fgermlinefiles);
ascat.bc = ascat.aspcf(ascat.bc, workdir, outputdir);

ascat.plotSegmentedData(ascat.bc, faspcffiles);
ascat.output = ascat.runAscat(ascat.bc,  sunrisefiles=fsunrisefiles, profilefiles=fprofilefiles, rawprofilefiles=frawprofilefiles);

tumour_content = ascat.output$aberrantcellfraction;
tumour_ploidy = ascat.output$ploidy;
write.table(tumour_content, file=filecont, sep="\t", quote=F);
write.table(tumour_ploidy, file=fileploi, sep="\t", quote=F);

#ASCAT fits its output to a model with a given ploidy, so tcn must be corrected for ploidy to determine gain and loss
#Extract ploidy calls
#First, add the allele specific copy numbers to get total copy number (tcn)

ascat.segments <- organize.ascat.segments(ascat.output, ascat.bc$SNPpos);
ascat.segments.tcn <- cbind(ascat.segments, ascat.segments[,6]+ascat.segments[,7]);
colnames(ascat.segments.tcn)[8] <- "tcn";
ploidy <- ascat.output$ploidy;
names(ploidy) <- colnames(ascat.output$nA);
ascat.segments.tcn <- cbind(ascat.segments.tcn,ploidy[match(ascat.segments.tcn[,1], names(ploidy))],ascat.segments.tcn[,8] - ploidy[match(ascat.segments.tcn[,1], names(ploidy))]);
colnames(ascat.segments.tcn)[9:10] <- c("ploidy","correctedTcnForPloidy");

##extract ASCAT summary output now and write segments out to a .txt file 

segments<-ascat.segments.tcn;
size<-segments$End-segments$Start;
segments<-data.frame(SampleID=segments$SampleID,chr=segments$Chr,start=segments$Start,end=segments$End,size=size,nSNP=segments$nProbes,segments[,-(1:5)]);
CNEventTypePloidyCorrected<-rep("None",dim(segments)[1]);
CNEventType<-rep("None",dim(segments)[1]);
LOH<-rep("None",dim(segments)[1]);
LOH[which(segments$nA==0 | segments$nB==0)]<-"LOH";
CNEventTypePloidyCorrected[which(segments$correctedTcnForPloidy<=(-fthrPloidy))]<-"Loss";
CNEventTypePloidyCorrected[which(segments$correctedTcnForPloidy>=fthrPloidy)]<-"Gain";
CNEventTypePloidyCorrected[which(segments$correctedTcnForPloidy<fthrPloidy & segments$correctedTcnForPloidy>(-fthrPloidy) & (segments$nA==0 | segments$nB==0))]<-"CN_LOH";
CNEventType[which((segments$tcn-2)<=(-fthr))]<-"Loss";
CNEventType[which((segments$tcn-2)>=fthr)]<-"Gain";
CNEventType[which((segments$tcn-2)<fthr & (segments$tcn-2)>(-fthr) & (segments$nA==0 | segments$nB==0))]<-"CN_LOH";
LOH<-rep("None",dim(segments)[1]);
LOH[which(segments$nA==0 | segments$nB==0)]<-"LOH";
segments_final<-data.frame(segments,CNEventType=CNEventType,CNEventTypePloidyCorrected=CNEventTypePloidyCorrected,LOH=LOH);

write.table(segments_final,file=filesegs, sep="\t", quote=F);
}