##########################################################################################
#
#File name: aqm.R
#
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#
#
#
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
##########################################################################################
##########################################################################################
#Description: normalisation of raw data from Affymetrix CEL files
#Arguments are:
#target_file = path to the target file 
#dat = affymetrix object created from CEL files
#analysis = paired/unpaired analysis
#project is the given project name
##########################################################################################
run_aqm <- function(target_file,dat,analysis,project) {
library("arrayQualityMetrics")

library("affy")
	library("simpleaffy")
	library("arrayMvout")
	library("affyPLM")
pd <- read.AnnotatedDataFrame(target_file, header = T, row.name = "Name", 
            sep = "\t")
#dat <- ReadAffy(filenames = paste(pd$FileName, sep = ""), 
 #           sampleNames = sampleNames(pd), phenoData = (pd))
  #              cat("Running arrayQualityMetrics",file="error.out",sep="\n",append=TRUE)
  outaqm = paste("ominer_results/",project,"/QC","/AQM",sep="")
                aqm=try(arrayQualityMetrics(expressionset = dat, outdir = outaqm, force = TRUE, do.logtransform = FALSE))
        
# QC summary data frame
qc.Adata <- qc(dat)
        
        pdtable <- read.table(target_file,sep="\t",as.is=T,header=T,row.names="Name")
     qcsummary<-data.frame(pdtable,avbg(qc.Adata),percent.present(qc.Adata),ratios(qc.Adata))
       write.table(qcsummary, paste("ominer_results/",project,"/QC","/qcsummary.txt",sep=""), sep = "\t",quote=F,row.names=FALSE)
       
}