##########################################################################################
#
#File name: illum_aqm.R
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
#Description: quality control on unnormalised data file 
#Arguments are:
#dataset = matrix of unnormalised data
#target_file = path to the target file 
#project is the given project name
##########################################################################################

illumina_aqm <- function(dataset,target_file,project)	{						
######make this as an option and arraymvout as standard for all 
library(arrayQualityMetrics)
library(lumi)
##### Customise arrayQualityMetrics reports
#target <- as.data.frame(read.table(paste(dataDir,targetFile,sep="/"),sep="\t",as.is=TRUE,header=TRUE))
#####think that File name is the normalised/unormalised expression matrix
library("affyPLM")
#pd <- read.(target_file, header = T, sep = "\t")
 pd <- as.data.frame(read.table(target_file, header = T, sep = "\t",as.is=TRUE))
rownames(pd) <- pd$Name
#fileName = pd[1,"FileName"]
#rownames(pd) <- read.table(target_file,sep="\t",as.is=TRUE,header=TRUE)[,"Name"]
colnames(pd) <- "Target"



#target <- as.data.frame(read.table(target_file,sep="\t",as.is=TRUE,header=TRUE))
#fileName = paste(dataDir, target[1,"FileName"], sep="/")
#fileName = paste(dataDir, target[1,"FileName"], sep="/")
#pd <- as.data.frame(read.table(target_file,sep="\t",as.is=TRUE,header=TRUE)[,"Target"])
#pd <- read.AnnotatedDataFrame(target_file, header = T, sep = "\t")
#target <- as.data.frame(read.table(paste(dataDir,targetFile,sep="/"),sep="\t",as.is=TRUE,header=TRUE)[,"Target"])
lumiData = lumiR( dataset, sep = "\t", detectionTh = 0.01, na.rm = TRUE, lib = NULL)
sampleNames(lumiData) <- rownames(pd)

#rownames(target) <- read.table(paste(dataDir,targetFile,sep="/"),sep="\t",as.is=TRUE,header=TRUE)[,"Name"]
#colnames(target) <- "Target"
#sampleNames(lumiData) <- rownames(target)
pd <- read.AnnotatedDataFrame(target_file, header = T, sep = "\t")
#data <- new("ExpressionSet", exprs = exprs(lumiData), phenoData = new("AnnotatedDataFrame", data=target)) 
data <- new("ExpressionSet", exprs = exprs(lumiData), phenoData = pd) 
#data <- new("ExpressionSet", exprs = exprs(lumiData), phenoData = pd) 

preparedData = prepdata(expressionset = data, intgroup = "Target", do.logtransform = TRUE)

QCboxplot <- aqm.boxplot(preparedData, subsample=20000, outlierMethod = "KS")
QCdensity <- aqm.density(preparedData)
QCheatmap <- aqm.heatmap(preparedData)
QCpca <- aqm.pca(preparedData)
QCmaplot <- aqm.maplot(preparedData, subsample=20000, Dthresh=0.15, maxNumArrays=8, nrColumns=4)
QCmeansd <- aqm.meansd(preparedData)

qm = list("Heatmap"=QCheatmap, "PCA"=QCpca, "Boxplot"=QCboxplot, "Density"=QCdensity, "MeanSD"=QCmeansd, "MAplot"=QCmaplot)
#out = ""
 outaqm = paste("ominer_results/",project,"/QC","/AQM",sep="")
aqm.writereport(modules = qm, reporttitle = "arrayQualityMetrics report", outdir = outaqm, arrayTable = pData(data))

#aqm.writereport(modules = qm, reporttitle = paste("arrayQualityMetrics report for", studyID, sep=" "), outdir = QCdir, arrayTable = pData(data))

}