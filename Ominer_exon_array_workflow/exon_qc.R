##########################################################################################
#
#File name: exon_qc.R
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
#Description: Run initial assessment of exon array CEL file data using aroma.affymetrix
#Arguments are:platform this refers to the type of Affymetrix exon array platform that is being analysed i.e. Genechip Human Exon 1.0ST
#project this is the given project name
#aromadir - this is the full path to the aroma.affymetrix directory 
#mydir - this is the full path to the directory where the master R script is called from 
##########################################################################################


###chipType = "HuEx-1_0-st-v2"
###dataset_name = "pancreas"
#annotation.data.exon <- "annotationData/chipTypes/HuEx-1_0-st-v2"
     
     
     
qc <- function(platform,project,aromadir,mydir){
setwd(aromadir)
library(aroma.affymetrix)
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
###Get annotation files
tag = "coreR2,A20070914,EP"
annot_level = "coreR2"
chipType <- platform
cdf <- AffymetrixCdfFile$byChipType(chipType,tags =tag)
print(cdf)



cs <- AffymetrixCelSet$byName(project,cdf=cdf)
print(cs) 

setCdf(cs,cdf) 
##Background adjustment and normalization
bc <- RmaBackgroundCorrection(cs,tag=annot_level) 
print(bc)
#RAM: 0.00MB
###Do background correction
csBC <- process(bc,verbose=verbose)
print(bc)
csBC <- process(bc,verbose=verbose)
print(bc)


####set up quantile normalization method
qn <- QuantileNormalization(csBC,typesToUpdate="pm")
print(qn)

csN <- process(qn, verbose=verbose)
####Summarization
getCdf(csN)


######Fit summary model transcript expression
#######fit a summary of the entire transcript (estimate overall expression for the transcript)

plmTr <- ExonRmaPlm(csN,mergeGroups=TRUE,tag="coreProbesetsGeneExpression")
print(plmTr)



####run the model
fit(plmTr,verbose=verbose)  
####generate data frame 
cesTr <- getChipEffectSet(plmTr)
trFit <- extractDataFrame(cesTr,units=NULL,addNames=TRUE)
setwd(mydir)
write.table(trFit, paste("ominer_results/",project,"/","transcript/norm","/","normalised.txt", 
            sep = ""), sep = "\t", col.names=TRUE,row.names = TRUE, quote = FALSE)
            
#plmTrtest <- ExonRmaPlm(csN,mergeGroups=TRUE)
#print(plmTrtest)



####run the model
#fit(plmTrtest,verbose=verbose)  
####generate data frame 
#cesTrtest <- getChipEffectSet(plmTrtest)
#trFittest <- extractDataFrame(cesTrtest,units=NULL,addNames=TRUE)
setwd(mydir)
#write.table(trFittest, paste("ominer_results/",project,"/","norm","/","trfittest.txt", 
            #sep = ""), sep = "\t", col.names=TRUE,row.names = TRUE, quote = FALSE)            
            
            
            
            
####Exon expression fit exon by exon
plmEx <- ExonRmaPlm(csN,mergeGroups=FALSE,tag="coreProbesetsExonExpression")
print(plmEx)

fit(plmEx, verbose=verbose)  ###got to here
###generate a data frame
cesEx <- getChipEffectSet(plmEx)
ExFit <- extractDataFrame(cesEx, units=NULL,addNames=TRUE)
####Quality assessment of PLM fit
write.table(ExFit, paste("ominer_results/",project,"/","exon","/norm/","normalised.txt", 
            sep = ""), sep = "\t", col.names=TRUE,row.names = TRUE, quote = FALSE)
            
            

            
####do some alternative splicing analysis by calculating FIRMA scores:
firma <- FirmaModel(plmTr)
fit(firma, verbose=verbose)
#fs <- getChipEffectSet(firma)
fs <- getFirmaScores(firma)
#FmFit <- log2(extractDataFrame(fs))
FmFit <- extractDataFrame(fs)

firma.result  <- merge(FmFit, ExFit[, c(1,2,5)])
####Quality assessment of PLM fit
write.table(firma.result, paste("ominer_results/",project,"/","splicing","/norm/","normalised.txt", sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)



####NUSE and RLE plots
qamTr <- QualityAssessmentModel(plmTr)
 pdf(paste("ominer_results/",project,"/transcript/QC","/NUSE.pdf",sep=""))

par(mar=c(4,4,1,1)+0.1)
plotNuse(qamTr)
dev.off()
 pdf(paste("ominer_results/",project,"/transcript/QC","/RLE.pdf",sep=""))
#pdf("RLE.pdf")
par(mar=c(4,4,1,1)+0.1)
plotRle(qamTr)
dev.off()
rs <- calculateResidualSet(plmTr, verbose=verbose)
#rsdata <- extractDataFrame(rs)
#write.table(rsdata, paste("ominer_results/",project,"/QC","rsdata.txt", sep=""),
 #       sep = "\t", row.names = TRUE,col.names=TRUE, quote = FALSE)
}
