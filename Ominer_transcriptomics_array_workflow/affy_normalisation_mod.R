##########################################################################################
#
#File name: affy_normalisation_mod.R
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
#Arguments are:dataType this refers to the type of data input i.e CEL for raw data files and normalised for the normalised expression matrix
#normalisation - this refers to the normalisation type either (1) rma, (2) gcrma or (3) trma
#target_file is the full path to the target file
#dat this is the affy object created from the raw CEL files
#project is the given project name
##########################################################################################
Normalisation <- function (normalisation, target,dat,project) 
{
   
        Norm <- NULL
    if (normalisation == "rma") {
        Norm <- rma(dat)
    }
    if (normalisation == "gcrma") {
library("gcrma")       
 Norm <- gcrma(dat)
    }
    if (normalisation == "trma") {
        library("trma")
        Norm <- trma(dat)
    }
#filename <- paste(normalisation,".exp",sep="")
 #   write.exprs(Norm, paste("ominer_results/",project,"/",norm, "/",filename, sep = ""), 
  #      sep = "\t")
   #}
   write.exprs(Norm,file=paste("ominer_results/",project,"/","norm","/","normalised",".txt",sep=""),sep="\t")
   }