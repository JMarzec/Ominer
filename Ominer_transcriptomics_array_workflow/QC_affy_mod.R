##########################################################################################
#
#File name: QC_affy_mod.R
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
#Description: Running quality control analysis using the R pacakge arrayMvout
#Arguments are:dataType this refers to the type of data input i.e CEL for raw data files and normalised for the normalised expression matrix
#dat this is the affy object created from the raw CEL files
#target_file is the full path to the target file
#project is the given project name
##########################################################################################
runQC <-function (dataType, dat,target_file,project) 
{
	
	
	pd <- read.AnnotatedDataFrame(target_file, header = T, row.name = "Name", 
            sep = "\t")
	pdtable<-read.table(target_file,sep="\t",as.is=T,header=T,row.names="Name")

     qc.Adata <- qc(dat)
   
   
    pdf(paste("ominer_results/",project,"/QC","/qc.pdf",sep=""))
    plot(qc.Adata)
    plot(qc.Adata, usemid = T)
    dev.off()
    write.table(cbind(avbg(qc.Adata), sfs(qc.Adata), percent.present(qc.Adata)), 
        paste("ominer_results/",project,"/QC","/qc.txt",sep=""), sep = "\t")
    library(affyPLM)
    Aset <- fitPLM(dat)
    pdf(paste("ominer_results/",project,"/QC","/PLM.pdf",sep=""))
    RLE(Aset, main = "RLE")
    NUSE(Aset, main = "NUSE")
    boxplot(Aset)
    Mbox(Aset)
    dev.off()
    write.table(NUSE(Aset, type = "stats"), paste("ominer_results/",project,"/QC","/NUSE.txt",sep=""), sep = "\t")
    write.table(RLE(Aset, type = "stats"), paste("ominer_results/",project,"/QC","/RLE.txt",sep=""), sep = "\t")
    deg <- AffyRNAdeg(dat)
    write.table(summaryAffyRNAdeg(deg), paste("ominer_results/",project,"/QC","/deg.txt",sep=""), sep = "\t")
    pdf(paste("ominer_results/",project,"/QC","/deg.pdf",sep=""))
    plotAffyRNAdeg(deg)
    dev.off()
    ratios(qc.Adata)
    avbg(qc.Adata)
    maxbg(qc.Adata)
    minbg(qc.Adata)
    spikeInProbes(qc.Adata)
    qcProbes(qc.Adata)
    percent.present(qc.Adata)
    plot(qc.Adata)
    sfs(qc.Adata)
    #target_file(qc.Adata)
    library(arrayMvout)
    ii = ArrayOutliers(dat, alpha = 0.01, qcOut = qc.Adata, plmOut = Aset, 
        degOut = deg)
    write.table(ii[[1]],  paste("ominer_results/",project,"/QC","/outliers.txt",sep=""), sep = "\t")
    write.table(ii[[2]],  paste("ominer_results/",project,"/QC","/summary.txt",sep=""), sep = "\t", quote = F)
    pdf(paste("ominer_results/",project,"/QC","/outliers.dat.pdf",sep=""))
    plot(ii, choices = c(1, 3))
    dev.off()
    out <- read.table(paste("ominer_results/",project,"/QC","/outliers.txt",sep=""), sep = "\t", as.is = T, 
        header = T)
    colnames(out) = c("Sample", "averageBG", "ScaleFactor", "Present", 
        "HSAC07", "GAPDH", "NUSE", "RLE", "RLE_IQR", "RNAslope")
    if (nrow(out) > 0) {
    	revised = "1"
    	}
    		else {
    			revised = "0"
    		}
    	
    	
        f <- data.frame(pdtable[out$Sample, ]$FileName, pdtable[out$Sample, 
            ]$Target, out)
        colnames(f) = c("FileName", "pdtable", "Sample", "averageBG", 
            "ScaleFactor", "Present", "HSAC07", "GAPDH", "NUSE", "RLE", "RLE_IQR", "RNAslope")
        write.table(f, paste("ominer_results/",project,"/QC","/outliers.fn.txt",sep=""), sep = "\t", quote = F, 
            row.names = FALSE)
    
      good <- read.table(paste("ominer_results/",project,"/QC","/summary.txt",sep=""), sep = "\t", header = T, 
        as.is = T, row.names = NULL)
    colnames(good) = c("Sample", "averageBG", "ScaleFactor","Present", "HSAC07", "GAPDH", "NUSE", "RLE", "RLE_IQR", "RNAslope")
    f <- data.frame(pdtable[good$Sample, ]$FileName, pdtable[good$Sample,]$Target, good)
    colnames(f) = c("FileName", "Target", "Sample", "averageBG", "ScaleFactor", "Present", "HSAC07", "GAPDH", "NUSE", "RLE", "RLE_IQR", "RNAslope")
    write.table(f,  paste("ominer_results/",project,"/QC","/good.fn.txt",sep=""), sep = "\t", quote = F, row.names = F)
   
   
   ####Remove arrays not passing QC
   if(analysis=="paired") {
                        pairstoremove=pairs=pd[as.vector(ii[[1]]$samp)]$Pairs
                        unique(pairstoremove)
                        samp=NULL
                        for(pair in 1:length(pairstoremove)) 
                        samp=c(samp,sampleNames(pd[pd$Pairs==pairs[pair]]))
                        pd=pd[setdiff(sampleNames(pd),samp)]
                        sel<-setdiff(rownames(pdtable),samp)            
                        
                } else {
                        pd=pd[setdiff(sampleNames(pd),as.vector(ii[[1]]$samp))] # Remove outliers from pdtable file
                        sel<-setdiff(rownames(pdtable),as.vector(ii[[1]]$samp))
                                       }
                
    pdtable<-pdtable[sel,]
       pdtable$Name<-rownames(pdtable)
       filenames <- pdtable$FileName
            names_file <- pdtable$Name
              Targets_new <- pdtable$Target
allnew_pdtable <- cbind(filenames,names_file,Targets_new)
colnames(allnew_pdtable) <- c("FileName","Name","Target")

       if (revised ==1 ) {       
        write.table(allnew_pdtable,paste("ominer_results/",project,"/QC","/target_qc.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
        }
        #restart with outliers removed
      

}
