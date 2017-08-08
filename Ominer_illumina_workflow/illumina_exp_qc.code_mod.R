##########################################################################################
#
#File name: illumina_exp_qc.code_mod.R
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
#dataset this is the matrix of unnormalised data
#project is the given project name
#target_file is the full path to the target file
#analysis - run either a paired or an unpaired analysis i.e. paired unpaired
#normalisation - the normalisation method chosen by the user either: vsn,rsn,ssn,quantile or loess
##########################################################################################
illumina_qc <- function(dataset,project,target_file,analysis,normalisation) {
library(lumi)
#target = "target_GSE29650.txt"
#pd <- read.table(target,header=T,row.name="Name",sep="\t")
###read in data from GEO GSE32571
#dataset = "GSE32571_non_normalized_mod.txt"
	pdtable<-read.table(target_file,sep="\t",header=T)
	library("affyPLM")
pd <- read.AnnotatedDataFrame(target_file, header = T, sep = "\t")

data = lumiR( dataset, sep = NULL, detectionTh = 0.01, na.rm = TRUE, lib = NULL)

###Quality control

QCsummary <- data@QC$sampleSummary
write.table(QCsummary, paste("ominer_results/",project,"/QC","/QC_summary.txt",sep=""),sep="\t", quote=FALSE,row.names=FALSE)
 #paste("ominer_results/",project,"/QC","/qc.txt",sep=""), sep = "\t")
#####Try out the prepare2write function for formatting data for output

prepare2write <- function (x) {
        
        x <- cbind(rownames(x), x)
        return(x)
}


#write.table(prepare2write(QCsummary),  paste("ominer_results/",project,"/QC","/QC_summary.txt",sep=""),sep="\t", row.names=c("Sample","mean","standard deviation","detection_rate","distance_to_sample_mean"))

##### signal intensity before normalisation
pdf( paste("ominer_results/",project,"/QC","/density.pdf",sep=""))
plot(data, what='density')
dev.off()


##### Sample relations using hierarchical clustering before normalisation
pdf(paste("ominer_results/",project,"/QC","/Sample_relations_cluster.pdf",sep=""), pointsize = 8 ,width = 0.2*length(sampleNames(data)), height = 6)
plot(data, what='sampleRelation', hang=-1)
dev.off()

##### Sample relations using multidimensional scaling before normalisation
pdf(paste("ominer_results/",project,"/QC","/Sample_relations_MDS.pdf",sep=""))
plot(data, what='sampleRelation', method='mds')
dev.off()

##### multivariate outlier detection
library(arrayMvout)
#library(affyPLM)
ii = ArrayOutliers(data, alpha = 0.001, pc2use=1:3)
if(is.null(ii[1]$outl))  {
	revised = "0"
	
	
	}
else {
#cnames <- c("Sample","mean","standard deviation","detection_rate","distance_to_sample_mean")
write.table(prepare2write(ii[[1]]), paste("ominer_results/",project,"/QC","/outliers.txt",sep=""),sep="\t", row.names=FALSE,col.names=TRUE,quote=F)
#write.table(ii[1],paste("ominer_results/",project,"/QC","/outliers.txt",sep=""),sep="\t", row.names=FALSE,quote=F)
#i#i[[1]]
revised = "1"

# write.table(ii[[1]],  paste("ominer_results/",project,"/QC","/outliers.txt",sep=""),sep = "\t",col.names= TRUE)

out <- read.table(paste("ominer_results/",project,"/QC","/outliers.txt",sep=""), sep = "\t", as.is = T, 
        header = T)
   
   
    	    	 colnames(out) = c("Sample", "mean","standard deviation","detection_rate","distance_to_sample_mean")
    	  f <- data.frame(pdtable[out$Sample, ]$Name, pdtable[out$Sample, 
            ]$Target, out)
        colnames(f) = c("pdtable", "Sample", "mean", 
            "standard_deviation", "detection_rate", "distance_to_sample_mean")
        write.table(f, paste("ominer_results/",project,"/QC","/outliers.fn.txt",sep=""), sep = "\t", quote = F, 
            row.names = FALSE)

    	}
    		
    		
    		
    		    	
        #f <- data.frame(pdtable[out$Sample, ]$Name, pdtable[out$Sample, 
         #   ]$Target, out)
        #colnames(f) = c("pdtable", "Sample", "mean", 
         #   "standard_deviation", "detection_rate", "distance_to_sample_mean")
        #write.table(f, paste("ominer_results/",project,"/QC","/outliers.fn.txt",sep=""), sep = "\t", quote = F, 
         #   row.names = FALSE)
    write.table(prepare2write(ii[[2]]),  paste("ominer_results/",project,"/QC","/summary.txt",sep=""),sep="\t", row.names=FALSE,quote= F)
      good <- read.table(paste("ominer_results/",project,"/QC","/summary.txt",sep=""), sep = "\t", header = T, 
        as.is = T, row.names = NULL)
    colnames(good) = c("Sample", "mean", "standard_deviation","detection_rate", "distance_to_sample_mean")
    
    #f <- data.frame(pdtable[good$Sample, ]$Name, pdtable[good$Sample,]$Target, good)
    #colnames(f) = c("Target", "Sample", "mean","standard_deviation","detection_rate","distance_to_sample_mean")
    write.table(good,  paste("ominer_results/",project,"/QC","/good.fn.txt",sep=""), sep = "\t", quote = F, row.names = F)
    
    
   if(analysis=="paired") {
                        pairstoremove=pairs=pd[as.vector(ii[[1]]$sample)]$Pairs
                        unique(pairstoremove)
                        samp=NULL
                        for(pair in 1:length(pairstoremove)) 
                        samp=c(samp,sampleNames(pd[pd$Pairs==pairs[pair]]))
                        pd=pd[setdiff(sampleNames(pd),samp)]
                        sel<-setdiff(rownames(pdtable),samp)            
                        
                } else {
                	
                	if (revised == 1) {
                        #pd=pd[setdiff(sampleNames(pd),as.vector(ii[[1]]$Sample))] # Remove outliers from pdtable file
                         pd=pd[setdiff(sampleNames(pd),out$Sample)] # Remove outliers from pdtable file
                        #sel<-setdiff(rownames(pdtable),as.vector(ii[[1]]$sample))
                        sel<-setdiff(pdtable$Name,out$Sample)
                                       
                                       
                                       #####here pdtable already ahas the bad sample removed
 newpdtable <- read.table(target_file,sep="\t",as.is=T,header=T,row.name ="Name")
               
    newpdtable<-newpdtable[sel,]
    new_names <- rownames(newpdtable)
    new_target <- newpdtable$Target
    
    
       allnew_pdtable <- cbind(new_names,new_target)
colnames(allnew_pdtable) <- c("Name","Target")
 write.table(allnew_pdtable,paste("ominer_results/",project,"/QC","/target_qc.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
 }
}

       if (revised ==1 ) {      #####need to write out a new targets file with bad quality one being removed AND also need to removed it from the normalised/unnormalised matrix
       	  tr_dataset <- read.table(dataset, sep = "\t", header = T, as.is = T)
       	  #####keep good$Sample for the data matrix 
       	  keep = NULL
       	 #good <- read.table(paste("ominer_results/",project,"/QC","/summary.txt",sep=""), sep = "\t", header = T, 
        #as.is = T)

       	  x <- length(good$Sample)
       	  #for (i in length(out$Sample)) {
       	  for (i in seq(from=1,to=x,by=1)) {
       	  	print (i)
       	  	wanted <- paste("\\",good$Sample[i],"\\b",sep="")
       	  	
       	  	keep <- c(keep,grep(wanted,colnames(tr_dataset)))
       	  	print (keep)
       	  	}
       	  	
       	  	new_dat <- subset(tr_dataset[keep])
       	  	TargetID <- tr_dataset$TargetID
new_matrix <- cbind(TargetID,new_dat)
       	  #	print (i)
       	  	#detection_names <- c(detection_names,grep(t_file$Name[i],colnames(n_matrix)))
       	  	#new_dat <- tr_dataset[, -grep(out$Sample[i],colnames(tr_dataset))]
       	  #	new_dat <- c(new_dat[,-grep(out$Sample[i],colnames(tr_dataset))])       	  
       	   #     	  }
       	        	  write.table(new_matrix,paste("ominer_results/",project,"/QC","/new_dataset",sep=""), sep="\t", quote=FALSE, row.names=FALSE)

rm(data)
new_data_path = paste("ominer_results",project,"/QC","new_dataset",sep="/")
       	  	data = lumiR(new_data_path, sep = NULL, detectionTh = 0.01, na.rm = TRUE, lib = NULL)

       
        
        }
        #restart with outliers removed



##### Report detected outliers in the project workspace
#samples2exclude = data.frame(c(studyID, dataDir, paste(rownames(ii[[1]]), collapse = ',')))

#write.table( t(samples2exclude),paste("ominer_results/",project,"/QC",/outliers.txt",sep="") , row.names = FALSE, col.names = c("DatasetName", "DataDir", "Samples2exclude"),  quote = FALSE, sep="\t" )

pdf(paste("ominer_results/",project,"/QC","/outliers.pdf",sep=""), pointsize = 6)
plot(ii, choices = c(1, 3))
dev.off()




replace0s <- function (lumiBatch) {
        
        for (i in 1:ncol(exprs(lumiBatch))) {
                
                exprs.ordered <- sort(exprs(lumiBatch)[,i])
                exprs.ordered <- exprs.ordered[ exprs.ordered != 0 ]
                                
                exprs(lumiBatch)[exprs(lumiBatch)[,i]==0] <- exprs.ordered[1]
        }       
        return(lumiBatch)
}
##### Replace signal intensities of '0', which introduce errors in RSN, with minimal value detected on the array
if (normalisation == "rsn") {
data  <- replace0s(data)
data.N <- lumiN(data, method='rsn')
}
if (normalisation == "vsn") {
data  <- replace0s(data)
data.N <- lumiN(data, method='vsn')
}
if (normalisation == "ssn") {
data  <- replace0s(data)
data.N <- lumiN(data, method='ssn')
}
if (normalisation == "quantile") {
data  <- replace0s(data)
data.N <- lumiN(data, method='quantile')
}
if (normalisation == "loess") {
data  <- replace0s(data)
data.N <- lumiN(data, method='quantile')
}

#####Quality control normalised data
data.N.Q <- lumiQ(data.N)

QCsummary <- data.N.Q@QC$sampleSummary

##### Write QC summary
write.table(prepare2write(QCsummary), paste("ominer_results/",project,"/QC","/QC_summary.rsn.txt",sep=""), sep="\t", row.names=FALSE)

##### signal intensity after RSN normalisation
pdf(paste("ominer_results/",project,"density.rsn.pdf",sep=""))
plot(data.N.Q, what='density')
dev.off()

##### Sample relations using hierarchical clustering after RSN normalisation
#pdf(paste("ominer_results/",project,"/QC","Sample_relations_cluster.rsn.pdf",sep=""),pointsize = 8 ,width = 0.2*length(sampleNames(data)), height = 6)
#plot(data.N.Q, what='sampleRelation', hang=-1)
#dev.off()

##### Sample relations using multidimensional scaling after RSN normalisation
pdf(paste("ominer_results/",project,"/QC/","Sample_relations_MDS.rsn.pdf",sep=""))
plot(data.N.Q, what='sampleRelation', method='mds')
dev.off()

##### Box-and-whisker plot
#pdf(paste("ominer_results/"),project,"/QC","boxplot.rsn.pdf", sep=""),pointsize = 8 ,width = 0.2*length(sampleNames(data)), height = 6)
#par(mar=c(13, 4, 3, 2))
#boxplot(data,col="red", las = 2) # Generates boxplot of non-normalized intensity values.
#boxplot(data.frame(log2(exprs(data.N.Q))),col="blue", main="Normalised RSN data", las = 2) # Generates boxplot of normalized log intensity values.
#dev.off()

##### Write expression file
write.table(prepare2write(log2(exprs(data.N.Q))), paste("ominer_results/",project,"/norm","/normalised.txt",sep=""),sep="\t", quote=FALSE,row.names=FALSE)
 }
 
 
 
 
 
 
 
 
 
 
 
 
 