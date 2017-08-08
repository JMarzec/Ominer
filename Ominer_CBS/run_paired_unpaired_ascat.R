##########################################This function runs ASCAT it takes as input files of BAF values and LRR values ###################
###############input arguments are (1) input - this takes as input a .txt file called input_ascat.txt, this file contains the LRR file, normalised_raw.txt
########and BAFR_input.txt this is the file containing BAF values. Each file should appear on a separate newline  (2) analysisType - paired/unpaired, (3) platform - type of ######Affymetrix platform of user data ###### (4)  type - CEL or log2ratio (ASCAT can run from log2ratio data - howver this option was not incorporated into O-miner , users can only ######use the ASCAT pipeline if they upload raw CEL file data. 
my_run_paired_unpaired_ascat <- function (input, analysisType, platform,type) 
{
	library("ominer")
	print ("This is where I am now")
	where <- getwd()
	print (where)
	print ("This is my input")
	print (input)
    inputs <- read.table(input, sep = "\t", as.is = TRUE)
    print ("this is my analysis type")
    print (analysisType)
    print ("this is my platform")
    print (platform)
    
    if (platform == "50xba") {
		platform = "Affy100k"
	}
        if (platform == "100") {
            platform = "Affy100k"
            
        }
	if (platform == "250sty") {
		platform <- "Affy250k_sty"
	}
	if (platform == "250nsp") {
            platform <- "Affy250k_nsp"
        }
        if (platform == "500k") {
            plaform <- "Affy500k"
        }
        if (platform == "SNP6") {
            platform <- "AffySNP6"
        }
print ("platform")
print (platform)
	
    
    
    if (analysisType == "unpaired") {
        LRR_all = inputs[1, ]
        BAF_all = inputs[2, ]
        ascat.bc = ascat.loadData(LRR_all, BAF_all)
      
        ascat.gg = ascat.predictGermlineGenotypes(ascat.bc, platform)
        here <- getwd()
             
        if (type == "CEL") {
        source("../../../../aspcf.R") 
        #print ("I found aspcf.R")
        source("../../../../ascat.aspc.R")
        }
        if (type == "log2ratio") {
            source("aspcf.R")
            source("ascat.aspc.R")
        }
        ascat.bc = ascat.aspc(ascat.bc, ascat.gg = ascat.gg)   ###the problem is here- aspcf.R is being called as a source within the ascat.aspcf function in the ominer library - take this function out of here and call as a source?
    }
    if (analysisType == "paired") {
        LRR_tumour = inputs[1, ]
        BAF_tumour = inputs[2, ]
        LRR_normal = inputs[3, ]
        BAF_normal = inputs[4, ]
        ascat.bc = ascat.loadData(LRR_tumour, BAF_tumour, LRR_normal, 
            BAF_normal)
        ascat.plotRawData(ascat.bc)
        ascat.plotRawData(ascat.bc)
        if (type == "CEL") {
        source("../../../../ascat.aspc.R")
        source("../../../../aspcf.R")
        }
        if (type == "log2ratio") {
         source("aspcf.R")
            source("ascat.aspc.R")
        }
        ascat.bc = ascat.aspc(ascat.bc)
        
    }
    ascat.plotSegmentedData(ascat.bc)
    ascat.output = ascat.runAscat(ascat.bc)
    str(ascat.output)
    tumour_content = ascat.output$aberrantcellfraction
    tumour_content
    tumour_ploidy = ascat.output$ploidy
    tumour_ploidy
    
    tumour_content_file <- paste ("../../../../../www/cgi-bin/onlinetool/version_2/ominer_results/",project, "/", "output", "/tumour_content.txt",sep="")
           
     write.table(tumour_content, file = tumour_content_file, 
        sep = "\t", quote = F)    
        
      
      tumour_ploidy_file <-   paste ("../../../../../www/cgi-bin/onlinetool/version_2/ominer_results/",project, "/", "output", "/tumour_ploidy.txt",sep="")
             write.table(tumour_ploidy, file = tumour_ploidy_file, sep = "\t", 
        quote = F)

        tumour_content_plot <- paste ("../../../../../www/cgi-bin/onlinetool/version_2/ominer_results/",project, "/", "output", "/tumour_content_plot.pdf",sep="")
        
   

############Added code not to do tumour content plots OR tumour ploidy plots here if the platform used is SNP6.......	
	if (platform == "SNP6") {
       pdf(file = tumour_content_plot, width = 6, height = 6, 
        pointsize = 11, bg = "white")
  
        
        
    plot(sort(ascat.output$aberrantcellfraction))
    dev.off()
    }
    
    if (platform == "SNP6") {
    tumour_ploidy_plot  <- paste ("../../../../../www/cgi-bin/onlinetool/version_2/ominer_results/",project, "/", "output", "/tumour_ploidy_plot.pdf",sep="")

    
    
  
        
       pdf(file = tumour_ploidy_plot, width = 6, height = 6, 
        pointsize = 11, bg = "white")
  
        
        
        
    plot(density(ascat.output$ploidy))
    dev.off()
    }
    
    
    ascat.segments <- organize.ascat.segments(ascat.output, ascat.bc$SNPpos)
    ascat.segments.tcn <- cbind(ascat.segments, ascat.segments[, 
        6] + ascat.segments[, 7])
    colnames(ascat.segments.tcn)[8] <- "tcn"
    ploidy <- ascat.output$ploidy
    names(ploidy) <- colnames(ascat.output$nA)
    ascat.segments.tcn <- cbind(ascat.segments.tcn, ploidy[match(ascat.segments.tcn[, 
        1], names(ploidy))], ascat.segments.tcn[, 6] - ploidy[match(ascat.segments.tcn[, 
        1], names(ploidy))])
    colnames(ascat.segments.tcn)[9:10] <- c("ploidy", "correctedTcnForPloidy")
    thrPloidy = 0.6
   if(is.na(threshold)) {
 		thr=1
                }
                else {
                    print ("Using user defined thresholds")
                    thr =as.numeric(threshold)
                    
                    }

    segments <- ascat.segments.tcn
    size <- segments$End - segments$Start
    segments <- data.frame(SampleID = segments$SampleID, chr = segments$Chr, 
        start = segments$Start, end = segments$End, size = size, 
        nSNP = segments$nProbes, segments[, -(1:5)])
    CNEventTypePloidyCorrected <- rep("None", dim(segments)[1])
    CNEventType <- rep("None", dim(segments)[1])
    LOH <- rep("None", dim(segments)[1])
    LOH[which(segments$nA == 0 | segments$nB == 0)] <- "LOH"
    CNEventTypePloidyCorrected[which(segments$correctedTcnForPloidy <= 
        (-thrPloidy))] <- "Loss"
    CNEventTypePloidyCorrected[which(segments$correctedTcnForPloidy >= 
        thrPloidy)] <- "Gain"
    CNEventTypePloidyCorrected[which(segments$correctedTcnForPloidy < 
        thrPloidy & segments$correctedTcnForPloidy > (-thrPloidy) & 
        (segments$nA == 0 | segments$nB == 0))] <- "CN_LOH"
    CNEventType[which((segments$tcn - 2) <= (-thr))] <- "Loss"
    CNEventType[which((segments$tcn - 2) >= thr)] <- "Gain"
    CNEventType[which((segments$tcn - 2) < thr & (segments$tcn - 
        2) > (-thr) & (segments$nA == 0 | segments$nB == 0))] <- "CN_LOH"
    LOH <- rep("None", dim(segments)[1])
    LOH[which(segments$nA == 0 | segments$nB == 0)] <- "LOH"
    segments_final <- data.frame(segments, CNEventType = CNEventType, 
        CNEventTypePloidyCorrected = CNEventTypePloidyCorrected, 
        LOH = LOH)
        
        segments_out_file <- paste("../../../../../www/cgi-bin/onlinetool/version_2/ominer_results/",project, "/", "Regions", "/segments_final_paired.out.txt",sep="")
    if (analysisType == "paired") {
        
            
             write.table(segments_final, file = segments_out_file, 
            sep = "\t", quote = F)

    }
    if (analysisType == "unpaired") {
    	
    	segments_out_file <- paste("../../../../../www/cgi-bin/onlinetool/version_2/ominer_results/",project, "/", "Regions", "/segments_final_unpaired.out.txt", sep="")
       
            write.table(segments_final, file = segments_out_file, 
            sep = "\t", quote = F)  
            
            
            
            
    }
    
    if (type == "CEL") {
    raw_LRR_dir <- paste("../../../../../www/cgi-bin/onlinetool/version_2/ominer_results/",project, "/","FP",sep="")
    #dir.create(raw_LRR_dir)
    current_dir <- getwd()
    system ("rm sunriseLRR_*.png")
    #system ("rm GermlineLRR_*.png")
    system ("rm rawprofileLRR_*.png")
    #system ("rm TumorLRR_*.png")
    
    
    
    output_files = list.files(".","*.png")
    for (outputfile in output_files)
    {
       	file.copy(paste(current_dir,outputfile,sep="/"),paste(raw_LRR_dir,outputfile,sep="/"))
           }

   
   }
   
   
   
    
    
    
        
    
}
