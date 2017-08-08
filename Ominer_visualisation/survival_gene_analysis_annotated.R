############################################################################################################
#
#survival_gene_analysis.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to calculate perform survival analysis for particular gene(s) of interest for each of the expression platforms supported by O-miner
###Arguments - (1) probeFile = a.txt file containing the probeid or gene symbol read in from the interface and written to a .txt file, (2) study - this is the project name
####(3) type = this is either symbol or probe, (4) platform = this indicates the platform that is used e.g. 450k or 27k for methylation (5) technology - this is the technol
####ogy that was used i.e. methylation or affy_expr, illumina_expr, affy_mirna and rna_seq, filename this is the name that is to be given to the output file
####Code to calculate a Kapla-meier (KM) plot for a gene of interest (picked by probe_id) against the expresiion of other probes on the array 
####KM_plot is generated in both .png and .pdf format - .png is displayed on the webpage and .pdf is available to downlaod from the webpage. 
####Code is dependent on (1) FC14.plotting.lib, (2) utility_methods.R and (3)  bcctb.utils.R 
####.tiff file and .pdf is generated of the KM plot - .tiff file is displayed within the webpage and .pdf is ab=vailable to download from the results page 
	
	for (e in commandArgs()) {
	ta = strsplit(e,"=",fixed=TRUE)
	ta
	if(! is.na(ta[[1]][2])) {
		temp = ta[[1]][2]
		assign(ta[[1]][1],temp)
		cat("assigned ",ta[[1]][1]," the value of |",temp,"|\n")
	} 
}
rootdir<-"/var/www/cgi-bin/onlinetool/version_2/ominer_results/"
output_dir <- "/var/www/html/onlinetool/temp/"

library("FC14.plotting.lib");
source("utility_methods.R");
source("bcctb.utils.R")

####Read in probeFile - this is the .txt file in which the probeid/gene symbol is read in and written to from the interface
gene_file<-read.table(probeFile,as.is=T,header=F)
gene.i <- gene_file$V1

####Read in the normalised matrix for each of the technologies suppoerted by O-miner
if (technology == "affy_expr") {
library(paste(platform,".db",sep=""),character.only = TRUE)
			
       	A.data <- read.table(paste(rootdir,study,"/","norm","/","normalised.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1) 
       
	  A.data.entrez<-getEntrez(A.data,platform,"SYMBOL")
	  A.data.entrez <- t( scale( t(A.data.entrez) ) );

        }
        
    if (technology == "illumina_expr") {
    	
       	A.data <- read.table(paste(rootdir,study,"/","norm","/","normalised.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1) 
        
                }

 if (technology == "affy_mirna") {
	
       	A.data <- read.table(paste(rootdir,study,"/","norm","/","filtered_data.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1) 
               }
        
          if (technology == "rna_seq") {
        	A.data <- read.table(paste(rootdir,study,"/","norm","/","normalised.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1) 
        	
				library("biomaRt")

ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
theFilters = c("ensembl_gene_id")   ######or illumina_humanht_v3
theAttributes = c("ensembl_gene_id","hgnc_symbol")
wanted_annot <- getBM(attributes = theAttributes, filters = theFilters, values = rownames(A.data), mart=ensembl)
rna_seq_gene_symbols <- wanted_annot[2]
colnames(wanted_annot) <- c("ensembl_gene_id","Symbol")

if (type == "symbol") {
 gene.index <- grep(gene.i,wanted_annot$Symbol)
		 }
           

        }

    
if (technology == "affy_exon") {
	A.data<-read.table(paste(rootdir,study,"/transcript/norm/normalised.txt",sep=""),sep="\t",as.is=T,header=T,row.names=1)   #####added for exon_array as A.data to be used is from teh transcript level
	  #####
 rownames(A.data) <- A.data$unitName

            		#gene.index <- grep(gene.of.interest,A.data$unitName)
            		 library("biomaRt")
			 mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
            #up <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  theFilters = c("affy_huex_1_0_st_v2")
  theAttributes =c("affy_huex_1_0_st_v2","hgnc_symbol")

ampGenes <- getBM(attributes = theAttributes, filters = theFilters, values = A.data$groupName, mart=mart)
gene_symbols <- ampGenes$hgnc_symbol
 u <- unique(gene_symbols)           		
            		
            			
            
            	            	colnames(ampGenes) <- c("probe_id","Symbol")
            	if (type == "symbol") {
            	gene.index <- grep(gene.i,ampGenes$Symbol)
     }
           
           }
           
           



	
	   if (technology == "illumina_expr") {   ###here need to get gene symbols for ALL of the probes on the array
	######map the illumina_probeids to gene symbol with biomaRt
	illum_probeids <- rownames(A.data)
	  probeid = unique(illum_probeids)
	  library("biomaRt")
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
if (platform == "ht12v3") {
	platform_annot = "illumina_humanht_12_v3"
}
if (platform == "ht12v4") {
	platform_annot = "illumina_humanht_12_v4"
}

theFilters = c(platform_annot)   ######or illumina_humanht_v3
theAttributes = c(platform_annot,"hgnc_symbol")
wanted_annot <- getBM(attributes = theAttributes, filters = theFilters, values = rownames(A.data), mart=ensembl)
		
		if (type == "symbol") {
		colnames(wanted_annot) <- c("Illum_probe_id","Symbol")
		
		  gene.index <- grep(gene.i,wanted_annot$Symbol)
		 }
                       }
print (" I got Illumina_gene_symbols")
             
 
####assigning functions used to calculate survival 	

# function used dichotomise exp data on median
local.dichotomise.dataset <- function(x, split_at = 99999) {
  if (split_at == 99999) { split_at = median(x, na.rm = TRUE); }
  return( as.numeric( x > split_at ) );
}

# function to round off values in first arg to specified dp
my.round <- function(number){
  rounded.value <- round(number, digits = 3);
  return(rounded.value);
}




  revised <- "0"
  
####For each of the technologies - read in the normalised expression matrix and traget file - note that if the user wishes to perform survival analysis on their data - two extra columns are required in the target file these 
####are (1) overall_survival and survival_status 0/1
  
  if (technology == "affy_mirna") {
			 
	exp.data <- as.matrix(read.table(file = paste("ominer_results/",study,"/norm/normalised.txt",sep=""),row.names = 1, header = T, sep = "\t"));
	ann.data <- as.matrix(read.table(file = paste("ominer_results/",study,"/QC/target.txt",sep=""),row.names = 1, header = T, sep = "\t"));

       	        }
       	        
      if (technology == "rna_seq") {
	
		  
	exp.data <- as.matrix(read.table(file = paste("ominer_results/",study,"/norm/normalised.txt",sep=""),row.names = 1, header = T, sep = "\t"));
	ann.data <- as.matrix(read.table(file = paste("ominer_results/",study,"/QC/target.txt",sep=""),row.names = 1, header = T, sep = "\t"));

       	        }
 	        
       	        
       	        
  if (technology == "affy_exon") {
  	exp.data <- as.matrix(read.table(file = paste("ominer_results/",study,"/transcript/norm/normalised.txt",sep=""),row.names = 1, header = T, sep = "\t"));
  	ann.data <- as.matrix(read.table(file = paste("ominer_results/",study,"/transcript/QC/target.txt",sep=""),row.names = 1, header = T, sep = "\t"));

totalcols <- ncol(exp.data)
	new_exp_dat <- exp.data[,c(6:totalcols)]
	exp.data <- new_exp_dat

  }
 if (technology == "affy_expr") {
out <- read.table(paste("ominer_results/",study,"/QC","/outliers.txt",sep=""), sep = "\t", as.is = T, 
        header = T)
    colnames(out) = c("Sample", "averageBG", "ScaleFactor", "Present", 
        "HSAC07", "GAPDH", "NUSE", "RLE", "RLE_IQR", "RNAslope")
    if (nrow(out) > 0) {
    	revised = "1"
    	}
    	 if (nrow(out) == 0){
    		revised = "0"
    	    		#file.copy(target_file,qcdir)
    	}
    	}
    	
   if (technology == "illumina_expr") {
out <- read.table(paste("ominer_results/",study,"/QC","/outliers.txt",sep=""), sep = "\t", as.is = T, 
        header = T)
    colnames(out) = c("mean", "standard deviation", "detection rate(0.01)", "distance to sample mean") 
       
    if (nrow(out) > 0) {
    	revised = "1"
    	}
    	 if (nrow(out) == 0){
    		revised = "0"
    	    		    	}
    	}
    	 	
    	
   
####For both Illumina_expression and Affymetrix_expression platforms need to take into account whether any samples did not pass quality control - if they did then the target file is not target.txt  - it is target_qc.txt 
  ####that needs to be read in  	

  if (technology == "affy_expr"){
  	exp.data <- as.matrix(read.table(file = paste("ominer_results/",study,"/norm/normalised.txt",sep=""),row.names = 1, header = T, sep = "\t"));

		if (revised == "1") {
    	    	ann.data <- as.matrix(read.table(file = paste("ominer_results/",study,"/QC/target_qc.txt",sep=""),row.names = 1, header = T, sep = "\t"));
    	}
    	
    	
    		if (revised == "0") {
    			ann.data <- as.matrix(read.table(file = paste("ominer_results/",study,"/QC/target.txt",sep=""),row.names = 1, header = T, sep = "\t"));

}
}
 if (technology == "illumina_expr") {
  	exp.data <- as.matrix(read.table(file = paste("ominer_results/",study,"/norm/normalised.txt",sep=""),row.names = 1, header = T, sep = "\t"));

		if (revised == "1") {
    	    	ann.data <- as.matrix(read.table(file = paste("ominer_results/",study,"/QC/target_qc.txt",sep=""),row.names = 1, header = T, sep = "\t"));
    	}
    	
    	
    		if (revised == "0") {
    			ann.data <- as.matrix(read.table(file = paste("ominer_results/",study,"/QC/target.txt",sep=""),row.names = 1, header = T, sep = "\t"));

}
}
 
    	

# scale exp.data
exp.data <- t(scale(t(exp.data)));
  

# UNIVARIATE SURVIVAL MODELING

    

# extract surv data
surv.time <- as.numeric( ann.data[ , "Surv_Period"] );
surv.stat <- as.numeric( ann.data[ , "Surv_Status"] );
        

if (type == "probe") {
gene.index <- grep(gene.i,rownames(exp.data))	
if (technology == "affy_exon") {
gene.index <- grep(gene.i,A.data$unitName)	
}
}
  
####Create the risk groups for survival analysis using local.dichtomise dataset 
if (type=="symbol" && technology =="affy_expr") {

rg <- local.dichotomise.dataset( A.data.entrez[ gene.i, ] );

}


if (type=="probe" && technology=="affy_expr") {
A.data <- as.matrix(A.data)

rg <- local.dichotomise.dataset( A.data[ gene.i, ] );

}

if (type=="probe" && technology=="affy_mirna") {
A.data <- as.matrix(A.data)
rg <- local.dichotomise.dataset( A.data[ gene.i, ] );
}

if (type=="symbol" && technology=="affy_mirna") {
A.data <- as.matrix(A.data)

gene.b <- paste(gene.i,"_st",sep="")
rg <- local.dichotomise.dataset( A.data[ gene.b, ] );

}

if (type=="probe" && technology=="affy_exon") {
A.data <- data.frame(A.data)

gene.index <- grep(gene.i,A.data$groupName)
rg <- local.dichotomise.dataset(as.numeric( A.data[ gene.index, ]) );
}

  
####using the Cox model generate the KM plots for 5yrs 

cox.fit <- summary( coxph( Surv(surv.time, surv.stat) ~ rg) ); 
        all.results <- c(
	"HR" = my.round(cox.fit$conf.int[1,1]), "CI95L" = my.round(cox.fit$conf.int[1,3]), "CI95U" = my.round(cox.fit$conf.int[1,4]), "WaldP" = cox.fit$coef[1,5] 
);

fn <- paste(output_dir,study,"/",study,"_",gene.i,"survplot.tiff",sep="")
# KM PLOTS
print ("I am NOW doing this")
time = "5" ####this if for 5 years
make.KM.curve(
                riskgroup = rg,
                survtime = surv.time,
                survstat = surv.stat,
                file.name = fn, 
                truncate.survival = as.numeric(time),
                cex.lab = 1.1,
                cex.axis = 1.0,
                main.title = paste(gene.i,sep = " " ),
                xaxis.label = "Survival time",
                yaxis.label = "Survival probability",
                resolution = 900
        );



fg_name <- paste(output_dir,study,"/",study,"_",gene.i,"survplot",sep="")

change_image_dir <- paste ("/var/www/html/onlinetool/temp/",study,sep="")
setwd(change_image_dir)
loc_dir <- getwd()


oldfn <- paste(fg_name,".tiff",sep="")
newfn <- paste(fg_name,".png",sep="")
try <-  paste(fg_name,".pdf",sep="")

command_png <- paste("convert",oldfn,"-resize 25%",newfn,sep=" ")
command_pdf <- paste("convert",oldfn,"-resize 25%",try,sep=" ")
system(command_png)
system(command_pdf)
setwd("../../../")

