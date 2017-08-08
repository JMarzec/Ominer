############################################################################################################
#
#affy_correlation_code.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to calculate the Pearson's correlation coefficient from either a gene symbol or a probeid for each of the expression platforms supported by O-miner
###Arguments - (1) probeFile = a.txt file containing the probeid or gene symbol read in from the interface and written to a .txt file, (2) study - this is the project name
####(3) type = this is either symbol or probe, (4) platform = this indicates the platform that is used e.g. 450k or 27k for methylation (5) technology - this is the technol
####ogy that was used i.e. methylation or affy_expr, illumina_expr, affy_mirna and rna_seq, filename this is the name that is to be given to the output file
####Code to calculate correlation  a gene of interest (picked by probe_id) against all other genes on the array 
###generates a matrix of three columns of gene = "affy_probeid", cor = "correlation coefficient" and pvalue = "probability _value)
####user enters probeid/gene of interest and results are shown to user as a table of the top twenty genes that are creeltaed with yours and also 
###reasults of all of the correlation of all of the genes in one table to downlaod as an .xls sheet.

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
library("annotate")
#library("hgu133plus2.db")
output_dir <- "/var/www/html/onlinetool/temp/"

####read in the .txt file read in from the interface with the idenity of the probeid/gene_symbol
gene_file<-read.table(probeFile,as.is=T,header=F)
####Assign variable to the gene of interest i.e. probeid or gene_symbol
gene.of.interest <- gene_file$V1
if (platform == "450k"||platform == "27k") {
technology <- "methylation"
print ("my technology is methylation")
}

print (technology)
if (technology == "affy_expr") {
	####R library according to the affy-expression platform that is used is loaded up - this is used when the user inputs a gene_symbol and the corresponding probeid is looked up.
	library(paste(platform,".db",sep=""),character.only = TRUE)
	}
	
		gene_file<-read.table(probeFile,as.is=T,header=F)
gene.of.interest <- gene_file$V1



 
        
  
       
if (technology == "affy_expr") {
			 
####read in the gene expression matrix of probes that pass the filter 
       	A.data <- read.table(paste(rootdir,study,"/","norm","/","filtered_data.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1) 
       
        
	ghj <- type
	      
        }
        
    if (technology == "illumina_expr") {
    	
       	A.data <- read.table(paste(rootdir,study,"/","norm","/","filtered_data.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1) 
        print ("I have got illumina expression matrix")
              }
	
   if (technology == "methylation") {
    	 print ("I am running correlation:methylation")

       	A.data <- read.table(paste(rootdir,study,"/","norm","/","normalised.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1) 
        print ("I have got methylation matrix")
           
        }

    	
    	
    	
	if (technology == "affy_mirna") {
		  print ("I am HERE")

       	A.data <- read.table(paste(rootdir,study,"/","norm","/","filtered_data.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1) 
        print ("I am HERE")
        print (A.data)
        }
    
        if (technology == "rna_seq") {
        	A.data <- read.table(paste(rootdir,study,"/","norm","/","normalised.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1) 
        	        }
        
    if (technology == "affy_mirna") {
		  print ("I am HERE")

       	A.data <- read.table(paste(rootdir,study,"/","norm","/","filtered_data.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1) 
        print ("I am HERE")
        print (A.data)
      
       
       if (type == "symbol") {
       	####Looking up the corresponding gene symbol if a probe is used is different for miRNA platforms as the symbol is the same as the probeid with a "_st" added on 
       gi <- paste(gene.of.interest, "_st",sep="")
        gene.index <- grep(gi,rownames(A.data))	
       }
        }
    if (type == "probe") {
    	####If a probe is used then this is looked up for in the rowid columns
            gene.index <- grep(gene.of.interest,rownames(A.data))	

            }

	if (technology == "affy_expr") {
	######If any affy_expression platform is used create a dataframe of probeid and gene symbol by using the probeids to map back to gene_symbol
	affy_probeids <- rownames(A.data)
	  probeid = unique(affy_probeids)
	  symbol <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
            "SYMBOL", sep = ""))), function(symbol) {
            return(paste(symbol, collapse = "; "))
        })))
               
                        
        
        NAME <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
            "GENENAME", sep = ""))), function(name) {
            return(paste(name, collapse = "; "))
        })))
####AnnotA dataframe is created with two columns (1) probeid, (2) symbol
        AnnotA <- data.frame(probeid, symbol, 
            NAME, row.names = NULL)
 

affy_wanted <- probeid
            if (type == "symbol") {
            	####Extract the indice of the gene_symbol of interest from dataframe AnnotA
            gene.index <- grep(gene.of.interest,AnnotA$symbol)

	    }
            else {
	     gene.index <- grep(gene.of.interest,rownames(A.data))	
	     ####If it is a probeid, extract the indice of the probeid of interest from dataframe AnnotA
	 print ("THIS IS THE GENE INDEX")
print (gene.index)   
	    }
            gene_symbols <- AnnotA$symbol
          
         
            
           	}
		
		
		if (technology == "methylation") {
		meth_probeids <- rownames(A.data)
		if (platform == "450k"){
        	
        
        }
               
####If dataset is from the methylation technology - read in a .txt file of probeids and gene symbols 

	     AnnotA <- read.table("methy_gs_lookup.txt",header=TRUE,sep="\t",as.is=T)
	    colnames(AnnotA) <- c("probe_id","symbol")
	    
	     if (type == "symbol") {
            gene.index <- grep(gene.of.interest,AnnotA$symbol)
            ####Extract the indice of the gene_symbol of interest from dataframe AnnotA
	    }
else {
	gene.index <- grep(gene.of.interest,rownames(A.data))	
	 ####If it is a probeid, extract the indice of the probeid of interest from dataframe AnnotA
RNA <- rownames(A.data)
            }
	   
	    
	    
	    
	    
	    
		}
				
####If dataset is from Illumina_expression technology need to create a dataframe of illumina probeids and gene symbol by using the probeids as input to BiomaRt  
      if (technology == "illumina_expr") {   ###here need to get gene symbols for ALL of teh probes onn the array
	######map the affyprobeids to gene symbol
	illum_probeids <- rownames(A.data)
	  probeid = unique(illum_probeids)
	  library("biomaRt")
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")

if (platform == "ht12v3") {
	platform_annot = "illumina_humanht_12_v3"
}
if (platform == "ht12v4") {
	platform_annot = "illumina_humanht_12_v4"
}

theFilters = c(platform_annot)   ######or illumina_humanht_v3
theAttributes = c(platform_annot,"hgnc_symbol")
wanted_annot <- getBM(attributes = theAttributes, filters = theFilters, values = rownames(A.data), mart=mart)
		
		
		illum_genesymbols <- wanted_annot[2]
		colnames(wanted_annot) <- c("Illum_probe_id","Symbol")
		  gene.index <- grep(gene.of.interest,wanted_annot$Symbol)
		 
            if (type == "probe") {
            gene.index <- grep(gene.of.interest,rownames(A.data))	

            }
	    
print (" I got Illumina_gene_symbols")
             }
             
             
 ####If dataset is RNA_Seq use the ENSEMBL gene_ids which are the rownames in filtered data to create a dataframe with two columns i.e. Ensembl gene id and HGNC Symbol          
	if (technology == "rna_seq") {
	print ("I am doing this extracting data from bioMart")
		#pl <- probeList$V1
	
		tester <- grepl("^ENSG",rownames(A.data))
     if (tester[1] == "TRUE") {
		

		library("biomaRt")
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
theFilters = c("ensembl_gene_id")   ######or illumina_humanht_v3
theAttributes = c("ensembl_gene_id","hgnc_symbol")
wanted_annot <- getBM(attributes = theAttributes, filters = theFilters, values = rownames(A.data), mart=mart)
rna_seq_gene_symbols <- wanted_annot[2]
colnames(wanted_annot) <- c("ensembl_gene_id","Symbol")
####If type is probe - extract the indice required from the "normalised" data 
if (type == "probe") {
 gene.index <- grep(gene.of.interest,rownames(A.data))
 }
 if (type == "symbol") {
####If type is symbol - extract the indice required from the dataframe - wanted_annot

 gene.index <- grep(gene.of.interest,wanted_annot$Symbol)
 }
			

}

###Code below id required if the reads counts matrix has gene symbols as rownames. For RNA-seq post-processing pipeline users can either upload a matrix with Ensembl gene_ids ad the rownames or of one with 
###gene symbols as rownames 
if (tester[1] == "FALSE") {   
library("biomaRt")
##wanted<-read.table("illumina_probes_annot.txt",header=T,sep="\t",as.is=T)
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
theFilters = c("hgnc_symbol")   ######or illumina_humanht_v3
theAttributes = c("ensembl_gene_id","hgnc_symbol")
wanted_annot <- getBM(attributes = theAttributes, filters = theFilters, values = rownames(A.data), mart=mart)
rna_seq_gene_symbols <- wanted_annot[2]
colnames(wanted_annot) <- c("ensembl_gene_id","Symbol")
###If symbol is input - extract the indice of the symbol from the "normalised matrix"
if (type == "symbol") {
 gene.index <- grep(gene.of.interest,rownames(A.data))
 
 }
 
 ###If probe id input extract the required indice from the data farme containing enesembl_gene_id and HGNC symbol
 if (type == "probe") {
 gene.index <- grep(gene.of.interest,wanted_annot$ensembl_gene_id)
 }





}


  }
             
                   	
            
           
            		if(technology == "affy_exon") {
            		            		 library("biomaRt")
			 mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
              theFilters = c("affy_huex_1_0_st_v2")
  
 theAttributes =c("affy_huex_1_0_st_v2","hgnc_symbol")

ampGenes <- getBM(attributes = theAttributes, filters = theFilters, values = A.data$groupName, mart=mart)
gene_symbols <- ampGenes$hgnc_symbol
 u <- unique(gene_symbols)           		
            		
            			
           
            	colnames(ampGenes) <- c("probe_id","Symbol")
            	gene.index <- grep(gene.of.interest,ampGenes$Symbol)
            	if (type == "probe") {
            	gene.index <- grep(gene.of.interest,A.data$unitName)	
      }      	
      }
            	 if (technology == "affy_expr") {
            	 gene_symbols <- AnnotA$symbol
            	 }
            	              
sorted.data <- A.data[ sort(rownames(A.data) ) , ];
all_affy_probe_ids <- rownames(sorted.data)
	mat_length <- nrow(A.data)
	if (technology == "affy_exon") {
		#mat_length <- length(u)
		#n <-A.data[1:18704,]
	u <- rownames(A.data)
	
	
data.output <- matrix( 
	
	data = 0, nrow = mat_length, ncol = 3, dimnames = list( u, c("probeid","Correlation", "P-value") )
);

}

###Setting up the structure to contain the gene_symbols, probeids and correlation coefficients for each of the different technologies
if (technology == "affy_mirna") {
	###Get the number of probeids/symbols present in the normalised matrix - this allows us to set the size of the correlation matrix that is to be generated.
		mat_length <- nrow(A.data)
		###Assign the probeids to a value - this is to annotate the correlation matrix that will be generated.
	u <- rownames(A.data)
	
	###Set the dimnensions of the correlation matrix - with the number of rows equal to the number of rows in the "normalised" matrix
	data.output <- matrix( 
	
###Assign the column names to the correlation matrix i.e. probeid, correlation and p-value 
	data = 0, nrow = mat_length, ncol = 3, dimnames = list( u, c("probeid","Correlation", "P-value") )
);
###Assign gene_symbols to the first column of the "correlation matrix"
data.output[,1] <- rownames(sorted.data)
}


###Repeat for correlation matrix for data derived from affy_expression platforms 

if (technology == "affy_expr") {
	data.output <- matrix( 
	data = 0, nrow = length(all_affy_probe_ids), ncol = 3, dimnames = list( all_affy_probe_ids, c("probeid", "Correlation", "P-value") )

	);

print ("THIS IS PROBEIDS")

data.output[,1] <- affy_wanted
###Assign gene_symbols to the first column of the "correlation matrix"

rownames(data.output) <- gene_symbols		
		
		
		
	
	
}

###Repeat for correlation matrix for data derived from methylation platforms 

if (technology == "methylation") {
	
data.output <- matrix( 
	
data=0,nrow=1000,ncol=3,dimnames=list(all_affy_probe_ids[1:1000],c("Gene","Correlation","P-value"))
	);
row_all <- rownames(sorted.data)
row_100 <- row_all[1:1000]
data.output[,1] <- row_100

}



###Repeat for correlation matrix for data derived from Illumina_expression platforms 

if (technology == "illumina_expr") {
			u <- rownames(A.data)
		mat_length <- nrow(A.data)
data.output <- matrix( 
	data = 0, nrow = mat_length, ncol = 3, dimnames = list( u, c("Gene", "Correlation", "P-value") )

	);
data.output[,1] <- rownames(sorted.data)
		
		
		
		
	#}
	print ("I generate data output matrix")
}




###Repeat for correlation matrix for data derived for the post-processing of RNA-Seq data  
if (technology == "rna_seq") {
			u <- rownames(A.data)
		mat_length <- nrow(A.data)
data.output <- matrix( 
	data = 0, nrow = mat_length, ncol = 3, dimnames = list( u, c("Gene", "Correlation", "P-value") )

	);
data.output[,1] <- rownames(sorted.data)
		
		
		
		
	#}
	print ("I generate data output matrix")
}




gl_length <- length(gene.index)
if (gl_length > 1) {
gene.index <- gene.index[1]
}
if (technology == "affy_exon") {
	x <- as.matrix(A.data)
	n <- A.data[1:18704,]
	m <- log2(n)
	for ( i in 1:nrow(data.output) ) {

	Obsv <- cor.test( as.numeric(A.data[ gene.index, ]), as.numeric(log2(sorted.data[ i, ])), method = "pearson" );

	data.output[i, "Correlation"] <- Obsv$estimate;
	data.output[i, "P-value"] <- Obsv$p.value;
	}
data.output[,1] <- A.data$groupName
}

if (technology == "affy_expr") {
	

for ( i in 1:nrow(data.output) ) {

	Obsv <- cor.test( as.numeric(A.data[ gene.index, ]), as.numeric(sorted.data[ i, ]), method = "pearson" );

	data.output[i, "Correlation"] <- Obsv$estimate;
	data.output[i, "P-value"] <- Obsv$p.value;
	}
	}
	

###Calculate the Pearson's correlation coefficient for each of the probeids 
if (technology == "illumina_expr") {
	
###for each of the gene symbols in the empty correlation table 
for ( i in 1:nrow(data.output) ) {
###Calculate the correlation_coefficient for each of the probeids present in the dataset against the probeid/gene _symbol of interest and populate the correlation table 
	Obsv <- cor.test( as.numeric(A.data[ gene.index, ]), as.numeric(sorted.data[ i, ]), method = "pearson" );

	data.output[i, "Correlation"] <- Obsv$estimate;
	data.output[i, "P-value"] <- Obsv$p.value;
	}
	print ("I have filled the matrix")
	}	
	
###Repeat for the methylation technology
	
	if (technology == "methylation") {
	

for ( i in 1:nrow(data.output) ) {

	Obsv <- cor.test( as.numeric(A.data[ gene.index, ]), as.numeric(sorted.data[ i, ]), method = "pearson" );

	data.output[i, "Correlation"] <- Obsv$estimate;
	data.output[i, "P-value"] <- Obsv$p.value;
	}
	print ("I have filled the matrix")
	}	
	
###Repeat for miRNA
	if (technology == "affy_mirna") {
	

for ( i in 1:nrow(data.output) ) {

	Obsv <- cor.test( as.numeric(A.data[ gene.index, ]), as.numeric(sorted.data[ i, ]), method = "pearson" );

	data.output[i, "Correlation"] <- Obsv$estimate;
	data.output[i, "P-value"] <- Obsv$p.value;
	}
	}
	
###Repeat for RNA-Seq

	if (technology == "rna_seq") {
	
for ( i in 1:nrow(data.output) ) {

	Obsv <- cor.test( as.numeric(A.data[ gene.index, ]), as.numeric(sorted.data[ i, ]), method = "pearson" );

	data.output[i, "Correlation"] <- Obsv$estimate;
	data.output[i, "P-value"] <- Obsv$p.value;
	}
	}

	
	
###Repeat for affy_expression 

	if (technology == "affy_expr") {
	rownames(data.output) <- gene_symbols
	gs <- as.character(gene_symbols)
	new_data_matrix <- cbind(data.output[,1],gs,data.output[,2],data.output[,3])
	colnames(new_data_matrix) <- c("probe_id","symbol","Correlation","P-value")
	}
	
	
###Repeat for RNA-post-processing pipeline	
	
	if (technology == "rna_seq") {
		colnames(wanted_annot) <- c("probe_id","symbol")
		colnames(data.output) <- c("probe_id","correlation","pvalue")
				new_data_matrix <- merge(wanted_annot,data.output,by.x="probe_id",by.y= "probe_id")
colnames(new_data_matrix) <- c("probe_id","symbol","Correlation","P-value")
	print (" I have the data matrix")
	
	}
	
###Repeat for methylation 

	if (technology == "methylation") {
		colnames(AnnotA) <- c("probe_id","symbol")
		colnames(data.output) <- c("probe_id","correlation","pvalue")
				new_data_matrix <- merge(AnnotA,data.output,by.x="probe_id",by.y= "probe_id")
colnames(new_data_matrix) <- c("probe_id","symbol","Correlation","P-value")
	print (" I have the data matrix")
	
	}
	
###Add a column containing the probeids to each of the correlation matrices generated from each of the technologies.

	if (technology == "affy_mirna") {
	
		fgh <- strsplit(all_affy_probe_ids,"_st")
		ap <- unlist(fgh)
		new_data_matrix <- cbind(data.output[,1],ap,data.output[,2],data.output[,3])
				colnames(new_data_matrix) <- c("probe_id","symbol","Correlation","P-value")
			}
	if (technology == "affy_exon") {
			colnames(ampGenes) <- c("probe_id","symbol")
		colnames(data.output) <- c("probe_id","correlation","pvalue")
		
		new_data_matrix <- merge(ampGenes,data.output,by.x="probe_id",by.y= "probe_id")	
	}
	if (technology == "illumina_expr") {
		colnames(wanted_annot) <- c("probe_id","symbol")
		colnames(data.output) <- c("probe_id","correlation","pvalue")
		
		new_data_matrix <- merge(wanted_annot,data.output,by.x="probe_id",by.y= "probe_id")
		
		print ("I have annotated the matrix")
			colnames(new_data_matrix) <- c("probe_id","symbol","Correlation","P-value")
			}
	
	
	
			
 ###write the correlation matrix for each of the technologies to a .txt file

write.table(new_data_matrix,paste(output_dir,study,"/",filename,sep=""), sep = "\t", quote = FALSE, row.names=FALSE,
         col.names = TRUE)
###Sort the correlation table by p-values to generate the top ten and write these to a .txt file - top ten are those that are displayed on the interface
         
new_sort <- new_data_matrix[order(new_data_matrix[,4]),]    ######sort by p-value and show the top ten on webpage

sorted_10 <- head(new_sort,10)
colnames(sorted_10) <- c("probe_id","symbol","Correlation","P-value")
top_10_filename <- paste(filename,"_top_10.txt",sep="")
write.table(sorted_10,paste(output_dir,study,"/",top_10_filename,sep=""), sep = "\t", quote = FALSE, row.names=FALSE,
         col.names = TRUE) #}
print ("I am printing OUT the correlations NOW")
