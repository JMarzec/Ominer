##########################################################################################
#
#File name: miRNA_ominer.R
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
#Description: Read in CEL files from Affymetrix miRNA array platforms
#Arguments are:
#target_file is the full path to the target file
#project is the given project name
##########################################################################################

run_miRNA <- function(target,project) {


	#library("simpleaffy")
	# rma normalisation
#if (norm == "rma") {
eset <- rma(affybatch)
#}

ID <- featureNames(eset)

# select human miRNA and create expressionset
#grep( "hsa", featureNames(eset) ) -> hsa; 
#eset.hsa <- eset[hsa , ];

#hsa.ID <- featureNames(eset.hsa);

exp.data <- exprs(eset);

# save rma file
pd<-read.table(target_file,header=T,sep="\t")
colnames(exp.data) <- pd$Name
write.table(exp.data,paste("ominer_results/",project,"/","norm","/normalised.txt",sep="") ,sep="\t",row.names = TRUE, quote = FALSE);

# ANNOTATION
# using biomaRt lib specify which "mart" to use
#ensembl <- useMart("ensembl");
#ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl);

# attributes available for download
# listAttributes(ensembl);

# retrieve user-specified attributes from BioMart db, in this instance miRNA annotations
#miRNA <- getBM( 
 #       attributes = c("mirbase_id", "ensembl_gene_id", "start_position", "chromosome_name"), 
 #       filters = c("with_mirbase"), 
   #     values = list(TRUE), 
  #      mart = ensembl
#);


}