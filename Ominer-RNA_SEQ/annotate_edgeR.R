##########################################################################################
#
#File name: annotate_edgeR.R
#
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#
#
#
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
##########################################################################################
###########################################################################################
#####Input arguments:
#####comp - this is the .txt file containing information regarding the comparisons to be made for each dataset analysed
#####project - name of the analysis
#####This function annotates lists of differentailly expressed genes from analysis with edgeR using either HGNC symbols as input or ENSEMBL gene ids (whichever is found/used as the identifier in the list(s) of differentially expressed genes)
annotated_edgeR <- function(comp,project) {
comparisons <- read.table(comp, sep = "\t", as.is = T, header = F, 
        strip.white = T)
 for (i in 1:length(comparisons)) {
    	print(comparisons[i])
    	        fg <- strsplit(comparisons[,i],"=")
        ci <- fg[[i]][1]

path <- paste("ominer_results/",project,"/","DifferentialExpression/",ci,"_unannotated.txt",sep="")
de_results  <- read.table(path, sep="\t", as.is=TRUE,strip.white=TRUE,header=T);
tester <- grepl("^ENSG",rownames(de_results))
     if (tester[1] == "TRUE") {
ensgids <- rownames(de_results)
cds <- cbind(ensgids,de_results)
library("biomaRt")
up = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
print ("I am here")
            #up <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
	    print ("I am here now 2")
  theFilters = c("ensembl_gene_id")
theAttributes =c("ensembl_gene_id","chromosome_name", "band","hgnc_symbol")
ampGenes <- getBM(attributes = theAttributes, filters = theFilters, values = cds$ensgids, mart=up)
#new <- cbind(ampGenes$chromosome_name,ampGenes$band)
new_col <- paste(ampGenes$chromosome_name, ampGenes$band,sep="")
together_new <- cbind(ampGenes$ensembl_gene_id,new_col,ampGenes$hgnc_symbol)
colnames(together_new) <- c("ensembl_gene_id","locus","hgnc_symbol")

gene_filtered = merge(together_new,cds,by.x="ensembl_gene_id",by.y="ensgids",all=TRUE)
#selected_filtered <- cbind(gene_filtered$ensembl_gene_id, gene_filtered$locus,gene_filtered$hgnc_symbol,gene_filtered$logFC, gene_filtered##$logCPM,gene_filtered$FDR)
selected_filtered <- gene_filtered[,c(1,2,3,4,8)]

unique_all <-selected_filtered[!duplicated(selected_filtered$ensembl_gene_id), ] 

write.table(unique_all, paste("ominer_results/",project,"/","DifferentialExpression/",ci,"_Filtered.txt",sep=""), sep = "\t", row.names = FALSE, 
        quote = FALSE)
	
	write.table(unique_all, paste("ominer_results/",project,"/","DifferentialExpression/",ci,"annotated.txt",sep=""), sep = "\t", row.names = FALSE, 
        quote = FALSE)
        }
	#####Annotate genes with the gene symbol as rownames identifier
	
	else {                                      
	
	gene_symbols <- rownames(de_results)
gds <- cbind(gene_symbols,de_results)
library("biomaRt")
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
theAttributes = c("hgnc_symbol","ensembl_gene_id","band")



theFilters = c("hgnc_symbol")   

wanted_annot <- getBM(attributes = theAttributes, filters = theFilters, values = x$ID, mart=mart)
        allA_symbol <- merge(x,wanted_annot,by.x="ID",by.y="hgnc_symbol")
    format_AllA <- cbind(allA_symbol$ensembl_gene_id,allA_symbol$band,allA_symbol$ID,allA_symbol$logFC,allA_symbol$adj.P.Val)    
        colnames(format_AllA) <-c("probe_id","locus","symbol","logFC","adjusted_p_value")
unique_all <-selected_filtered[!duplicated(format_AllA$probe_id), ] 

write.table(unique_all, paste("ominer_results/",project,"/","DifferentialExpression/",ci,"_Filtered.txt",sep=""), sep = "\t", row.names = FALSE, 
        quote = FALSE)
	
write.table(unique_all, paste("ominer_results/",project,"/","DifferentialExpression/",ci,"annotated.txt",sep=""), sep = "\t", row.names = FALSE, 
        quote = FALSE)	
	}
	
	}
        }