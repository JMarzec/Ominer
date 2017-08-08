##########################################################################################
#
#File name: GOStats_function.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#
#
#
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
##########################################################################################
##########################################################################################
#Description: function to perform GOstats analysis on the results of differential expression analysis 
#project - name given to the analysis
#comp - full path to the comparisons file

################################################################################################KM plots - overall survival from a target file generates .txt file and KM plots for 5,10 and 15 years
######







#####GOstats function uses as an input the affyprobeids of the probes that pass the logfc and pvalue cutoffs
######takes as input platform
#####normalisation method that is used
#######project name
#####path to comp file - can get this from the master script i.e.miRNA_GO <- function(project,comp_file) {
library("biomaRt")
library("annotate")
library("org.Hs.eg.db")
library("GOstats")
library("xtable")
#wanted<-read.table("miRNA_toannotate.R",header=T,sep="\t",as.is=T)
miRNA_GO.R <- function(project,comp_file) {
comps <- read.table(comp_file, sep = "\t", as.is = T, header = F, 
        strip.white = T)
        for (i in 1:length(comps)) {
                print(comps[i])
        fg <- strsplit(comps[,i],"=")
        ci <- fg[[i]][1]
####need to do this for all comparisons
        #filteredx <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_annotated_filtered_DE.txt",sep=""),row.names = 1, header = T, sep = "\t") 
#filteredx <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_annotated_filtered_DE.txt",sep=""), header = T, sep = "\t",as.is=T) 
filteredx <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_all.txt",sep=""), header = T, sep = "\t",as.is=T) 



#sigGene=rownames(filteredx)
all_x_probeids <- gsub('_st',"",filteredx$ID)

#ensembl=useMart("ensembl")
#ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl=useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')

miRNA <- getBM( 
        attributes = c("mirbase_id", "entrezgene"), 
        filters = c("mirbase_id"), 
        values = all_x_probeids,
        #values = list(TRUE), 
        mart = ensembl
);
 

#wanted_annot <- allx[allx$mirbase_id %in% all_x_probeids, ]
	#filtered



#sigLL <- miRNA$entrezgene
sigLL <- miRNA$entrezgene
#if (length(sigLL) > 1) {

#sigLL <- unique(unlist(lookUp(sigGene,'lumiHumanAll.db','ENTREZID')))
sigLL <- as.character(sigLL[!is.na(sigLL)]) 
goTypes=c("BP","CC","MF")
for(go in (goTypes))
{
	print(go)
	params <- new("GOHyperGParams",geneIds= sigLL, annotation="org.Hs.eg", ontology=go, pvalueCutoff= 0.005, conditional=FALSE, testDirection="over")

	hgOver <- hyperGTest(params)
	df <- summary(hgOver)
	row <- nrow(df)
	print (row)
	if (row == 0) next
	for (i in 1:length(df[,7])){
		df[i,7] <- paste('<a href="http://amigo.geneontology.org/amigo/term/',df[i,1],'" target="blank">',df[i,7],'</a>',sep='')
		}
	print(xtable(df, caption=paste("Gene to GO ",go," test for over representation"),align="llrrrrrl",digits=3,display=c("s","s","f","f","f","d","d","s")),type="html",file=paste("ominer_results/",project,"/GO/",ci,go,"_over.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)   
	## Get the p-values of the test 
	gGhyp.pv <- pvalues(hgOver)
	## Adjust p-values for multiple test (FDR) 
	gGhyp.fdr <- p.adjust(gGhyp.pv, 'fdr')
	## select the Go terms with adjusted p-value less than 0.01 
	sigGO.ID <- names(gGhyp.fdr[gGhyp.fdr < 0.005]) ####p is the cutoff threshold for finding significant GO terms associated with significant DE genes
	## Here only show the significant GO terms of BP (Molecular Function) ##For other categories, just follow the same procedure. 
	if (length(sigGO.ID > 1)) {
	
		sigGO.Term <- getGOTerm(sigGO.ID)[[go]]
		gg.counts <- geneCounts(hgOver)[sigGO.ID]
		total.counts <- universeCounts(hgOver)[sigGO.ID]
		 
		ggt <- unlist(sigGO.Term)
		numCh <- nchar(ggt)
		ggt2 <- substr(ggt, 1, 17)
		ggt3 <- paste(ggt2, ifelse(numCh > 17, "...", ""), sep="")
			
		## 	## output the significant GO categories as a table
		ggMat <- matrix(c(names(sigGO.Term), ggt3, signif(gGhyp.pv[sigGO.ID],5), gg.counts, total.counts),
				byrow=FALSE, nc=5, dimnames=list(1:length(sigGO.Term), c("GO ID",
				"Term", "p-value","Significant Genes No.", "Total Genes No.")))
		
		write.table(ggMat, paste("ominer_results/",project,"/","GO/",ci,go,"over_results.txt", 
			   sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
			   
		#file.create(paste("ominer_results/",project,"/GO/",ci,go,"_over.html",sep=""))
		#htmlReport(hgOver,file=paste("ominer_results/",project,"/GO/",ci,go,"_over.html",sep=""))
	} 
}


for(go in (goTypes))
{
	print(go)
	params <- new("GOHyperGParams",geneIds= sigLL, annotation="org.Hs.eg", ontology=go, pvalueCutoff= 0.005, conditional=FALSE, testDirection="under")
	
	hgOver <- hyperGTest(params)
	df <- summary(hgOver)
	row <- nrow(df)
	print (row)
	if (row == 0) next
	for (i in 1:length(df[,7])){
		df[i,7] <- paste('<a href="http://amigo.geneontology.org/amigo/term/',df[i,1],'" target="blank">',df[i,7],'</a>',sep='')
		}
	print(xtable(df, caption=paste("Gene to GO ",go," test for under representation"),align="llrrrrrl",digits=3,display=c("s","s","f","f","f","d","d","s")),type="html",file=paste("ominer_results/",project,"/GO/",ci,go,"_under.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)   

	## Get the p-values of the test 
	gGhyp.pv <- pvalues(hgOver)
	## Adjust p-values for multiple test (FDR) 
	gGhyp.fdr <- p.adjust(gGhyp.pv, 'fdr')
	## select the Go terms with adjusted p-value less than 0.01 
	sigGO.ID <- names(gGhyp.fdr[gGhyp.fdr < 0.005]) ####p is the cutoff threshold for finding significant GO terms associated with significant DE genes
	if (length(sigGO.ID > 1)) {
	
		## Here only show the significant GO terms of BP (Molecular Function) ##For other categories, just follow the same procedure. 
		sigGO.Term <- getGOTerm(sigGO.ID)[[go]]
		gg.counts <- geneCounts(hgOver)[sigGO.ID]
		total.counts <- universeCounts(hgOver)[sigGO.ID]
		ggt <- unlist(sigGO.Term)
		numCh <- nchar(ggt)
		ggt2 <- substr(ggt, 1, 17)
		ggt3 <- paste(ggt2, ifelse(numCh > 17, "...", ""), sep="")
		
		## 	## output the significant GO categories as a table
		ggMat <- matrix(c(names(sigGO.Term), ggt3, signif(gGhyp.pv[sigGO.ID],5), gg.counts, total.counts),
			byrow=FALSE, nc=5, dimnames=list(1:length(sigGO.Term), c("GO ID",
			"Term", "p-value","Significant Genes No.", "Total Genes No.")))
	
		write.table(ggMat, paste("ominer_results/",project,"/","GO/",ci,go,"under_results.txt", 
		   sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
		   
		#file.create(paste("ominer_results/",project,"/GO/",ci,go,"_under.html",sep=""))
		#hgOver
		#htmlReport(hgOver,file=paste("ominer_results/",project,"/GO/",ci,go,"_under.html",sep=""))
	}
}
}

}
}

