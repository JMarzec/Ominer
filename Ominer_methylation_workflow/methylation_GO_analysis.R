############################################################################################################
#
#File name methylation_GO_analysis.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#######DESCRIPTION:GOstats analysis for Illumina methylation 450K and 27K chips try doing teh same way as illumina_expression gene ontology analysis
#######Arguments: (1) platform - 27k or 450k, project - name given to project, comp_file - full path to comparisons file 
#######output - significant Gene Ontology (GO ) terms associated to dataset 



gostats_analysis <- function(platform,project,comp_file) {
if (platform == "450k"){
platform <- "IlluminaHumanMethylation450k"
}
library("annotate")
library(paste(platform,".db",sep=""),character.only = TRUE)
#library("hgu133plus2.db")
library("GOstats")
library("xtable")
####Filter out records that do not have any annotation
#norm_file <- paste(norm,".exp",sep="")
#A.data <- read.table(norm_file,row.names = 1, header = T, sep = "\t")
A.data <- read.table(paste("ominer_results/",project,"/","norm","/","normalised.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1)
	dblookup <- paste(platform,".db",sep="")
                entrezIds <- sapply(lookUp(rownames(A.data), dblookup, "ENTREZID"), function(x) x[1])
               entrezIds <- as.character(entrezIds[!is.na(entrezIds)])
	       
comps <- read.table(comp_file, sep = "\t", as.is = T, header = F, 
        strip.white = T)
for (i in 1:length(comps$V1)) {

                print(comps$V1[i])
        fg <- strsplit(comps[,1],"=")
        ci <- fg[[i]][1]
        print ("This is the comparison")
        print(ci)
	 filteredx <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_Filtered.txt",sep=""),row.names = 1, header = T, sep = "\t")
         #filteredx <- read.table(paste(ci,"_annotated.txt",sep=""),row.names = 1, header = T, sep = "\t")
probeid=rownames(filteredx)
#platform = "hgu133plus2"
####collect entrez IDs for all of the affy_probeids passing the FC and pval cutoffs


dblookup <- paste(platform,".db",sep="")
                filtered_entrezIds <- sapply(lookUp(rownames(filteredx), dblookup, "ENTREZID"), function(x) x[1])
               all_filtered_entrezIds <- as.character(filtered_entrezIds[!is.na(entrezIds)])

	       
	       
goTypes=c("BP","CC","MF")
for(go in (goTypes))
{
	print(go)
	params <- new("GOHyperGParams",geneIds=all_filtered_entrezIds,universeGeneIds=entrezIds,annotation=dblookup, ontology=go, pvalueCutoff= 0.005, conditional=FALSE, testDirection="over")
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

sigGO.ID <- names(gGhyp.fdr[gGhyp.fdr < 0.05]) ####p is the cutoff threshold for finding significant GO terms associated with significant DE genes
if (length(sigGO.ID > 1)) {
## Here only show the significant GO terms of BP (Molecular Function) ##For other categories, just follow the same procedure. 
sigGO.Term <- getGOTerm(sigGO.ID)[[go]]

#sigGO.Term.mf <- getGOTerm(sigGO.ID)[["MF"]]
#sigGO.Term.cc <- getGOTerm(sigGO.ID)[["CC"]]
##get gene counts at each GO category
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

#write.table(ggMat, paste("ominer_results/",project,"/","GO/",ci,go,"over_results.txt", 
          # sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
           
           #file.create(paste("ominer_results/",project,"/GO/",ci,go,"_over.html",sep=""))
                                #htmlReport(hgOver,file=paste("ominer_results/",project,"/GO/",ci,go,"_over.html",sep=""))

}
#file.create(paste("ominer_results/",project,"/GO/",ci,go,"_over.html",sep=""))
}
 goTypes=c("BP","CC","MF")
for(go in (goTypes))
{
	print(go)
	params <- new("GOHyperGParams",geneIds=all_filtered_entrezIds,universeGeneIds=entrezIds, annotation=dblookup, ontology=go, pvalueCutoff= 0.005, conditional=FALSE, testDirection="under")
        hgUnder <- hyperGTest(params)
	df <- summary(hgUnder)
	row <- nrow(df)
	print (row)
	if (row == 0) next
	for (i in 1:length(df[,7])){
		df[i,7] <- paste('<a href="http://amigo.geneontology.org/amigo/term/',df[i,1],'" target="blank">',df[i,7],'</a>',sep='')
		}
	print(xtable(df, caption=paste("Gene to GO ",go," test for under representation"),align="llrrrrrl",digits=3,display=c("s","s","f","f","f","d","d","s")),type="html",file=paste("ominer_results/",project,"/GO/",ci,go,"_under.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)   

## Get the p-values of the test 
gGhyp.pv <- pvalues(hgUnder)
## Adjust p-values for multiple test (FDR) 
gGhyp.fdr <- p.adjust(gGhyp.pv, 'fdr')
## select the Go terms with adjusted p-value less than 0.01 

sigGO.ID <- names(gGhyp.fdr[gGhyp.fdr < 0.05]) ####p is the cutoff threshold for finding significant GO terms associated with significant DE genes
if (length(sigGO.ID > 1)) {
## Here only show the significant GO terms of BP (Molecular Function) ##For other categories, just follow the same procedure. 
sigGO.Term <- getGOTerm(sigGO.ID)[[go]]

#sigGO.Term.mf <- getGOTerm(sigGO.ID)[["MF"]]
#sigGO.Term.cc <- getGOTerm(sigGO.ID)[["CC"]]
##get gene counts at each GO category
	gg.counts <- geneCounts(hgUnder)[sigGO.ID]
 	total.counts <- universeCounts(hgUnder)[sigGO.ID]
 
 	ggt <- unlist(sigGO.Term)
 	numCh <- nchar(ggt)
 	ggt2 <- substr(ggt, 1, 17)
 	ggt3 <- paste(ggt2, ifelse(numCh > 17, "...", ""), sep="")
 	
## 	## output the significant GO categories as a table
 	ggMat <- matrix(c(names(sigGO.Term), ggt3, signif(gGhyp.pv[sigGO.ID],5), gg.counts, total.counts),
     		byrow=FALSE, nc=5, dimnames=list(1:length(sigGO.Term), c("GO ID",
    		"Term", "p-value","Significant Genes No.", "Total Genes No.")))

#
write.table(ggMat, paste("ominer_results/",project,"/","GO/",ci,go,"under_results.txt", 
           sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
           
           #file.create(paste("ominer_results/",project,"/GO/",ci,go,"_under.html",sep=""))
                                #htmlReport(hgOver,file=paste("ominer_results/",project,"/GO/",ci,go,"_under.html",sep=""))

}
#file.create(paste("ominer_results/",project,"/GO/",ci,go,"_over.html",sep=""))
}




}
}
                                         