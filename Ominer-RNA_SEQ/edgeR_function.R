##########################################################################################
#
#File name:edgeR_function.R
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
#####data - filtered normalised expression matrix
#####target - .txt file containing information about the samples to be analysed. i.e. a simple analysis will have a traget file containing three columns one for Filename, samplename and target (target column contains the biological/phenotypic information for each sample)
#####analysis - whether a paired or an unpaired analysis is run
#####comp - this is the .txt file containing information regarding the comparisons to be made for each dataset analysed
#####replicates - yes/no to whether the datasets contains technical replicates
#####adjust - method to adjust the P-value/FDR e.g. BH - Benjamini-Hochberg, or BY
#####pvalue - this is the adjusted p value cutoff used when the list of differentially expressed genes is prepared/displayed
#####project - name of analysis that is run 
#####This code performs differential expression analysis from the filtered expression matrix using edgeR the output are lists of differentially expressed genes for each of the comparisons 

run_edgeR <- function (data,target,analysis,comp,replicates,adjust,pvalue,project) {
	library("edgeR")
	readdir <- paste("ominer_results",project,sep="/")
	ReadsCountTable <- read.delim(paste(readdir,"/norm/",data,sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1)
	#ReadsCountTable = read.delim(data,header=TRUE,row.names=1)
	pd<-read.table(target,header=T,sep="\t")
	     #f <- factor(pd$Target, levels = lev)
	A.data <- DGEList (counts=ReadsCountTable,group=pd$Target)
        A.data <- calcNormFactors(A.data)
	ng <- nrow(A.data)
	 comparisons <- read.table(comp, sep = "\t", as.is = T, header = F, 
        strip.white = T)
	group <- factor(pd$Target)

if (analysis == "unpaired") {
        design <- model.matrix(~0 + group)
        colnames(design) <- levels(group)
    }
    else {
        #levp = unique(pd$Pairs)
        p = factor(pd$Pairs)
        design <- model.matrix( ~group + p)
        #colnames(design) <- sub("group", "", colnames(design))
        #colnames(design) <- levels(group)
    }
    if (replicates == "yes") {
        block <- pd$Replicates
        dupcor <- duplicateCorrelation(A.data, design = design, 
            block = block)
        tryCatch(fit <- lmFit(A.data, design, block = block, 
            correlation = dupcor$consensus), error = function(err) {
            writeLines(err$message, fileErr)
            cat("Error with lmfit of replicates", file = errorfile, 
                sep = "\n", append = TRUE)
                 A.data <- estimateGLMCommonDisp(A.data,design)
A.data <- estimateGLMTrendedDisp(A.data,design)
A.data <- estimateGLMTagwiseDisp(A.data,design)

        fit <- glmFit(A.data, design, block = block, correlation = dupcor$consensus)
        })
       }
    
    else {
    	 A.data <- estimateGLMCommonDisp(A.data,design)
A.data <- estimateGLMTrendedDisp(A.data,design)
A.data <- estimateGLMTagwiseDisp(A.data,design)
        fit <- glmFit(A.data, design)
    }
 myContrasts <- paste(c(comparisons$V1), collapse = ",")
    prestr = "makeContrasts("
    poststr = ",levels=design)"
    commandstr = paste(prestr, myContrasts, poststr, sep = "")
    cont.dif <- eval(parse(text = commandstr))
    lrt <- glmLRT(fit,contrast=cont.dif)
    
pval <- as.numeric(pvalue)    
de <- decideTestsDGE(lrt,adjust.method=adjust,p.val=pval)
dim(de@.Data)
summary(de)
sum <- summary(de)
comp <- colnames(de)   
    
    for (i in 1:length(comparisons)) {
    	print(comparisons[i])
    	    	fg <- strsplit(comparisons[,i],"=")
    	    	ci <- fg[[i]][1]
lrt <- glmLRT(fit,contrast=cont.dif[,i])
#for (i in 1:length(comps)) {
        pdf(paste("ominer_results/",project,"/","DifferentialExpression","/hist", ci, ".pdf", sep = ""))
        hist(lrt$AveLogCPM, breaks = 100, col = "orange", 
            main = paste("Histogram for ", ci, sep = ""), 
            xlab = "pvalue")
        dev.off()
    #}

de <- decideTestsDGE(lrt,adjust.method=adjust,p.val=pvalue)
comps <- colnames(de)
sum <- summary(de)
  write.table(de, paste("ominer_results/",project,"/","DifferentialExpression","/ebA_allresults.txt",sep=""), 
        sep = "\t", row.names = TRUE, quote = FALSE)
summary(de)
sum <- summary(de)
colnames(sum) <- ci
 write.table(sum, paste("ominer_results/",project,"/","DifferentialExpression","/decideTestsSummary.txt",sep=""),sep = "\t", 
        row.names = TRUE,col.names=TRUE, quote = FALSE)



filtered_all <- topTags(lrt,adjust.met=adjust,n=ng)
#ensembl_geneids <- rownames(filtered_all)
#all_table <- cbind(ensembl_geneids,filtered_all)
#colnames(all_table) <- c("ensg_id","logFC","logCPM","LR","PValue","FDR")

   comp_file_name <- strsplit(comparisons[,i],"=")
    
    ci <- comp_file_name[[i]][1]
  #library("biomaRt")
   #         up <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  #theFilters = c("ensembl_gene_id")
#theAttributes =c("ensembl_gene_id","chromosome_name", "start_position","end_position","hgnc_symbol")
#ampGenes <- getBM(attributes = theAttributes, filters = theFilters, values = all_table$ensg_id, mart=up)
#new <- cbind(ampGenes$chromosome_name,ampGenes$start_position,ampGenes$end_position)
#new_c <- paste("chr",ampGenes$chromosome_name,":",ampGenes$start_position,"-",ampGenes$end_position,sep="")
#together <- cbind(ampGenes$ensembl_gene_id,new_c,ampGenes$hgnc_symbol)
#colnames(together) <- c("ensembl_gene_id","locus","symbol")
#gene_filtered = merge(together,all_table,by.x="ensembl_gene_id",by.y="ensg_id",all=TRUE)

#unique_all <-filtered_all[!duplicated(filtered_all$ensembl_gene_id), ] 
	
	write.table(filtered_all, paste("ominer_results/",project,"/","DifferentialExpression/",ci,"_unannotated.txt",sep=""), sep = "\t", row.names = TRUE, 
        quote = FALSE)
	
	}
	
}