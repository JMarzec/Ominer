##########################################################################################
#
#File name:run_limma_affy_RNAseq_modS.R
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
#####target - a simple analysis will have a target file consisting of three columns filename, samplename and target
#####analysis - whether a paired or an unpaired analysis is requested inputs to this argument are:paired/unpaired
#####data - filtered_data.txt - this is the filtered normalised expression matrix
#####comp - this is the .txt file containing information regarding the comparisons to be made for each dataset analysed
#####replicates - yes/no to whether the datasets contains technical replicates
#####limmamethod - method to be used within limma to compare biological groups
#####adjust - method to adjust the P-value/FDR e.g. BH - Benjamini-Hochberg, or BY
#####pvalue - this is the adjusted p value cutoff used when the list of differentially expressed genes is prepared/displayed
#####foldchange - this is the log2 fold change cutoff used when the list of differentially expressed genes is prepared/displayed
#####project - name of analysis that is run 
#####dataType - whether data is normalised or unnormalised

library("limma")
library("edgeR")
library("affy")

run_limma <- function (target, analysis, data, comp, replicates, limmamethod, 
    adjust, pvalue, foldchange,project,dataType) {
    pd <- read.AnnotatedDataFrame(target, header = T, row.name = "Name", 
        sep = "\t")
          readdir <- paste("ominer_results",project,sep="/")
    Normdata <- read.table(paste(readdir,"/norm/",data,sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1)
   ad = read.delim(paste(readdir,"/norm/",data,sep=""),header=TRUE,row.names=1)
A.data <- DGEList(counts=ad,genes=rownames(ad))
A.data <- calcNormFactors(A.data)
if (dataType == "normalised") {
    A.data <- new("ExpressionSet", exprs = data.matrix(Normdata))
    }
    
    comparisons <- read.table(comp, sep = "\t", as.is = T, header = F, 
        strip.white = T)
    lev <- unique(pd$Target)
    lev
    f <- factor(pd$Target, levels = lev)
    if (analysis == "unpaired") {
        design <- model.matrix(~0 + f)
        colnames(design) <- sub("f", "", colnames(design))
        
         if (dataType == "unnormalised") {
        v <- voom(A.data,design,plot=TRUE)
        colnames(design) <- make.names(lev)
     }
     
     
    }else {
        levp = unique(pd$Pairs)
        p = factor(pd$Pairs, levels = levp)
        design <- model.matrix(~-1 + f + p)
        if (dataType == "unnormalised") {
         v <- voom(A.data,design,plot=TRUE)
        
         }
          colnames(design) <- sub("f", "", colnames(design))
        
    }
    if (replicates == "yes") {
        block <- pd$Replicates
        dupcor <- duplicateCorrelation(A.data, design = design, 
            block = block)
               v <- voom(A.data,design,plot=TRUE)
        tryCatch(fit <- lmFit(v, design, block = block, 
            correlation = dupcor$consensus), error = function(err) {
            writeLines(err$message, fileErr)
             cat("Error with lmfit of replicates", file = errorfile, 
                sep = "\n", append = TRUE)
        })
        if (dataType == "unnormalised") {
        fit <- lmFit(v, design, block = block, correlation = dupcor$consensus)
        }
        if (dataType == "normalised") {
        fit <- lmFit(A.data, design, block = block, correlation = dupcor$consensus)
        }
    }
    if (dataType == "unnormalised") {
        fit <- lmFit(v, design)
        }
        if (dataType == "normalised") {
        fit <- lmFit(A.data,design)
        }
    
    myContrasts <- paste(c(comparisons$V1), collapse = ",")
    prestr = "makeContrasts("
    poststr = ",levels=design)"
    commandstr = paste(prestr, myContrasts, poststr, sep = "")
   cont.dif <- eval(parse(text = commandstr))
    fit2 <- contrasts.fit(fit, cont.dif)
    ebA <- eBayes(fit2)
    pv <- as.numeric(pvalue)
    lfval <- as.numeric(foldchange)
    ebA.decideTests <- decideTests(ebA, method = limmamethod, 
      adjust.method = adjust, p.value = pv, lfc = lfval)
    write.table(ebA.decideTests, paste("ominer_results/",project,"/","DifferentialExpression","/ebA_allresults.txt",sep=""), 
   sep = "\t", row.names = TRUE, quote = FALSE)
     
   dim(ebA.decideTests@.Data)
    summary(ebA.decideTests)
    sum <- summary(ebA.decideTests)
    write.table(sum, paste("ominer_results/",project,"/","DifferentialExpression","/decideTestsSummary.txt",sep=""),sep = "\t", 
      row.names = TRUE, quote = FALSE)
    write.table(ebA, paste("ominer_results/",project,"/","DifferentialExpression","/ebA.txt",sep=""), sep = "\t", row.names = TRUE, 
        quote = FALSE)
    comps = colnames(ebA.decideTests)
    for (i in 1:length(comps)) {
      pdf(paste("ominer_results/",project,"/","DifferentialExpression","/hist", comps[i], ".pdf", sep = ""))
      hist(ebA$p.value[, i], breaks = 100, col = "orange", 
            main = paste("Histogram for ", comps[i], sep = ""), 
           xlab = "pvalue")
       dev.off()
    }
    comps = colnames(ebA.decideTests)
 for (i in 1:length(comps)) {
        print(comps[i])
        if (dataType == "unnormalised") {
       x <- topTable(ebA, coef = comps[i], adjust.method = adjust, 
            sort.by = "logFC", number = nrow(v))
          }
          if (dataType == "normalised") {
          x <- topTable(ebA, coef = comps[i], adjust.method = adjust, 
            sort.by = "logFC", number = nrow(A.data))
          }
            
            #####annotate Ensembl gene IDS
            
            ##### Create 'not in' operator
write.table(x, paste("ominer_results/",project,"/","DifferentialExpression","/",comps[i],"_unannotated.txt",sep=""), sep = "\t", row.names = FALSE, 
        quote = FALSE)

##### Retrieve gene annotation information input genes - list 

     #Test to see whether user has enterd counts matrix with ENSG IDS as rownames
     
     tester <- grepl("^ENSG",rownames(Normdata))
     if (tester[1] == "TRUE") {
     
     
      #####annotate Ensembl gene IDS
     
    library("biomaRt")
    mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
    #mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
    
    theFilters = c("ensembl_gene_id")
 theAttributes =c("ensembl_gene_id","chromosome_name", "band","hgnc_symbol")
ampGenes <- getBM(attributes = theAttributes, filters = theFilters, values = x$ID,mart=mart)

    
    
    
allA <- merge(x,ampGenes,by.x="ID",by.y="ensembl_gene_id",all =T)

write.table(allA, paste("ominer_results/",project,"/","DifferentialExpression","/",comps[i],"annotated.txt",sep=""), sep = "\t", row.names = FALSE, 
        quote = FALSE)
                       
        filtered_x <- topTable(ebA, coef = comps[i],  adjust.method = adjust, 
            p.value = pv, lfc = lfval, sort.by = "P")
            
           

            number_filtered <- length(filtered_x)
if (number_filtered > 0) {
 #filtered_x_annot <- getGene( id=filtered_x$genes, type="ensembl_gene_id", mart = mart)
 
 filtered_x_annot <- getBM(attributes = theAttributes, filters = theFilters, values = filtered_x$ID,mart=mart)
            #new_c <- paste(allA$chromosome_name,allA$band,sep="")
            #new_table <- cbind(allA$genes,new_c,allA$hgnc_symbol,allA$logFC,allA$adj.P.Val)
           
   filtered_x_annot <- merge(filtered_x,filtered_x_annot,by.x="ID",by.y="ensembl_gene_id",all=T)
#colnames(new_table) <- c("ensembl_gene_id","locus","hgnc_symbol","logFC","P.Value")

             csome_band <- paste(filtered_x_annot$chromosome,filtered_x_annot$band,sep="")
new_filtered_A <- cbind(filtered_x$ID,csome_band,filtered_x_annot$hgnc_symbol,filtered_x_annot$logFC,filtered_x_annot$adj.P.Val,filtered_x_annot$description,filtered_x_annot$logFC,filtered_x_annot$logFC,filtered_x_annot$P.Value,filtered_x_annot$B,filtered_x_annot$adj.P.Val)
colnames(new_filtered_A) <- c("probeid","MAP","symbol","B","adj.P.Val","logFC","P.Value","B_values","strand","adj.P.Val")
 write.table(new_filtered_A, paste("ominer_results/",project,"/","DifferentialExpression","/",comps[i],"_Filtered.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE,col.names=TRUE)  ###this should be filtered_A problem with filtered_A so

}
             
            
                             ######Need to annotate if gene symbol is used as rowname(s)
        }
        else {
        library("biomaRt")
     mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
    #mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
   
    
	
	theAttributes = c("hgnc_symbol","ensembl_gene_id","band","chromosome_name")



theFilters = c("hgnc_symbol")   ######or illumina_humanht_v3

wanted_annot <- getBM(attributes = theAttributes, filters = theFilters, values = x$ID, mart=mart)
        allA_symbol <- merge(x,wanted_annot,by.x="ID",by.y="hgnc_symbol")
         all_datasym <- paste(allA_symbol$chromosome_name,allA_symbol$band,sep="")
    format_AllA <- cbind(allA_symbol$ensembl_gene_id,all_datasym,allA_symbol$ID,allA_symbol$logFC,allA_symbol$adj.P.Val)
    
    for_GO <- cbind(allA_symbol$ensembl_gene_id,all_datasym,allA_symbol$ID,allA_symbol$logFC,allA_symbol$P.Value)
        colnames(format_AllA) <-c("probe_id","locus","symbol","logFC","adjusted_p_value")
        colnames(for_GO) <-c("probe_id","locus","symbol","logFC","p_value")
         write.table(format_AllA, paste("ominer_results/",project,"/","DifferentialExpression","/",comps[i],"annotated.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE,col.names=TRUE)  ###this should be filtered_A problem with filtered_A so
 write.table(for_GO, paste("ominer_results/",project,"/","DifferentialExpression","/",comps[i],"_all_GO.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE,col.names=TRUE)  ###this should be filtered_A problem with filtered_A so



                
        filtered_x <- topTable(ebA, coef = comps[i],  adjust.method = adjust, 
            p.value = pvalue, lfc = foldchange, sort.by = "P")
            
           

            number_filtered <- length(filtered_x)
if (number_filtered > 0) {
filtered_wanted_annot <- getBM(attributes = theAttributes, filters = theFilters, values = filtered_x$ID, mart=mart)
        allA_filtered <- merge(filtered_x,filtered_wanted_annot,by.x="ID",by.y="hgnc_symbol")
        all_cl <- paste(allA_filtered$chromosome_name,allA_filtered$band,sep="")

       
    filtered_AllA <- cbind(allA_filtered$ensembl_gene_id,all_cl,allA_filtered$ID,allA_filtered$logFC,allA_filtered$adj.P.Val)    
        colnames(filtered_AllA) <-c("probe_id","locus","symbol","logFC","adjusted_p_value")
         write.table(filtered_AllA, paste("ominer_results/",project,"/","DifferentialExpression","/",comps[i],"_Filtered.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE,col.names=TRUE)  ###this should be filtered_A problem with filtered_A so


        }
        
        
        }
    
    
    }

}
