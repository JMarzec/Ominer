##########################################################################################
#
#File name: miRNA_limma_annotation.R
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
#Description: normalisation of raw data from Affymetrix CEL files
#Arguments are:dataType this refers to the type of data input i.e CEL for raw data files and normalised for the normalised expression matrix
#targets = path to the target file 
#analysis = paired/unpaired for a (1) paired or an (2) unpaired analysis
#data = "filtered_data.txt" this is the output of the function filtering
#comp = full path to the comaprisons file
#replicates = yes/no yes if the dataset contains repliactes and no if there are no replicates
#limmamethod = method used within limma to compare between different groups e.g. separate
#adjust = method to adjust the false discovery rate (FDR) e.g BH, BY etc.
#pvalue = pvalue threshold used to report genes that are differentially expressed and meet the cutoff threshold(s)
#foldchange = log2 fold change value used to report genes that meet the cutoff threshold(s)
#project - given name for analysis
##########################################################################################
run_limma <- function (target, analysis, data, comp, replicates, limmamethod, 
    adjust, pvalue, foldchange, project) 
{
    pd <- read.AnnotatedDataFrame(target, header = T, row.name = "Name", 
        sep = "\t")
        #library("biomaRt")
       # ensembl <- useMart("ensembl");
#ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl);

# attributes available for download
# listAttributes(ensembl);

# retrieve user-specified attributes from BioMart db, in this instance miRNA annotations


        readdir <- paste("ominer_results",project,sep="/")
    Normdata <- read.table(paste(readdir,"/norm/",data,sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1)
    A.data <- new("ExpressionSet", exprs = data.matrix(Normdata))
    comparisons <- read.table(comp, sep = "\t", as.is = T, header = F, 
        strip.white = T)
    lev <- unique(pd$Target)
    lev
    f <- factor(pd$Target, levels = lev)
    if (analysis == "unpaired") {
        design <- model.matrix(~0 + f)
        colnames(design) <- make.names(lev)
    }
    else {
        levp = unique(pd$Pairs)
        p = factor(pd$Pairs, levels = levp)
        design <- model.matrix(~-1 + f + p)
        colnames(design) <- sub("f", "", colnames(design))
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
        })
        fit <- lmFit(A.data, design, block = block, correlation = dupcor$consensus)
    }
    else {
        fit <- lmFit(A.data, design)
    }
    myContrasts <- paste(c(comparisons$V1), collapse = ",")
    prestr = "makeContrasts("
    poststr = ",levels=design)"
    commandstr = paste(prestr, myContrasts, poststr, sep = "")
    cont.dif <- eval(parse(text = commandstr))
    fit2 <- contrasts.fit(fit, cont.dif)
    ebA <- eBayes(fit2)
    ebA.decideTests <- decideTests(ebA, method = limmamethod, 
        adjust.method = adjust, p.value = pvalue, lfc = foldchange)
    write.table(ebA.decideTests, paste("ominer_results/",project,"/","DifferentialExpression","/ebA_allresults.txt", sep=""),
        sep = "\t", row.names = TRUE, quote = FALSE)
    dim(ebA.decideTests@.Data)
    summary(ebA.decideTests)
    sum <- summary(ebA.decideTests)
    write.table(sum, paste("ominer_results/",project, "/","DifferentialExpression","/decideTestsSummary.txt",sep=""), sep = "\t", 
        row.names = TRUE, quote = FALSE)
    write.table(ebA, paste("ominer_results/",project, "/","DifferentialExpression","/ebA.txt", sep=""),sep = "\t", row.names = TRUE, 
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
        x <- topTable(ebA, coef = comps[i], adjust.method = adjust, 
            sort.by = "logFC", number = nrow(A.data))
            pval <- as.numeric(pvalue)
            fc <- as.numeric(foldchange)
        filtered_x <- topTable(ebA, coef = comps[i], n = 1000, 
            genelist = fit$genes[, 1], adjust.method = adjust, 
            p.value = pval, lfc = fc, sort.by = "P")
        probeid = unique(x$ID)
        probeids <- gsub('_st',"",x$ID)
          write.table(x, paste("ominer_results/",project,"/","DifferentialExpression/",comps[i],"annotated.txt", 
            sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)


 #write.table(filtered_x, paste("ominer_results/",project,"/","DifferentialExpression/",comps[i],"_Filtered.txt", 
  #          sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)



        #ensembl <- useMart("ensembl");
#ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl);

# attributes available for download
# listAttributes(ensembl);

# retrieve user-specified attributes from BioMart db, in this instance miRNA annotations
#miRNA <- getBM( 
   #     attributes = c("mirbase_id", "band", "chromosome_name","entrezgene","hgnc_symbol"), 
    #    filters = c("mirbase_id"), 
    #    #values = probeid,
     #   values = all_x_probeids, 
    #    mart = ensembl
#);

         #all_x_probeids <- gsub('_st',"",probeid)
#all_annotated_list <- miRNA[miRNA$mirbase_id %in% all_x_probeids, ]


   #write.table(all_annot, paste("ominer_results/",project,"/","DifferentialExpression/",comps[i], "_","all_DE.txt", 
    #       sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)

if (length(filtered_x$ID) > 0) {
		probeids <- gsub('_st',"",filtered_x$ID)
		filtered_annot <- cbind(probeids,filtered_x)
		 write.table(filtered_annot, paste("ominer_results/",project,"/","DifferentialExpression/",comps[i], "_","Filtered.txt", 
            sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
		#miRNA_annotation <- readLines("miRNA_changed.txt")
		#miRNA_annotation <- read.table("miRNA.txt", sep = "\t", as.is = T, header = T, strip.white = TRUE)


		#miRNA_accession <- miRNA_annotation[miRNA_annotation$ID %in% filtered_x_probeids]
		#miRNA_probe_annot <- cbind(miRNA_accession,filtered_x_probeids)
		#miRNA_annot <- miRNA_annotation$Accession
		#mature2_ID <- miRNA$Mature2_ID
		#new_annot <- cbind(miRNA_annot,mature2_ID)
		#new_names = NULL
		#get_lines = NULL
		#for (i in 1:length(filtered_x_probeids)) {
		#prot_string <- filtered_x_probeids[i]
		#new_prot_string <- tolower(c(prot_string))
		#all_names <- c(new_prot_string)
		#print (prot_string)
		#print (new_prot_string)
			#new_names <- c(new_names,grep(new_prot_string,miRNA_annotation$ID))
			#new_names <- c(new_names,grep(wanted,miRNA_annotation$Accession))
			#get_lines <- miRNA_annotation[new_names,]
			#new_lines <- c(get_lines)
			
			
			
#all_mirbase_IDS <- cbind(new_lines$Accession,new_lines$ID)
		
#filtered_annotated_list <- miRNA[miRNA$mirbase_id %in% filtered_x_probeids, ]
#new_filtered_merged <- cbind(filtered_x_probeids,filtered_annotated_list)
#all_filtered_A <- merge(new_filtered_merged,filtered_x_probeids,by.x = "probeid",by.y="ID")
#write.table(all_mirbase_IDS, paste("ominer_results/",project,"/","DifferentialExpression/",comps[i], "_","probe_annot_DE_getlines.txt", 
  #          sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)



       #write.table(filtered_annot, paste("ominer_results/",project,"/","DifferentialExpression/",comps[i], "_","annotated_filtered_DE.txt", 
        #    sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)

      #all_annotated_list <- miRNA[miRNA$mirbase_id %in% probeid, ]
                             
    }
    }
    
}
