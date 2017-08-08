##########################################################################################
#
#File name: run_limma_illumina.R
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
#Arguments are:
#target = full path to the target file 
#analysis = paired/unpaired for a (1) paired or an (2) unpaired analysis
#data = "filtered_data.txt" this is the output of the function filtering
#comp = full path to the comaprisons file
#replicates = yes/no yes if the dataset contains repliactes and no if there are no replicates
#limmamethod = method used within limma to compare between different groups e.g. separate
#adjust = method to adjust the false discovery rate (FDR) e.g BH, BY etc.
#pvalue = pvalue threshold used to report genes that are differentially expressed and meet the cutoff threshold(s)
#foldchange = log2 fold change value used to report genes that meet the cutoff threshold(s)
#platform = platform of array data that is under analysis e.g. ht12v3 is for illuminaHumanv3 platform
#project - given name for analysis
##########################################################################################

run_limma <- function (target, analysis, data, comp, replicates, limmamethod, 
    adjust, pvalue, foldchange, platform,project) 
    
{
	library("limma")
	library("affy")
	library("illuminaHumanv2.db")
	library("annotate")
	library("lumi")
    pd <- read.AnnotatedDataFrame(target, header = T, row.name = "Name", 
        sep = "\t")
          readdir <- paste("ominer_results",project,sep="/")
    Normdata <- read.table(paste(readdir,"/","norm/",data,sep=""), sep = "\t", as.is = T, header = T, 
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
    write.table(ebA.decideTests, paste("ominer_results/",project,"/DifferentialExpression","/ebA_allresults.txt", sep=""),
        sep = "\t", row.names = TRUE, quote = FALSE)
        dim(ebA.decideTests@.Data)
    summary(ebA.decideTests)
    sum <- summary(ebA.decideTests)
     write.table(sum, paste("ominer_results/",project, "/DifferentialExpression","/decideTestsSummary.txt",sep=""), sep = "\t", 
        row.names = TRUE, quote = FALSE)
    write.table(ebA, paste("ominer_results/",project, "/DifferentialExpression","/ebA.txt", sep=""),sep = "\t", row.names = TRUE, 
        quote = FALSE)
        comps = colnames(ebA.decideTests)
    for (i in 1:length(comps)) {
    	 pdf(paste("ominer_results/",project,"/DifferentialExpression","/hist", comps[i], ".pdf", sep = ""))
       
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
        filtered_x <- topTable(ebA, coef = comps[i], n = 1000, 
            genelist = fit$genes[, 1], adjust.method = adjust, 
            p.value = pvalue, lfc = foldchange, sort.by = "P")
        probeid = unique(x$ID)
        filtered_id = unique(filtered_x$ID)
       write.table(x, paste("ominer_results/",project,"/","DifferentialExpression/",comps[i],"all.txt", 
         sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
        library("biomaRt")
#wanted<-read.table("illumina_probes_annot.txt",header=T,sep="\t",as.is=T)
#ensembl=useMart("ensembl")
 mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
if (platform == "ht12v3") {
	platform_annot = "illumina_humanht_12_v3"
	ensembl = useDataset("hsapiens_gene_ensembl",mart=mart)
	theAttributes = c(platform_annot,"refseq_mrna_predicted","chromosome_name","band","hgnc_symbol", "description")
}
if (platform == "ht12v4") {
	platform_annot = "illumina_humanht_12_v4"
	ensembl = useDataset("hsapiens_gene_ensembl",mart=mart)
	theAttributes = c(platform_annot,"refseq_mrna_predicted","chromosome_name","band","hgnc_symbol", "description")
}
if (platform == "MO_8_v2") {
	platform_annot = "illumina_mouseref_8_v2"
	ensembl = useDataset("mmusculus_gene_ensembl",mart=mart)
	plat_filt <- "illumina_mouseref_8_v2"
	theAttributes = c(plat_filt,"refseq_mrna_predicted","chromosome_name","band","hgnc_symbol", "description")
	
}

theFilters = c(platform_annot)   ######or illumina_humanht_v3

wanted_annot <- getBM(attributes = theAttributes, filters = theFilters, values = x$ID, mart=mart)
colnames(wanted_annot) <- c("probe_id","AGI","chromosome","band","symbol","gene_name")
annot_all <- merge(wanted_annot,x, by.x="probe_id",by.y="ID",all=T)
write.table(annot_all, paste("ominer_results/",project,"/","DifferentialExpression/",comps[i],"annotated.txt", 
           sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
if (length(filtered_id) > 0) { 
filtered_annot <- getBM(attributes = theAttributes, filters = theFilters, values = filtered_x$ID, mart=mart)
colnames(filtered_annot) <- c("probe_id","AGI","chromosome","band","symbol","gene_name")
filtered_all <- merge(filtered_annot,filtered_x, by.x="probe_id",by.y="ID",all=T)
write.table(filtered_all, paste("ominer_results/",project,"/","DifferentialExpression/",comps[i],"_Filtered.txt", 
            sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
#write.table(filtered_x, paste("ominer_results/",project,"/","DifferentialExpression/",comps[i],"FilteredX.txt", 
            #sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
}
 # write.table(annot_all,"all_annotated.txt",sep="\t",quote=FALSE,row.names=FALSE)

        
        
            
        
        #library("lumiHumanAll.db")
        #if (require(lumiHumanAll.db) & require(annotate)) {
        #geneSymbol <- getSYMBOL(probeid, 'illuminaHumanv2.db')
        #geneName <- sapply(lookUp(probeid, 'lumiHumanAll.db', 'GENENAME'),function)
       #annot_all <- data.frame(ID= probeid, geneSymbol=geneSymbol)
       }
 
           # }
  }
