##########################################################################################
#
#File name: run_limma_affy.R
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
#platform = platform of array data that is under analysis e.g. hgu133plus2 is for Affymetrix hgu133plus2 array platform
#project - given name for analysis
##########################################################################################
library("affy")
library("simpleaffy")
library("affyPLM")
library("limma")
library("annotate")


run_limma <- function (target, analysis, data, comp, replicates, limmamethod, 
    adjust, pvalue, foldchange, platform,project) 
{
    pd <- read.AnnotatedDataFrame(target, header = T, row.name = "Name", 
        sep = "\t")
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
        x <- topTable(ebA, coef = comps[i], adjust.method = adjust, 
            sort.by = "logFC", number = nrow(A.data))
            write.table(x,paste("ominer_results/",project,"/","DifferentialExpression","/",comps[i], "un_annotated.txt", 
            sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)
            
        filtered_x <- topTable(ebA, coef = comps[i],
             adjust.method = adjust, 
            p.value = pv, lfc = lfval, sort.by = "P")
            probeid <- rownames(x)
        AGI <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
            "ACCNUM", sep = ""))), function(AGI) {
            return(paste(AGI, collapse = "; "))
        })))
        Chr <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
            "CHR", sep = ""))), function(Chr) {
            return(paste(Chr, collapse = "; "))
        })))
        symbol <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
            "SYMBOL", sep = ""))), function(symbol) {
            return(paste(symbol, collapse = "; "))
        })))
        NAME <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
            "GENENAME", sep = ""))), function(name) {
            return(paste(name, collapse = "; "))
        })))
        if (platform != "mouse4302") {
        MAP <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
            "MAP", sep = ""))), function(MAP) {
            return(paste(MAP, collapse = "; "))
        })))
        }
        if (platform == "mouse4302") {
         MAP <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
            "CHRLOC", sep = ""))), function(MAP) {
            return(paste(MAP, collapse = "; "))
        })))
        
        }
        AnnotA <- data.frame(probeid, AGI, Chr, MAP, symbol, 
            NAME, row.names = NULL)
            ID <- rownames(x)
            dx <- cbind(ID,x)
           x <- dx
           print ("I got to before merging here 1")
        allA <- merge(AnnotA, x, by.x = "probeid", by.y = "ID", 
            all = T)
            print ("I got to after here merging 2!")
        write.table(allA,paste("ominer_results/",project,"/","DifferentialExpression","/",comps[i], "annotated.txt", 
            sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
            ID = rownames(filtered_x)
            filt_x <- cbind(ID,filtered_x)
            number_filtered <- length(filt_x)
        fpid = unique(filt_x$x)
                write.table(filtered_x, paste("ominer_results/",project,"/","DifferentialExpression","/",comps[i],"_Filtered.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
         fp_length <- length(fpid)
        if (fp_length > 0) {
        print ("YES")
        AGI <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
           "ACCNUM", sep = ""))), function(AGI) {
            return(paste(AGI, collapse = "; "))
        })))
        Chr <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
         "CHR", sep = ""))), function(Chr) {
            return(paste(Chr, collapse = "; "))
        })))
        symbol <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
         "SYMBOL", sep = ""))), function(symbol) {
           return(paste(symbol, collapse = "; "))
        })))
        NAME <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
           "GENENAME", sep = ""))), function(name) {
           return(paste(name, collapse = "; "))
        })))
        if (platform != "mouse4302") {
        MAP <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
            "MAP", sep = ""))), function(MAP) {
            return(paste(MAP, collapse = "; "))
        })))
        }
        if (platform == "mouse4302") {
         MAP <- as.character(unlist(lapply(mget(probeid, env = get(paste(platform, 
            "CHRLOC", sep = ""))), function(MAP) {
            return(paste(MAP, collapse = "; "))
        })))
        }
        AnnotA <- data.frame(probeid, AGI, Chr, MAP, symbol, 
            NAME, row.names = NULL)
        allFiltA <- merge(AnnotA, filtered_x, by.x = "probeid", 
            by.y = "ID")
        write.table(allFiltA, paste("ominer_results/",project,"/","DifferentialExpression","/",comps[i],"_Filtered.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
    }
    }
    }
    

