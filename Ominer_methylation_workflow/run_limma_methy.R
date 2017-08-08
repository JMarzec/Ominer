############################################################################################################
#
#File name run_limma_methy.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
########Arguments - 1. limma_target_file = target_file, analysis - paired/unpaired analysis, data = normalised data file, comp_file = pathe to comparison file, replicates = yes/no, limammethod = separate etc, adjust = BH, BF etc,
########foldchange = 2.0, platform = 450k or 250k, project name of project
########outputs annotated lists of statistically signifiacntly methylated probes 
library("affy")
library("simpleaffy")
library("affyPLM")
library("limma")
library("annotate")


run_limma <- function (limma_target_file, analysis, data, comp_file, replicates, limmamethod, 
    adjust, pvalue, foldchange, platform,project) 
{
    pd <- read.table(limma_target_file, header = T,sep = "\t")
          readdir <- paste("ominer_results",project,sep="/")
    Normdata <- read.table(paste(readdir,"/norm/",data,sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1)
   
    A.data <- new("ExpressionSet", exprs = data.matrix(Normdata))
    comparisons <- read.table(comp_file, sep = "\t", as.is = T, header = F, 
        strip.white = T)
    lev <- unique(pd$Name)
    lev
    f <- factor(pd$Name, levels = lev)
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
             ID <- rownames(x)
            all_x <- cbind(ID,x)
        probeList = unique(all_x$ID)
        pl <- as.character(probeList)
        #library("IlluminaHumanMethylation27k.db")
           #library("IlluminaHumanMethylation450k.db")

        #library("annotate")
        #geneSymbol <- getSYMBOL(probeList, "IlluminaHumanMethylation450k.db")
        if (platform == "450k"){
        	annotation_database <- "IlluminaHumanMethylation450k.db"
        library("IlluminaHumanMethylation450k.db")
        }
        if (platform == "250k") {
        	annotation_database <- "IlluminaHumanMethylation250k.db"
                library("IlluminaHumanMethylation250k.db")
        }
        #library(annotation_database)
        geneSymbol <- sapply(lookUp(pl, annotation_database, 
            "SYMBOL"), function(x) x[1])
        geneName <- sapply(lookUp(pl, annotation_database, 
            "GENENAME"), function(x) x[1])
            chromosome <- sapply(lookUp(pl, annotation_database, 
            "CHR"), function(x) x[1])
            map <- sapply(lookUp(pl, annotation_database, 
            "MAP"), function(x) x[1])
            dmr <- sapply(lookUp(pl, annotation_database, 
            "DMR"), function(x) x[1])
            cpgi_location <-  sapply(lookUp(pl, annotation_database, 
            "CPGILOCATION"), function(x) x[1])
            #print ("GENESYMBOL")
            #print (geneSymbol)
            #print ("genename")
            #print (geneName)
            #print ("chromosome")
            #print (chromosome)
            #print ("map")
            #print (map)
            #print ("dmr")
            #print (dmr)
            #print ("cpgi_location")
            #print (cpgi_location)
            #print ("these are the probes")
            #print (pl)

        fit$genes <- data.frame(ID = pl, geneSymbol = geneSymbol, 
            geneName = geneName,map=map,dmr=dmr,cpgi_location = cpgi_location, stringsAsFactors = FALSE)
           # print ("This is fit$genes content")
            #print (fit$genes)
        ax <- topTable(ebA, coef = comps[i], adjust.method = adjust, 
            sort.by = "logFC", genelist = fit$genes, number = nrow(A.data))
        write.table(x, paste("ominer_results/", project, "/DifferentialExpression/",comps[i], "TopTable.txt", 
            sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)
             #write.table(sum, paste("ominer_results/",project,"/","DifferentialExpression","/decideTestsSummary.txt",sep=""),sep = "\t", 
        #row.names = TRUE, quote = FALSE)
        print(ax)
       # write.table(ax,  paste("ominer_results/",project, "/", "DifferentialExpression/", comps[i], "annotated_TopTable.txt", 
            #sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)
 write.table(ax,  paste("ominer_results/",project, "/", "DifferentialExpression/", comps[i], "annotated.txt", 
            sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)

        filtered_x <- topTable(ebA, coef = comps[i],  adjust.method = adjust, 
            p.value = pvalue, lfc = foldchange, sort.by = "P")
            if (length(filtered_x) > 0) {
      filtered_ids <- rownames(filtered_x)
      all_ids_filtered <- cbind(filtered_ids,filtered_x)
      filtered_probeList <- unique(all_ids_filtered)
      filtered_pl <- as.character(filtered_probeList)
      
        filtered_geneSymbol <- sapply(lookUp(filtered_pl, annotation_database, 
            "SYMBOL"), function(x) x[1])
        filtered_geneName <- sapply(lookUp(filtered_pl, annotation_database, 
            "GENENAME"), function(x) x[1])
            filtered_chromosome <- sapply(lookUp(filtered_pl, annotation_database, 
            "CHR"), function(x) x[1])
            filtered_map <- sapply(lookUp(filtered_pl, annotation_database, 
            "MAP"), function(x) x[1])
            filtered_dmr <- sapply(lookUp(filtered_pl, annotation_database, 
            "DMR"), function(x) x[1])
            filtered_cpgi_location <-  sapply(lookUp(filtered_pl, annotation_database, 
            "CPGILOCATION"), function(x) x[1])
            

        #fit$genes <- data.frame(ID = filtered_pl, geneSymbol = filtered_geneSymbol, 
           # geneName = filtered_geneName, chromosome=filtered_chromosome,map=filtered_map,dmr=filtered_dmr,cpgi_location = filtered_cpgi_location, stringsAsFactors = FALSE)
        filtered_ax <- topTable(ebA, coef = comps[i], adjust.method = adjust, 
            sort.by = "logFC", genelist = fit$genes, number = nrow(A.data))
        #write.table(filtered_x, paste("ominer_results/", project, "/DifferentialExpression/",comps[i], "filtereannotated.txt", sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)
        write.table(filtered_ax,paste("ominer_results/",project, "/DifferentialExpression/",comps[i], "_Filtered.txt",sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)

     }
     
        }
    }
    

