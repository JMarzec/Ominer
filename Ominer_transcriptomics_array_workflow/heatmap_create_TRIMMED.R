##########################################################################################
#
#File name: heatmap_create_TRIMMED.R
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
#Description: Generation of heatmaps from the list of differentially expressed genes that passed filters
#Arguments are:
#target_file is the full path to the target file
#project is the given project name
#comp_file is the full path to the comparisons file
##########################################################################################
heatmap_generate <- function (target, project,comp_file) 
{
    library("gplots")
    library("limma")
    library("hgu133plus2.db")
    readdir <- paste("ominer_results",project,sep="/")
    A.data <- read.table(paste(readdir,"/","norm","/filtered_data.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1)
        A.data <- as.matrix(A.data)
        number_rows <- nrow(A.data)
        print ("THIS IS THE NUMBER OF ROWS")
        print (number_rows)
        number_columns <- ncol(A.data)
        if (number_columns > 20) {
        nc <- 20
        trim_A.data <- subset(A.data,select = c(1:nc))
        print ("THIS IS THE TRIMMED MATRIX")
        #print (trim_A.data)
       }
       else {
       nc <- number_columns
       trim_A.data <- A.data
       }
    pd <- read.table(target, sep = "\t", as.is = TRUE, strip.white = TRUE, 
        header = T)
        trimmed_target <- pd$Target[1:nc]
    colnames(trim_A.data) = trimmed_target
    
    comps <- read.table(comp_file, sep = "\t", as.is = T, header = F, 
        strip.white = T)
         for (i in 1:length(comps$V1)) {
          print(comps$V1[i])
        fg <- strsplit(comps[,1],"=")
        ci <- fg[[i]][1]
####need to do this for all comparisons
        #filteredx <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_annotated_filtered_DE.txt",sep=""),row.names = 1, header = T, sep = "\t") 
#filteredx <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_annotated_filtered_DE.txt",sep=""), header = T, sep = "\t",as.is=T) 
x <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_Filtered.txt",sep=""), header = T, sep = "\t",as.is=T) 

    names <- trimmed_target
    
    #x <- read.table(paste(readdir,"/","DifferentialExpression/",de_genes,sep=""), blank.lines.skip = FALSE, sep = "\t", 
        #as.is = TRUE, header = TRUE)
    v = x$probeid
    length_v <- length(v) ####Need to impose a limit to the number of rows
    if (length_v > 25) {
    z <- v[1:20]
     v <- z
    }
    sym <- x$symbol
    length_sym <- length(sym)
    if (length_sym > 25) {
    sy <-sym[1:20]
    sym <- sy
    }
   esetSel = trim_A.data[v, ]  ###problem here extracting the relevant rows
    type_e <- esetSel
    print ("THIS IS THE SIZE OF TYPE_e")
    print (type_e)
    #trim_esetSel <- esetSel[number_rows:nc]
    trimmed_size <- dim(esetSel)
    print ("THIS IS THE SIZE OF trimmed_size")
    print (trimmed_size)
    colnames(esetSel) <- trimmed_target
    #trimmed_names <- colnames(trim_esetSel)
    #symbol <- as.character(unlist(lapply(mget(v, env = get(paste(platform, 
     #   "SYMBOL", sep = ""))), function(symbol) {
      #  return(paste(symbol, collapse = "; "))
    #})))
    es = as.matrix(esetSel)
    print ("THIS IS THE EXPRESSION MATRIX")
    #print (es)
    es_size <- dim(es)
    print ("THIS IS THE SIZE OF es")
    print (es_size)
    #es = as.numeric(esetSel)
    colnames(es) = trimmed_target
    print ("I GOT TO HERE TRIMMED ES")
    png(paste(readdir, "/DifferentialExpression/",ci,".heatmap.png",sep=""))
    res=900
    heatmap.2(es, col = bluered, scale = "row", key = "TRUE", 
        symkey = "FALSE", density.info = "none", trace = "none", 
        cexRow = 0.8, cexCol = 1.0, sepcolor="black",labRow = sym,labCol=trimmed_target)
    dev.off()
    
    pdf(paste(readdir, "/DifferentialExpression/",ci,".heatmap.pdf",sep=""))
    res=900
    heatmap.2(es, col = bluered(75), scale = "row", key = "TRUE", 
        symkey = "FALSE", density.info = "none", trace = "none", 
        cexRow = 0.8, cexCol = 1.0, labRow = sym,labCol=trimmed_target)
    dev.off()
    
}
}
