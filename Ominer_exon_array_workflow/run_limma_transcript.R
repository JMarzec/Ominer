##########################################################################################
#
#File name: run_limma_transcript.R
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
#Description: Script to identify differentially expressed transcripts
#Arguments are:
#target = path to the target file 
#analysis = paired/unpaired for a (1) paired or an (2) unpaired analysis
#data = "filtered_data.txt" this is the output of the function filtering
#comp = full path to the comaprisons file
#replicates = yes/no yes if the dataset contains repliactes and no if there are no replicates
#limmamethod = method used within limma to compare between different groups e.g. separate
#adjust = method to adjust the false discovery rate (FDR) e.g BH, BY etc.
#pvalue = pvalue threshold used to report genes that are differentially expressed and meet the cutoff threshold(s)
#foldchange = log2 fold change value used to report genes that meet the cutoff threshold(s)

##########################################################################################


#####target = targets file - txt file pd
#####comp = text file containing comparisons that need to be made
#### adjust = method of adjusting FDR
####output = name of output file
####data text file containing normalised data
####Also use this to analyse FIRMA scores from aroma ---load the matrix of FIRMA scores as input to the function instead
run_limma_transcript <- function(target,analysis,data,comp,replicates,limmamethod,adjust,pvalue,foldchange) {
	library("limma")
	library("affy")
   	pd<-read.AnnotatedDataFrame(target,header=T, row.name="Name",sep="\t") # read in targets file
	   readdir <- paste("ominer_results",project,sep="/")
	   Normdata <- read.table(paste(readdir,"/","transcript/norm/",data,sep=""), sep = "\t", as.is = T, header = T, strip.white = TRUE)
ncolumns <- ncol(Normdata)
newdata <- Normdata[6:ncolumns]
allnewdata <- log2(newdata)
#new_rowname <- paste(Normdata$groupName,"-",Normdata$unitName,sep="")
#new_rowname <- paste(Normdata$unitName,"_at",sep="")
#rownames(allnewdata) = new_rowname
to_keep <- cbind(Normdata$unitName,Normdata$groupName)
colnames(to_keep) <- c("probe","transcript")
rownames(allnewdata) = Normdata$groupName
 #write.table(newdata, paste(readdir,"/norm","/", data, sep=""),
        #sep = "\t", row.names = TRUE, quote = FALSE)
A.data<- new("ExpressionSet",exprs = as.matrix(allnewdata))
    sampleNames(A.data)
    names <- featureNames(A.data)
comparisons<-read.table(comp,sep="\t",as.is=T,header=F,strip.white=T)
	lev <- unique(pd$Target)
    lev
    f <- factor(pd$Target, levels = lev)
    # Allow for pairs and replicates
    if(analysis == "unpaired") {
        design <- model.matrix(~0 + f)
        colnames(design) <-  make.names(lev)
    }else {
        
        #if you have a paired analysis you need to add one parameter for the pairs grouping
        levp=unique(pd$Pairs)
        p=factor(pd$Pairs, levels = levp)
        design <- model.matrix(~ -1+f+p)
        colnames(design)<-sub("f","",colnames(design));
    }
    
    #       fit <- lmFit(A.data, design)
    
    if(replicates=="yes"){
        block <- pd$Replicates
        dupcor <- duplicateCorrelation(A.data,design=design,block=block)
        tryCatch(fit <- lmFit(A.data, design,block=block,correlation=dupcor$consensus)
        ,error=function(err) {writeLines(err$message,fileErr);
            cat("Error with lmfit of replicates",file=errorfile,sep="\n",append=TRUE);
            
        })
        fit <- lmFit(A.data, design,block=block,correlation=dupcor$consensus)
    } else {
        fit <- lmFit(A.data, design)
    }
    
    
    myContrasts<-paste(c(comparisons$V1),collapse=",")
    prestr="makeContrasts("
    poststr=",levels=design)"
    commandstr=paste(prestr,myContrasts,poststr,sep="")
    
    #cont.dif<-makeContrasts(contrasts=myContrasts,levels=design)
    cont.dif<-eval(parse(text=commandstr))
    
    #cont.dif <- makeContrasts(
    #comp[1,],
    #levels = design)
    fit2 <- contrasts.fit(fit , cont.dif)
    ebA <- eBayes(fit2)
    
    ebA.decideTests <- decideTests(ebA,method=limmamethod,adjust.method=adjust,p.value=pvalue,lfc=foldchange)
    write.table(ebA.decideTests, paste(readdir,"/transcript/DifferentialExpression","/ebA_allresults.txt", sep=""),
        sep = "\t", row.names = TRUE, quote = FALSE)

    dim(ebA.decideTests@.Data)
    summary(ebA.decideTests)
    sum<-summary(ebA.decideTests)
     write.table(sum, paste(readdir, "/transcript/DifferentialExpression","/decideTestsSummary.txt",sep=""), sep = "\t", 
        row.names = TRUE, quote = FALSE)
    write.table(ebA, paste(readdir, "/transcript/DifferentialExpression","/ebA.txt", sep=""),sep = "\t", row.names = TRUE, 
        quote = FALSE)

      
    comps=colnames(ebA.decideTests)
    # Draw histogram to decide if BH adjustment is valid   - do this per contrast if there is >1
    ##pdf("hist.pdf")
    
    for (i in 1:length(comps)) {
    	pdf(paste(readdir,"/transcript/DifferentialExpression","/hist", comps[i], ".pdf", sep = ""))
               hist(ebA$p.value[,i],breaks=100,col="orange",main=paste("Histogram for ",comps[i],sep=""),xlab='pvalue')
        dev.off()
    }
    
#annot <- read.table("Affy_HuEx1ST_annot_collapsed.txt",header=T,sep="\t")

    comps=colnames(ebA.decideTests)
    for (i in 1:length(comps)) {
        print(comps[i])
        x <- topTable(ebA, coef=comps[i], adjust.method=adjust, sort.by="logFC", number=nrow(A.data))
        probeid = unique(x$ID)
         filtered_x <- topTable(ebA, coef = comps[i], n = 1000, 
            genelist = fit$genes[, 1], adjust.method = adjust, 
            p.value = pvalue, lfc = foldchange, sort.by = "P")
            
            
           # collect_all <- NULL
            #tset_all <- NULL
   #test_id <- NULL
   #nr <- nrow(filtered_x)
   #for (i in 1:nr)  {
   	#test_id <- as.numeric(strsplit(as.character(filtered_x$ID[i]),"_at")[[1]])
   	#test_id <- as.numeric(strsplit(as.character(filtered_x$ID[i]),"-")[[1]])

   	#pset_id <- test_id[[1]][1]
   	#tcluster_id <- test_id[[2]][1] 
   	#collect_all <- c(tcluster_id,collect_all)
   	#collect_all <- c(tcluster_id,collect_all)   ######this was tcluster_id try changing to pset_id
   	#exon_id_annot <- annot[annot$probeset_id %in% pset_id, ]
   	#all_line <- cbind(all_filtered[i],exon_id_annot)
   	#print (all_line)
#collect_all <- c(pset_id,collect_all)
#tset_all <- c(tcluster_id,tset_all)
   #}
            	
           if (length(filtered_x > 0)) {

            
            library("biomaRt")
            up <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  theFilters = c("affy_huex_1_0_st_v2")
  #####need to add ucsc and ottg (vega annotation) to "theAttributes" if uscs is used the same region can be output many times if it has different UCSC annotation 
    #theAttributes =c("affy_huex_1_0_st_v2","chromosome_name","band","hgnc_symbol")
 theAttributes =c("affy_huex_1_0_st_v2","chromosome_name","start_position","end_position","ensembl_gene_id","hgnc_symbol")

ampGenes <- getBM(attributes = theAttributes, filters = theFilters, values = filtered_x$ID, mart=up)
#gene_all = merge(x,ampGenes,by.x="ID",by.y="affy_huex_1_0_st_v2",all=T)
new <- cbind(ampGenes$chromosome_name,ampGenes$start_position,ampGenes$end_position)
new_c <- paste("chr",ampGenes$chromosome_name,":",ampGenes$start_position,"-",ampGenes$end_position,sep="")
together <- cbind(ampGenes$affy_huex_1_0_st_v2,new_c,ampGenes$ensembl_gene_id,ampGenes$hgnc_symbol)
        #write.table(x, paste(readdir,"/transcript/DifferentialExpression/",comps[i],"TopTable_all.txt", 
           # sep = ""), sep = "\t", col.names=TRUE,row.names = TRUE, quote = FALSE)


	 #transcript_id_annot <- annot[annot$probeset %in% filtered_x$ID, ]
	 #theFilters = c("affy_huex_1_0_st_v2")
  #####need to add ucsc and ottg (vega annotation) to "theAttributes" if uscs is used the same region can be output many times if it has different UCSC annotation 
    #theAttributes =c("affy_huex_1_0_st_v2","ensembl_gene_id","chromosome_name","start_position","end_position","description","hgnc_symbol","band")
#theAttributes =c("affy_huex_1_0_st_v2","band","chromosome_name","hgnc_symbol","ensembl_gene_id")
#ampGenes <- getBM(attributes = theAttributes, filters = theFilters, values = filtered_x$ID, mart=up)
filtered_x_parts <- cbind(filtered_x$ID,filtered_x$logFC,filtered_x$adj.P.Val)
colnames(filtered_x_parts) <- c("ID","logFC","adj.P.Val")
#transcript_id_annot <- annot[annot$probeset %in% filtered_x$ID, ]

#gene_filtered = merge(filtered_x_parts,ampGenes,by.x="ID",by.y="affy_huex_1_0_st_v2",all=T)
colnames(together) <- c("probeset_id","locus","ensembl_gene_id","symbol")
gene_filtered = merge(together,filtered_x_parts,by.x="probeset_id",by.y="ID",all=TRUE)

######Replace probeset _id with transcript_id

#t_ids <- to_keep[to_keep[,1] %in% gene_filtered$probeset_id, ]
annot <- read.table("anno_tr.txt",header=T,sep="\t")
t_ids <- annot[annot$probeset_id %in% gene_filtered$probeset_id, ]
final_annot <- merge(t_ids,gene_filtered,by.x="probeset_id",by.y="probeset_id")

f_annot_table <- final_annot[,c(2,3,4,5,6,7)]
filtered_table <- subset(f_annot_table,!duplicated(transcript_cluster_id))


#all_new_filtered <- gene_filtered[,c(1,2,4,5,6,8,9)]
#colnames(all_new_filtered) <- c("probeset","locus","ensembl_transcript_id","ensembl_gene_id","symbol","logFC","P.Value")
	        #all_new_filtered <- cbind(gene_filtered$probeset,gene_filtered$locus,gene_filtered$transcript_id,gene_filtered$gene_id,gene_filtered$gene_short_name,gene_filtered$logFC,gene_filtered$adj.P.Val)
	        c#olnames(all_new_filtered) <- c("probeset","locus","transcript_id","gene_id","gene_name","logFC","adj.P.Val")
	        #all_new <- cbind(gene_filtered$chromosome_no,gene_filtered$band)
	        
	        #cols <- c('band','chromosome_no')
	        #df <- apply([ ,cols] ,1,paste,collapse =":")
	        
	        
	        
	        #new_format_table <- cbind(gene_filtered$probeset,df,gene_filtered$description,gene_filtered$name,gene_filtered$logFC,gene_filtered$adjP.val)
	        #colnames(new_format_table) <- c("probeset","locus","description","name","logFC","adjP.val")
             write.table(filtered_x, paste(readdir,"/transcript/DifferentialExpression/", comps[i],"_Filtered_unannotated.txt", 
            sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
            write.table(filtered_table, paste(readdir,"/transcript/DifferentialExpression/", comps[i],"_Filtered.txt", 
            sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)

           }
}


	
	
}
