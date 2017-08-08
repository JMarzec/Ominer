##########################################################################################
#
#File name: run_limma_exon.R
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
#Description: Script to identify differentially expressed exons
#Arguments are:
#targets = path to the target file 
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
run_limma_exon <- function(target,analysis,data,comp,replicates,limmamethod,adjust,pvalue,foldchange) {
	library("limma")
	library("affy")
   	pd<-read.AnnotatedDataFrame(target,header=T, row.name="Name",sep="\t") # read in targets file
	   readdir <- paste("ominer_results",project,sep="/")
	   Normdata <- read.table(paste(readdir,"/","exon/norm/",data,sep=""), sep = "\t", as.is = T, header = T, strip.white = TRUE)
ncolumns <- ncol(Normdata)
newdata <- Normdata[6:ncolumns]
allnewdata <- log2(newdata)
new_rowname <- paste(Normdata$groupName,"-",Normdata$unitName,sep="")
#new_rowname <- paste(Normdata$unitName,"_at",sep="")
rownames(allnewdata) = new_rowname
 #write.table(allnewdata, paste(readdir,"/norm","/newexFit.txt", sep=""),
        #sep = "\t", row.names = TRUE, quote = FALSE)
A.data<- new("ExpressionSet",exprs = data.matrix(allnewdata))
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
    write.table(ebA.decideTests, paste(readdir,"/exon/DifferentialExpression","/ebA_allresults.txt", sep=""),
        sep = "\t", row.names = TRUE, quote = FALSE)

    dim(ebA.decideTests@.Data)
    summary(ebA.decideTests)
    sum<-summary(ebA.decideTests)
     write.table(sum, paste(readdir, "/exon/DifferentialExpression","/decideTestsSummary.txt",sep=""), sep = "\t", 
        row.names = TRUE, quote = FALSE)
    write.table(ebA, paste(readdir, "/exon/DifferentialExpression","/ebA.txt", sep=""),sep = "\t", row.names = TRUE, 
        quote = FALSE)

      
    comps=colnames(ebA.decideTests)
    # Draw histogram to decide if BH adjustment is valid   - do this per contrast if there is >1
    ##pdf("hist.pdf")
    
    for (i in 1:length(comps)) {
    	pdf(paste(readdir,"/exon/DifferentialExpression","/hist", comps[i], ".pdf", sep = ""))
               hist(ebA$p.value[,i],breaks=100,col="orange",main=paste("Histogram for ",comps[i],sep=""),xlab='pvalue')
        dev.off()
    }
    
#annot <- read.table("HuEx-1_0-st-v2.na32.hg19.probeset.csv",header=T,sep=",")

#######Annotation - current annotation problem is that the file being read in for annot has  foreach probeset ID four different exons but only the same exon ids listed in this table....
#######positions these are the positions of the probeset(s) not the gene/transcript or exon - this table in contrast to the table produced by Alex regarding positions list the start
#######location as the start co-ordinates of the first probeset with the stop position listed as the last position of the last probeset.
########Use this table for now but if it is used need to reformat the output table wanted to remove the locations that are listed as these are not for each individual exon and make sure that each exon ID is
########listed only once..................




annot <- read.table("Affy_HuEx1ST_annot_collapsed.txt",header=T,sep="\t")
new_annot <- annot[1:6]
    comps=colnames(ebA.decideTests)
    for (i in 1:length(comps)) {
        print(comps[i])
        x <- topTable(ebA, coef=comps[i], adjust.method=adjust, sort.by="logFC", number=nrow(A.data))
        probeid = unique(x$ID)
        #colnames(x) <- c("ID","logFC","AveExpr","t","P.Value","adj.P.Val","B")
        
        
         filtered_x <- topTable(ebA, coef = comps[i], n = 1000, 
            genelist = fit$genes[, 1], adjust.method = adjust, 
            p.value = pvalue, lfc = foldchange, sort.by = "P")
            
       #exon_annotation <-  read.table("hg19_mod.txt",header = T,sep="\t")        

   #exon_query <- filtered_x$ID
   #exon_id_annot <- exon_annotation[exon_annotation$probe %in% exon_query,]
  
  collect_all <- NULL
   #test_id <- NULL
   nr <- nrow(filtered_x)
   for (i in 1:nr)  {
   	test_id <- as.numeric(strsplit(as.character(filtered_x$ID[i]),"-")[[1]])
   	tcluster_id <- test_id[[1]][1]
   	pset_id <- test_id[[2]][1] 
   	#collect_all <- c(tcluster_id,collect_all)
   	collect_all <- c(tcluster_id,collect_all)   ######this was tcluster_id try changing to pset_id
   	#exon_id_annot <- annot[annot$probeset_id %in% pset_id, ]
   	#all_line <- cbind(all_filtered[i],exon_id_annot)
   	#print (all_line)
#collect_all <- c(pset_id,collect_all)
   }
            	#


#exon_id_annot <- annot[annot$probeset_id %in% collect_all, ]

exon_id_annot <- new_annot[new_annot$probeset %in% collect_all, ]

#allg <- gsub("//",",",exon_id_annot$gene_assignment)
#new_g <- gsub("/","",allg)

#all_exon_annotated <- cbind(exon_id_annot$transcript_cluster_id,exon_id_annot$seqname,exon_id$start,exon_id$stop,exon_id_annot$strand,new_g)
#all_exon_annotated <-cbind(exon_id_annot$cell,exon_id_annot$locus,exon_id_annot$gene_id,exon_id_annot$gene_short_name)
 #new_merged <- merge(exon_id_annot,filtered_x,by.x="cell",by.y="ID",all=TRUE)

#






new_matrix <- cbind(collect_all,filtered_x$logFC,filtered_x$adj.P.Val)
nm <- as.data.frame(new_matrix)

colnames(nm) <- c("ID","logFC","adj.P.Val")

#new_merged <- merge(exon_id_annot,nm,by.x="cell",by.y="ID",all=TRUE)
new_merged <- merge(exon_id_annot,nm,by.x="probeset",by.y="ID",all=TRUE)
colnames(new_merged) <- c("probeset","locus","ensembl_exon_id","ensembl_transcript_id","ensembl_gene_id","symbol","logFC","P.Value")
annot <- read.table("anno_tr.txt",header=T,sep="\t")
t_ids <- annot[annot$probeset_id %in% new_merged$probeset, ]
final_annot <- merge(t_ids,new_merged,by.x="probeset_id",by.y="probeset")
output_table <- final_annot[,c(1,2,3,4,6,7,8,9)]


#wanted <- cbind(new_merged[1],new_merged[3],new_merged[4],new_merged[5],new_merged[8],new_merged[9])
#colnames(wanted) <- c("ID","locus","ensembl_id","HGNC_symbol","logFC","adj.P.Val")

#colnames(all_exon_annotated) <- c("ID","chromosome_no","start","stop","strand","gene_information")
#colnames(all_exon_annotated) <- c("ID","locus","ENSEMBL_ID","gene_symbol")

#sorted_new_matrix <- new_matrix[order(new_matrix[,1]),]

#all_filtered_annot <- cbind(all_exon_annotated,sorted_new_matrix)
 for (i in 1:length(comps)) {
        print(comps[i])

        #colnames(filtered_x) <- c("ID","logFC","AveExpr","t","P.Value","adj.P.Val","B")   
            #modify_rownames_x <- rownames(x)
           # split_rownames_x <- strsplit(modify_rownames_x,"-")
            #nr <- nrow(A.data)
            #pkeep <- NULL
            #for (i in length(nr)) {
            #	pkeep <- c(pkeep,split_rownames_x[[i]][1])
            #}
            #print (comps[i])
            
            #names_keep <- split_rownames[1]
             write.table(output_table, paste(readdir,"/exon/DifferentialExpression/", comps[i],"_Filtered.txt", 
            sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)

                  # write.table(x, paste(readdir,"/DifferentialExpression/","exon/", comps[i],"TopTable_all.txt", 
            #sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
        
             write.table(filtered_x, paste(readdir,"/exon/DifferentialExpression/",comps[i],"_Filtered_unannotated.txt", 
            sep = ""), sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
            }
            }



	
	
}
