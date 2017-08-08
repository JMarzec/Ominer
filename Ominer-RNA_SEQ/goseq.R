##########################################################################################
#
#File name:goseq.R
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
#####project - this is the name given to the analysis
#####comp - this is the .txt file containing information regarding the different comparisons to be made: sample comp_file is comp.txt
#####diffmethod - the differential expression method used in the analysis 'limma' if LIMMA was used and 'edge' if edgeR was used
#####data - filtered_data.txt - this is the filtered normalised expression matrix
#####This code takes in the ebA object generated from limma/edgeR analysis 
#####Output generated are txt/html files of statistically significant Gene Ontology terms that are both over and under represented in EACH of the three ontologies these are:Cellular Component (CC), Biological Process (BP) and Molecular Function (MF)
run_rna_seq_GO <-function(project,comp,diffmethod,data) {
readdir <- paste("ominer_results",project,sep="/")
Normdata <- read.table(paste(readdir,"/norm/",data,sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1)
comps <- read.table(comp, sep = "\t", as.is = T, header = F, 
        strip.white = T)
 for (i in 1:length(comps$V1)) {
                print(comps$V1[i])
        fg <- strsplit(comps[,1],"=")
        ci <- fg[[i]][1]
	print(ci)
####need to do this for all comparisons
tester <- grepl("^ENSG",rownames(Normdata))
     if (tester[1] == "TRUE") {
        tested <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_unannotated.txt",sep=""),row.names = 1, header = T, sep = "\t") 
}
else {
 tested <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"_all_GO.txt",sep=""), header = T, sep = "\t")
 bn <- tested[!duplicated(tested$probe_id),]
tested <- cbind(bn$locus,bn$symbol,bn$logFC,bn$p_value)
rownames(tested) <- bn$probe_id
ty <- as.data.frame(tested)
 genes=as.integer(p.adjust(ty$V4[ty$V3!=0],
 method="BH")<.05)
 }
library("goseq")
library("GO.db") 
library("xtable")
#tested <- read.table("TumourvsNormal_Filtered.txt",sep = "\t", as.is = T, header = T, 
        #strip.white = TRUE, row.names = 1) ####limma output
      
      if (tester[1] == "TRUE") {
 if (diffmethod == "limma") {
  genes=as.integer(p.adjust(tested$P.Value[tested$logFC!=0],
 method="BH")<.05)
 }else {
 	  genes=as.integer(p.adjust(tested$PValue[tested$logFC!=0],
 method="BH")<.05)
 	
 }
 }
 tester <- grepl("^ENSG",rownames(Normdata))
     if (tester[1] == "TRUE") {
    names(genes)=row.names(tested[tested$logFC!=0,])
    }
    else {
      names(genes)=row.names(ty[ty$V3!=0,])
    }
 #  table(genes)      
pwf=nullp(genes,"hg18","ensGene",plot.fit=FALSE)
#####calculate the over and under expreed GO categories among DE genes
#G#O.wall=goseq(pwf,"hg18","ensGene")
#sigCats <-GO.wall[which(GO.wall[,2] < 1.1),]
#cats <- sigCats$category
#terms <- stack(lapply(mget(cats, GOTERM, ifnotfound=NA), Term))
#sigCats$Term <- with(sigCats, terms$values[match(terms$ind, sigCats$category)] )
#
#allGOS <-stack(getgo(rownames(tested),'hg18',"ensGene"))
#onlySigCats <- allGOS[allGOS$values %in% sigCats$category,]
#onlySigCats$Term <- with( onlySigCats, terms$value[match(onlySigCats$values, terms$ind)] )
#onlySigCats$Symbol <- with( onlySigCats, tested[,7][match(onlySigCats$ind, rownames(tested) )] )
#onlySigCats$logFC <- with( onlySigCats, tested$logFC[match(onlySigCats$ind, rownames(tested) )] )

####Now can write out the terms with all information similar to current GO output with affy pipeline for ominer


# write.table(sigCats,paste("ominer_results/",project,"/GO/",ci,"_significant_GO_categories.txt",sep=""), sep = "\t", row.names = FALSE, #quote = FALSE)
 ######need to find a wayto create a html report from this 
#write.table(onlySigCats,paste("ominer_results/",project,"/GO/",ci,"_only_GO_categories.txt", sep=""),sep = "\t", row.names = FALSE, quote = #FALSE)


######Test according to category AND filter the results by over AND under represneted terms in EACH of the three categories
#######Run interpretation of results for molecular function analysis
GO.MF= goseq(pwf,"hg19","ensGene",test.cats=c("GO:MF"))
sigCats.MF <-GO.MF[which(GO.MF[,2] < 1.1),]
cats <- sigCats.MF$category
terms <- stack(lapply(mget(cats, GOTERM, ifnotfound=NA), Term))
sigCats.MF$Term <- with(sigCats.MF, terms$values[match(terms$ind, sigCats.MF$category)] )
over.sig.MF <- subset(sigCats.MF, over_represented_pvalue < 0.05,select = c(category,over_represented_pvalue,Term))
reformat_oversig_MF <- cbind(over.sig.MF$category,over.sig.MF$Term,over.sig.MF$over_represented_pvalue)
colnames(reformat_oversig_MF) <- c("GO ID","Term","p-value")


write.table(reformat_oversig_MF,paste("ominer_results/",project,"/GO/",ci,"MFover_results.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
#print.xtable(reformat_oversig_MF, caption="Gene to GO MF  test for under-representation",file=paste("ominer_results/",project,"/GO/",ci,"MF_over.html",sep=""),caption.placement="top")
#print.xtable(reformat_oversig_MF, caption="Gene to GO MF  test for under-representation",file="MF_under.html",caption.placement="top")
#file.create(paste("ominer_results/",project,"/GO/",ci,"MFover.html",sep=""))
                                #htmlReport(reformat_oversig_MF,file=paste("ominer_results/",project,"/GO/",ci,"MFover.html",sep=""))
                                
#total_rows <- nrow(reformat_oversig_MF)  
if (length(reformat_oversig_MF[,2]) > 2) {                              
  for (i in 1:length(reformat_oversig_MF[,2])){
  	print(i)
reformat_oversig_MF[i,2] <- paste('<a href="http://amigo.geneontology.org/amigo/term/',reformat_oversig_MF[i,1],'">',reformat_oversig_MF[i,2],'</a>',sep='')
#http://amigo.geneontology.org/amigo/term/GO:0004869
}
print(xtable(reformat_oversig_MF, caption="Gene to GO MF  test for over-representation",align="llrr",display=c("s","s","e","e")),type="html",file=paste("ominer_results/",project,"/GO/",ci,"MF_over.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)   ####this one

}
 #if (length(reformat_oversig_MF[,2]) < 1) {   
		#print(xtable(reformat_oversig_MF, caption="Gene to GO MF  test for over-representation"),type="html",file=paste("ominer_results/",project,"/GO/",ci,"MF_over.html",sep=""),caption.placement="top")

	
#}   
   
   #print(xtable(reformat_oversig_MF, caption="Gene to GO MF  test for over-representation",type="html",file=paste("ominer_results/",project,"/GO/",ci,"MF_over.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)
                           
#print(xtable(reformat_oversig_MF, caption="Gene to GO MF  test for over-representation"),type="html",file=paste("ominer_results/",project,"/GO/",ci,"MF_over.html",sep=""),caption.placement="top")
under.sig.MF <- subset(sigCats.MF, under_represented_pvalue < 0.05,select = c(category,under_represented_pvalue,Term))
reformat_undersig_MF <- cbind(under.sig.MF$category,under.sig.MF$Term,under.sig.MF$under_represented_pvalue)
colnames(reformat_undersig_MF) <- c("GO ID","Term","p-value")


write.table(reformat_undersig_MF,paste("ominer_results/",project,"/GO/",ci,"MFunder_results.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
#tr <- xtable(reformat_undersig_MF)
#print(xtable(reformat_undersig_MF, caption="Gene to GO MF  test for under-representation"),type="html",file=paste("ominer_results/",project,"/GO/",ci,"MF_under.html",sep=""),caption.placement="top")
#print(tr, caption="Gene to GO MF  test for under-representation",type="html",file=paste("ominer_results/",project,"/GO/",ci,"MF_under.html",sep=""),caption.placement="top")
 if (length(reformat_undersig_MF[,2]) > 2) {
 for (i in 1:length(reformat_undersig_MF[,2])){
  	print(i)
reformat_undersig_MF[i,2] <- paste('<a href="http://amigo.geneontology.org/amigo/term/',reformat_undersig_MF[i,1],'">',reformat_undersig_MF[i,2],'</a>',sep='')
}
print(xtable(reformat_undersig_MF, caption="Gene to GO MF  test for over-representation",align="llrr",display=c("s","s","e","e")),type="html",file=paste("ominer_results/",project,"/GO/",ci,"MF_under.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)   ####this one
}
#else {
	#print(xtable(reformat_undersig_MF, caption="Gene to GO MF  test for under-representation"),type="html",file=paste("ominer_results/",project,"/GO/",ci,"MF_under.html",sep=""),caption.placement="top")

	
#}

   #print(xtable(reformat_oversig_MF, caption="Gene to GO MF  test for over-representation",type="html",file=paste("ominer_results/",project,"/GO/",ci,"MF_over.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)
######Run interpretation for Cellular Component
GO.CC= goseq(pwf,"hg19","ensGene",test.cats=c("GO:CC"))
sigCats.CC <-GO.CC[which(GO.CC[,2] < 1.1),]
cats <- sigCats.CC$category
terms <- stack(lapply(mget(cats, GOTERM, ifnotfound=NA), Term))
sigCats.CC$Term <- with(sigCats.CC, terms$values[match(terms$ind, sigCats.CC$category)] )
over.sig.CC <- subset(sigCats.CC, over_represented_pvalue < 0.05,select = c(category,over_represented_pvalue,Term))
reformat_oversig_CC <- cbind(over.sig.CC$category,over.sig.CC$Term,over.sig.CC$over_represented_pvalue)
colnames(reformat_oversig_CC) <- c("GO ID","Term","p-value")

write.table(reformat_oversig_CC,paste("ominer_results/",project,"/GO/",ci,"CCover_results.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
if (length(reformat_oversig_CC[,2]) > 2) {
 for (i in 1:length(reformat_oversig_CC[,2])){
  	print(i)
reformat_oversig_CC[i,2] <- paste('<a href="http://amigo.geneontology.org/amigo/term/',reformat_oversig_CC[i,1],'">',reformat_oversig_CC[i,2],'</a>',sep='')
}
print(xtable(reformat_oversig_CC, caption="Gene to GO CC  test for over-representation",align="llrr",display=c("s","s","e","e")),type="html",file=paste("ominer_results/",project,"/GO/",ci,"CC_over.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)   ####this one
}
#else {
	#print(xtable(reformat_oversig_CC, caption="Gene to GO CC  test for over-representation"),type="html",file=paste("ominer_results/",project,"/GO/",ci,"CC_over.html",sep=""),caption.placement="top")

	
#}


#print(xtable(reformat_oversig_CC, caption="Gene to GO CC  test for over-representation"),type="html",file=paste("ominer_results/",project,"/GO/",ci,"CC_over.html",sep=""),caption.placement="top")
#print.xtable(reformat_oversig_CC, caption="Gene to GO CC  test for over-representation",file=paste("ominer_results/",project,"/GO/",ci,"CC_over.html",sep=""),caption.placement="top")

under.sig.CC <- subset(sigCats.CC, under_represented_pvalue < 0.05,select = c(category,under_represented_pvalue,Term))
reformat_undersig_CC <- cbind(under.sig.CC$category,under.sig.CC$Term,under.sig.CC$under_represented_pvalue)
colnames(reformat_undersig_CC) <- c("GO ID","Term","p-value")


write.table(reformat_undersig_CC,paste("ominer_results/",project,"/GO/",ci,"CCunder_results.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
#print.xtable(reformat_undersig_CC, caption="Gene to GO CC  test for under-representation",file=paste("ominer_results/",project,"/GO/",ci,"CC_under.html",sep=""),caption.placement="top")
if (length(reformat_undersig_CC[,2]) > 2) {
 for (i in 1:length(reformat_undersig_CC[,2])){
  	print(i)
reformat_undersig_CC[i,2] <- paste('<a href="http://amigo.geneontology.org/amigo/term/',reformat_undersig_CC[i,1],'">',reformat_undersig_CC[i,2],'</a>',sep='')
}
print(xtable(reformat_undersig_CC, caption="Gene to GO CC  test for under-representation",align="llrr",display=c("s","s","e","e")),type="html",file=paste("ominer_results/",project,"/GO/",ci,"CC_under.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)   ####this one
}
#else {
#print(xtable(reformat_undersig_CC, caption="Gene to GO CC  test for under-representation"),type="html",file=paste("ominer_results/",project,"/GO/",ci,"CC_under.html",sep=""),caption.placement="top")
#}





#######Run analysis for Biological Process
GO.BP= goseq(pwf,"hg19","ensGene",test.cats=c("GO:BP"))
sigCats.BP <-GO.BP[which(GO.BP[,2] < 1.1),]
cats <- sigCats.BP$category
terms <- stack(lapply(mget(cats, GOTERM, ifnotfound=NA), Term))
sigCats.BP$Term <- with(sigCats.BP, terms$values[match(terms$ind, sigCats.BP$category)] )
over.sig.BP <- subset(sigCats.BP, over_represented_pvalue < 0.05,select = c(category,over_represented_pvalue,Term))
reformat_oversig_BP <- cbind(over.sig.BP$category,over.sig.BP$Term,over.sig.BP$over_represented_pvalue)
colnames(reformat_oversig_BP) <- c("GO ID","Term","p-value")


write.table(reformat_oversig_BP,paste("ominer_results/",project,"/GO/",ci,"BPover_results.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
#print.xtable(reformat_oversig_BP, caption="Gene to GO BP  test for over-representation",file=paste("ominer_results/",project,"/GO/",ci,"BP_over.html",sep=""),caption.placement="top")

if (length(reformat_oversig_CC[,2]) > 2) {
 for (i in 1:length(reformat_oversig_BP[,2])){
  	print(i)
reformat_oversig_BP[i,2] <- paste('<a href="http://amigo.geneontology.org/amigo/term/',reformat_oversig_BP[i,1],'">',reformat_oversig_BP[i,2],'</a>',sep='')
}
print(xtable(reformat_oversig_BP, caption="Gene to GO BP  test for over-representation",align="llrr",display=c("s","s","e","e")),type="html",file=paste("ominer_results/",project,"/GO/",ci,"BP_over.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)   ####this one
}
#else {
#print(xtable(reformat_oversig_BP, caption="Gene to GO BP  test for over-representation"),type="html",file=paste("ominer_results/",project,"/GO/",ci,"BP_over.html",sep=""),caption.placement="top")
#}
under.sig.BP <- subset(sigCats.BP, under_represented_pvalue < 0.05,select = c(category,under_represented_pvalue,Term))
reformat_undersig_BP <- cbind(under.sig.BP$category,under.sig.BP$Term,under.sig.BP$under_represented_pvalue)
colnames(reformat_undersig_BP) <- c("GO ID","Term","p-value")


write.table(reformat_undersig_BP,paste("ominer_results/",project,"/GO/",ci,"BPunder_results.txt",sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
#print.xtable(reformat_undersig_BP,caption="Gene to GO BP  test for under-representation",file=paste("ominer_results/",project,"/GO/",ci,"BP_under.html",sep=""),caption.placement="top")
if (length(reformat_undersig_BP[,2]) > 2) {
 for (i in 1:length(reformat_undersig_BP[,2])){
  	print(i)
reformat_undersig_BP[i,2] <- paste('<a href="http://amigo.geneontology.org/amigo/term/',reformat_undersig_BP[i,1],'">',reformat_undersig_BP[i,2],'</a>',sep='')
}
print(xtable(reformat_undersig_BP, caption="Gene to GO BP  test for under-representation",align="llrr",display=c("s","s","e","e")),type="html",file=paste("ominer_results/",project,"/GO/",ci,"BP_under.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)   ####this one
}
#else {
#print(xtable(reformat_undersig_BP, caption="Gene to GO BP  test for under-representation"),type="html",file=paste("ominer_results/",project,"/GO/",ci,"BP_under.html",sep=""),caption.placement="top")
#}


		}
		}