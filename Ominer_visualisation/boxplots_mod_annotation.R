############################################################################################################
#
#boxplots_mod_annotation.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to generate boxplots from either a gene symbol or a probeid for each of the expression platforms supported by O-miner
###Arguments - (1) study - this is teh name of the project (2) probeFile - this is the name of the .txt file in which either the gene/symbol or probeid has been written out to from the interface
####(3) filename - this is the name given to the .txt file that is output containing the boxplot (4) type - this is either symbol or probeid, (5) platform  this indicates the platform that is used e.g. 450k or 27k for methylation (5) technology - this is the technoogy that was used i.e. methylation or affy_expr, illumina_expr, affy_mirna and rna_seq, filename this is the name that is to be given to the output file
####Code to generate expression boxplot(s) of the gene of interest from the normalised expression matrix 
####Depends on R libaries - (1) Hmisc and (2) affy (3) biomaRt
####.pdf and .png is generated of the expression boxplots - .png is displayed in the webpage and .pdf is available to download from the webpage - has to be .png that can be viewed from teh webpage 


# Parse the command line options
for (e in commandArgs()) {
	ta = strsplit(e,"=",fixed=TRUE)
	ta
	if(! is.na(ta[[1]][2])) {
		temp = ta[[1]][2]
		assign(ta[[1]][1],temp)
		cat("assigned ",ta[[1]][1]," the value of |",temp,"|\n")
	} 
}

rootdir<-"/var/www/cgi-bin/onlinetool/version_2/ominer_results/"
output_dir <- "/var/www/html/onlinetool/temp/"
library("Hmisc")    ##### downloaded from CRAN
library("affy")

####First thing is to read in the expression matrix for each of the technologies 
if (platform == "450k"||platform == "27k") {
technology <- "methylation"
print ("my technology is methylation")
}



if (technology == "affy_exon") {    #added fo exon_array as A.data to be used is the one from the transcript level
	A.data<-read.table(paste(rootdir,study,"/transcript/norm/normalised.txt",sep=""),sep="\t",as.is=T,header=T)   #####added for exon_array as A.data to be used is from teh transcript level
	rownames(A.data) <- A.data$groupName
e=data.frame(A.data)  #####

	
	
}else {
A.data<-read.table(paste(rootdir,study,"/norm/normalised.txt",sep=""),sep="\t",as.is=T,header=T,row.names=1)
e=data.frame(A.data)
print (" I have just read in the methyaltion matrix")
}
probeList<-read.table(probeFile,as.is=T,header=F)
if (platform == "hugene1.1") {
png(file=paste(output_dir,study,"/",filename,".png",sep=""))
	par(las=2,cex.axis=0.5)
if (type == "symbol") {
probesAnnotFile <- read.table("Affy_HuGene1ST_annot_collapsed.txt",sep="\t",as.is=TRUE,header=TRUE,quote=NULL, comment='')
#converted <- grep(probeList$V1,annot_1.0)
converted <- probesAnnotFile[probesAnnotFile$gene_short_name %in% probeList$V1, ]
print ("THIS IS THE CONVERTED VALUE FOR HUGENE1.1")
print (converted)
rev <- converted$probeset
print ("THIS IS THE AMH PROBESET")
print (rev)
}
}


if (technology == "methylation") {

orig_dir <- paste("/var/www/html/onlinetool/temp/",study,"/","target.txt",sep="")
print (orig_dir)
pd<-read.AnnotatedDataFrame(orig_dir,header=T,row.name="FileName",sep="\t")
print ("I just read this in")
}

if (technology != "affy_exon" && technology != "methylation") {
if(file.exists(paste(rootdir,study,"/QC/target_qc.txt",sep=""))) {
pd<-read.AnnotatedDataFrame(paste(rootdir,study,"/QC/target_qc.txt",sep=""),header=T, row.name="Name",sep="\t")
}
else {
pd<-read.AnnotatedDataFrame(paste(rootdir,study,"/QC/target.txt",sep=""),header=T, row.name="Name",sep="\t")
}
}

if (platform == "hugene1.1") {
pd<-read.AnnotatedDataFrame(paste(rootdir,study,"/QC/target.txt",sep=""),header=T, row.name="FileName",sep="\t")
e=as.data.frame(A.data)
}


if(technology == "affy_exon") {
print ("I am doing this")
pd<-read.AnnotatedDataFrame(paste(rootdir,study,"/transcript/QC/target.txt",sep=""),header=T, row.name="FileName",sep="\t")
}


if (technology != "methylation") {
lev=unique(pd$Target)
}
if (technology == "methylation") {
lev=unique(pd$Name)
}




####If the user has entered a probeid then this can be extracted rom the expesion matrix and an expression boxplot generated 
print ("I have reached here 1")
if (technology == "methylation" && type == "probe") {
pdf(file=paste(output_dir,study,"/",filename,".pdf",sep=""))
	par(las=2,cex.axis=0.5)
pl <- probeList$V1
	print ("this is probe")
	print (pl)
	if(pl %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[pl,sampleNames(pd[pd$Name==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
				print(lev)
					mylist<-c(mylist,list(as.numeric(e[pl,sampleNames(pd[pd$Name==lev[i]])])))
				}
                              boxplot(mylist,cex.axis=1,names=lev,main=pl,col=rainbow(length(lev)))
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}
	dev.off()
	
	
	png(file=paste(output_dir,study,"/",filename,".png",sep=""))
	par(las=2,cex.axis=0.5)
	
	if(pl %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[pl,sampleNames(pd[pd$Name==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
				print(lev)
					mylist<-c(mylist,list(as.numeric(e[pl,sampleNames(pd[pd$Name==lev[i]])])))
				}
                                #boxplot(mylist,xaxt="n",main=pl,names=lev)
boxplot(mylist,cex.axis=1,names=lev,main=pl,col=rainbow(length(lev)))
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}
	
	
	}


####If a symbol has been used and the technology is methylation - then this needs to be looked up - to extract the relevant probeid from teh normalised expression matrix in order to generate the boxplot 
if (technology == "methylation" && type == "symbol") {


pl <- probeList$V1
print(pl)

        	        
                       
       
 
         
methy_lookup <- read.table("methy_gs_lookup.txt",header=TRUE,sep="\t",as.is=T)
rev_list <- methy_lookup[which(methy_lookup$geneSymbol==pl),]
print ("THIS IS MY REVLIST")
print (rev_list)
rev_length <- length(rev_list)
if (rev_length >1) {
print ("LENGTH was > 1")
ty <- typeof(rev_length)
pid <- rev_list$ID[1]
print ("THIS IS MY PROBE ID")
print (pid)
}
else {
pid <- rev_list$ID
}
print("I tried to look up the probe for THS symbol")
print(pid)
#pid <- probe[1][[1]]
print("THIS is the PROBE I AM looking for")
print (pid)
png(file=paste(output_dir,study,"/",filename,".png",sep=""))
	par(las=2,cex.axis=0.5)
	
	if(pid %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[pid,sampleNames(pd[pd$Name==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
				print(lev)
					mylist<-c(mylist,list(as.numeric(e[pid,sampleNames(pd[pd$Name==lev[i]])])))
				}
                               #boxplot(mylist,xaxt="n",main=pid,names=lev)
boxplot(mylist,cex.axis=1,main=pid,names=lev,col=rainbow(length(lev)))
print("THESE ARE MY NAMES")
print (lev)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}
	
	
dev.off()	
pdf(file=paste(output_dir,study,"/",filename,".pdf",sep=""))
        par(las=2,cex.axis=0.5)

	
	
	if(pid %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[pid,sampleNames(pd[pd$Name==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
				print(lev)
					mylist<-c(mylist,list(as.numeric(e[pid,sampleNames(pd[pd$Name==lev[i]])])))
				}
                                #boxplot(mylist,xaxt="n",main=pid,names=lev)
boxplot(mylist,cex.axis=1,main=pid,names=lev,col=rainbow(length(lev)))
print ("THESE ARE MY NAMES PDF")
print (lev)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}

	
	
	}
	





####Repeat if symbol is used for affy_expression technology(s)

if(type == "symbol") {
	if (technology == "affy_expr") {
	if (platform != "hugene1.1") {
	if (platform != "mouse4302") {
library(annotate)
####load the annotation libray that is sepcific to the platform being used 
library(paste(platform,".db",sep=""),character.only=TRUE)
rev=revmap(getAnnMap("SYMBOL",paste(platform,".db",sep="")))




pdf(file=paste(output_dir,study,"/",filename,".pdf",sep=""))
par(las=2,cex.axis=0.5)
for (probe in probeList$V1) {
for(p in mget(probeList$V1,rev,ifnotfound = NA) [[probe]]) {


count = 0
if(p %in% rownames(A.data)) {
mylist <- list(as.numeric(e[p,sampleNames(pd[pd$Target==lev[1]])]))


for(i in 2:length(lev)){
mylist<- c(mylist,list(as.numeric(e[p,sampleNames(pd[pd$Target==lev[i]])])))
}
boxplot(mylist,xaxt="n",main=p,names=lev)
text(1:length(lev),par("usr")[3]-0.1,srt=45,adj=1,labels=lev,xpd=TRUE)
count=1
}

}

}

}
}
}

}

####Repeat for mouse_affymetrix platform and if symbol has been selected 
	if (platform == "mouse4302") {
if (type == "symbol") {
		pl <- probeList$V1
				



m_lookup<-read.table("/var/www/cgi-bin/onlinetool/version_2/mouse_4302.txt",sep="\t",as.is=T,header=T)
rev_list <- m_lookup[which(m_lookup$symbol==pl),]
print ("I HAVE DONE THIS")
print ("I AM GENERATING THE PDF")

pdf(file=paste(output_dir,study,"/",filename,".pdf",sep=""))
par(las=2,cex.axis=0.5)

for(i in 1:length(rev_list$probeid)) {
rev <- rev_list$probeid[i]
		
	
count=0
if(rev %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[i]])])))
				}
                                boxplot(mylist,cex.axis=1,names=lev,main=paste(rev,sep=""),las=2,col = rainbow(length(lev)))
				
				count=1
		}	

	
	
}

}
dev.off()

}

if (type == "symbol") {
	if (platform == "mouse4302") {
		pl <- probeList$V1
		

	

m_lookup<-read.table("/var/www/cgi-bin/onlinetool/version_2/mouse_4302.txt",sep="\t",as.is=T,header=T)
rev_list <- m_lookup[which(m_lookup$symbol==pl),]
print ("I HAVE DONE THIS")
print ("I AM DOING THE PNG")

png(file=paste(output_dir,study,"/",filename,".png",sep=""))
par(las=2,cex.axis=0.5)
i <- 0
for(rev in rev_list$probeid) {
	print ("THIS IS MOUSE REV")
	print (rev)
	
	
count=0
if(rev %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[i]])])))
				}
                                
				boxplot(mylist,cex.axis=1,names=lev,main=paste(rev,sep=""),las=2,col = rainbow(length(lev)))
                                boxplot(mylist,xaxt="n",main=rev,col= rainbow(length(lev)))
				text(1:length(lev), par("usr")[3]-0.02, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}	
	
	
	
	



}

}

dev.list()
graphics.off()
}
if (platform == "mouse4302") {
if (type == "symbol") {
d <- dev.list()
a <- dev.cur()

}
}

if (platform == "mouse4302") {
if (type == "probe") {
pdf(file=paste(output_dir,study,"/",filename,".pdf",sep=""))
par(las=2,cex.axis=0.5)
for (probe in probeList$V1) {
		if(probe %in% rownames(A.data)) {
		ghj <- grep("-",sampleNames(pd))
				k <- length(ghj)
				if (k > 0) {
				colnames(e) <- sampleNames(pd)
				}
				if (platform == "hugene1.1") {
				probe <-as.character(probe)
				}
			mylist<-list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[1]])]))
			pdt <- pd$Target
			print (pdt)
			for(i in 2:length(lev)){
			probe <- as.character(probe)
				mylist<-c(mylist,list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[i]])])))
			}
			boxplot(mylist,xaxt="n",main=probe,col=rainbow(length(lev)))
			text(1:length(lev), par("usr")[3]-0.02, srt = 45, adj = 1,labels = lev, xpd = TRUE)

		
			count=1
		}
	}
dev.off()
}


}
if (platform == "mouse4302") {
if (type == "probe") {
png(file=paste(output_dir,study,"/",filename,".png",sep=""))
par(las=2,cex.axis=0.5)
for (probe in probeList$V1) {
		if(probe %in% rownames(A.data)) {
		ghj <- grep("-",sampleNames(pd))
				k <- length(ghj)
				if (k > 0) {
				colnames(e) <- sampleNames(pd)
				}
				if (platform == "hugene1.1") {
				probe <-as.character(probe)
				}
			mylist<-list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[1]])]))
			pdt <- pd$Target
			print (pdt)
			for(i in 2:length(lev)){
			probe <- as.character(probe)
				mylist<-c(mylist,list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[i]])])))
			}
			boxplot(mylist,xaxt="n",main=probe,col=rainbow(length(lev)))
			text(1:length(lev), par("usr")[3]-0.02, srt = 45, adj = 1,labels = lev, xpd = TRUE)
			


			count=1
		}
	}
dev.off()
dev.list()
}

}


####If the technology was Illumina_expression and symbol was chosen 
if (type == "symbol") {
	if (technology == "illumina_expr") {
		pl <- probeList$V1
type_pl <- typeof(pl)
print ("THIS IS THE TYPE OF PL")
print (type_pl)
length_pl <- length(pl)
if (length_pl == 0 ) {
pl <- probeList$V3
type_pl <- typeof(pl)

}
	####Map the probe_id back to gene_symbol using BioMaRt this is to generate a .txt file with two columns - 1 probeid and the other being gene_symbol 	
		library("biomaRt")
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")

if (platform == "ht12v3") {
	platform_annot = "illumina_humanht_12_v3"
}
if (platform == "ht12v4") {
	platform_annot = "illumina_humanht_12_v4"
}

theFilters = c("hgnc_symbol")   ######or illumina_humanht_v3
theAttributes = c(platform_annot,"hgnc_symbol")
wanted_annot <- getBM(attributes = theAttributes, filters = theFilters, values = probeList$V1, mart=mart)
		print (wanted_annot)
		rev <- wanted_annot[,1][1]
		print (rev)
	
	png(file=paste(output_dir,study,"/",filename,".png",sep=""))
	par(las=2,cex.axis=0.5)
count=0
if(rev %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[i]])])))
				}
                                boxplot(mylist,xaxt="n",main=rev)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}	
	
	
	
	
}
}



if (type == "symbol") {
	if (technology == "illumina_expr") {
		pl <- probeList$V1
		print ("just doing illumina pdf NOW")
		print (pl)
		library("biomaRt")
##wanted<-read.table("illumina_probes_annot.txt",header=T,sep="\t",as.is=T)
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
#ensembl=useMart("ensembl")
#ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
if (platform == "ht12v3") {
	platform_annot = "illumina_humanht_12_v3"
}
if (platform == "ht12v4") {
	platform_annot = "illumina_humanht_12_v4"
}

theFilters = c("hgnc_symbol")   ######or illumina_humanht_v3
theAttributes = c(platform_annot,"hgnc_symbol")
wanted_annot <- getBM(attributes = theAttributes, filters = theFilters, values = probeList$V1, mart=mart)
		print (wanted_annot)
		rev <- wanted_annot[,1][1]
		print (rev)
	
	pdf(file=paste(output_dir,study,"/",filename,".pdf",sep=""))
	#par(las=2,cex.axis=0.5)
count=0
if(rev %in% rownames(A.data)) {
print ("I have opened SYMBOL PDF")
				mylist<-list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[i]])])))
				}
                                boxplot(mylist,xaxt="n",main=rev)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}	
	
	
	
}

}





allprobe<-NULL


####If the technology chosen was RNA-Seq post-processing generate the expression boxplots 

	if (technology == "rna_seq") {
	pdf(file=paste(output_dir,study,"/",filename,".pdf",sep=""))
	par(las=2,cex.axis=0.5)
count=0
	if (type == "symbol") {
			pl <- probeList$V1
print ("here I am")
	tester <- grepl("^ENSG",rownames(A.data))
     if (tester[1] == "TRUE") {
		
print (pl)
		library("biomaRt")
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
theFilters = c("hgnc_symbol")   ######or illumina_humanht_v3
theAttributes = c("ensembl_gene_id")
wanted_annot <- getBM(attributes = theAttributes, filters = theFilters, values = probeList$V1, mart=mart)
				rev <- wanted_annot[1][1]
		
		if(rev %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[1]])]))
				
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[i]])])))
				}
                                boxplot(mylist,xaxt="n",main=rev)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}
	}
	}
	else {
		pl <- probeList$V1
	if(pl %in% rownames(A.data)) {
	
				mylist<-list(as.numeric(e[pl,sampleNames(pd[pd$Target==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[pl,sampleNames(pd[pd$Target==lev[i]])])))
				}
                                boxplot(mylist,xaxt="n",main=pl)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}
	}
	dev.off()
}



allprobe<-NULL

count=0
if (technology == "rna_seq") {
	if (type == "symbol") {
png(file=paste(output_dir,study,"/",filename,".png",sep=""))
par(las=2,cex.axis=0.5)
	tester <- grepl("^ENSG",rownames(A.data))
     if (tester[1] == "TRUE") {
		pl <- probeList$V1
		print ("I am here now 6")
print (pl)
		library("biomaRt")
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org")
theFilters = c("hgnc_symbol")   ######or illumina_humanht_v3
theAttributes = c("ensembl_gene_id")
wanted_annot <- getBM(attributes = theAttributes, filters = theFilters, values = probeList$V1, mart=mart)
		print (wanted_annot)
		rev <- wanted_annot[1][1]
		print (rev)
				if(rev %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[i]])])))
				}
                                boxplot(mylist,xaxt="n",main=rev)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}

	}
	}
	else {
	pl <- probeList$V1
	print ("this is probe")
	print (pl)
	if(pl %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[pl,sampleNames(pd[pd$Target==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[pl,sampleNames(pd[pd$Target==lev[i]])])))
				}
                                boxplot(mylist,xaxt="n",main=pl)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}
	}
dev.off()
}


####If the technology was Illumina_expression and probeid is chsen generate boxplots 

for (probe in probeList$V1) {


if (technology == "illumina_expr") {
png(file=paste(output_dir,study,"/",filename,".png",sep=""))
	
	
	
	
		if (type == "probe") {
                		if(probe %in% rownames(A.data)) {
		print ("yes")
		ghj <- grep("-",sampleNames(pd))
				k <- length(ghj)
				if (k > 0) {
				colnames(e) <- sampleNames(pd)
				}
				mylist<-list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[1]])]))
			
				print (mylist)
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[i]])])))
				}

                                boxplot(mylist,xaxt="n",main=probe)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				
			count=1
		}
	}
        

dev.off()
}

}
####If the technology was Illumina_expression and probeid was chosen 

for (probe in probeList$V1) {

print(probe)
if (technology == "illumina_expr") {
pdf(file=paste(output_dir,study,"/",filename,".pdf",sep=""))
	
	
	
	
		if (type == "probe") {
                		if(probe %in% rownames(A.data)) {
		print ("yes")
		ghj <- grep("-",sampleNames(pd))
				k <- length(ghj)
				if (k > 0) {
				colnames(e) <- sampleNames(pd)
				}
				mylist<-list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[1]])]))
				print (mylist)
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[i]])])))
				}

                                boxplot(mylist,xaxt="n",main=probe)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
							count=1
		}
	}
        
dev.off()

}

}

####If the technology was Illumina_expression and symbol was chosen 

if (type == "symbol") {
	if (technology == "illumina_expr") {
png(file=paste(output_dir,study,"/",filename,".png",sep=""))		
if(rev %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[rev,sampleNames(pd[pd$Target==lev[i]])])))
				}
                                boxplot(mylist,xaxt="n",main=rev)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}
		}

		


	}



####If the technology was affymetrix expression and not mouse4302
for (probe in probeList$V1) {

print(probe)
if (technology == "affy_expr") {
if (platform != "mouse4302") {
pdf(file=paste(output_dir,study,"/",filename,".pdf",sep=""))
	if(type == "symbol") {
		
	
	if (platform != "hugene1.1") {
	if (platform != "mouse4302")
		for(probe in mget(probeList$V1,rev,ifnotfound = NA)[[probe]]) {
		
		if(probe %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[1]])]))
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[i]])])))
				}
                                boxplot(mylist,xaxt="n",main=probe)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
						}
		}
		
	}
	}
	####If the technology is hugene1.1 and type==probe

	if (platform == "hugene1.1") {
	type <- "probe"
	probe <- rev
	}
	}
	
	
		if (type == "probe") {
               
		KJ <- length(rownames(A.data))
		
		if(probe %in% rownames(A.data)) {
		
		ghj <- grep("-",sampleNames(pd))
				k <- length(ghj)
				if (k > 0) {
				colnames(e) <- sampleNames(pd)
				}
				if (platform == "hugene1.1") {
				probe <- as.character(probe)
				}
				mylist<-list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[1]])]))
								cnames <- colnames(e)
				type_e <- typeof(e)
				
				spd <- sampleNames(pd)
				
				pdt <- pd$Target
				
				
				for(i in 2:length(lev)){
				if (platform == "hugene1.1") {
				probe <- as.character(probe)
				}
					mylist<-c(mylist,list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[i]])])))
					
				}

                                boxplot(mylist,xaxt="n",main=probe)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				
			count=1
		}
	}
        

dev.off()
}

}




######If affy_expression is used and symbol - generate expression boxplots 




for (probe in probeList$V1) {
	
		
		if (technology == "affy_expr") {
		if (platform != "mouse4302") {
		png(file=paste(output_dir,study,"/",filename,".png",sep=""))
	par(las=2,cex.axis=0.5)
		if(type == "symbol") {
		for(probe in mget(probeList$V1,rev,ifnotfound = NA)[[probe]]) {
			print(probe)
			if(probe %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[1]])]))
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[i]])])))
				}

                                boxplot(mylist,xaxt="n",main=probe,col= rainbow(length(lev)))
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
				
			}
		}
	
	}
	
	######If affy_expression is used and probeid - generate expression boxplots 

		if (type == "probe") {
		if(probe %in% rownames(A.data)) {
		ghj <- grep("-",sampleNames(pd))
				k <- length(ghj)
				if (k > 0) {
				colnames(e) <- sampleNames(pd)
				}
				if (platform == "hugene1.1") {
				probe <-as.character(probe)
				}
			mylist<-list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[1]])]))
			pdt <- pd$Target
			print (pdt)
			for(i in 2:length(lev)){
			probe <- as.character(probe)
				mylist<-c(mylist,list(as.numeric(e[probe,sampleNames(pd[pd$Target==lev[i]])])))
			}
			boxplot(mylist,xaxt="n",main=probe,col=rainbow(length(lev)))
			text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)



			count=1
		}
	}

dev.off()
}

}

}

######For miRNA and probe & symbol - generate expression boxplots

if (technology == "miRNA") {
if (type == "probe") {
png(file=paste(output_dir,study,"/",filename,".png",sep=""))
par(las=2,cex.axis=0.5)
pl <- probeList$V1
p  <- pl
}
if (type == "symbol"){
pl <- probeList$V1
p <- paste(pl,"_st",sep="")
}
	if(pl %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[pl,sampleNames(pd[pd$Target==lev[1]])]))
				
				for(i in 2:length(lev)){
				
					
					mylist<-c(mylist,list(as.numeric(e[pl,sampleNames(pd[pd$Target==lev[i]])])))
					LEV <- lev[i]
					print (LEV)
					
				}
				boxplot(mylist,cex.axis=1,names=lev,main=paste(p," ",probe,sep=""),las=2)
                                boxplot(mylist,xaxt="n",main=probe)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}

}


	if (technology == "miRNA_old") {
if (type == "symbol") {

png(file=paste(output_dir,study,"/",filename,".png",sep=""))
	par(las=2,cex.axis=0.5)
pl <- probeList$V1
wanted <- paste(pl,"_st",sep="")
if(wanted %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[wanted,sampleNames(pd[pd$Target==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[wanted,sampleNames(pd[pd$Target==lev[i]])])))
				}
                                boxplot(mylist,xaxt="n",main=pl)
				boxplot(mylist,cex.axis=1,names=lev,main=probe)
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}



}
dev.off()
}


	if (technology == "miRNA_old") {
if (type == "symbol") {
pdf(file=paste(output_dir,study,"/",filename,".pdf",sep=""))
	par(las=2,cex.axis=0.5)
pl <- probeList$V1
wanted <- paste(pl,"_st",sep="")
if(wanted %in% rownames(A.data)) {
				mylist<-list(as.numeric(e[wanted,sampleNames(pd[pd$Target==lev[1]])]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
					mylist<-c(mylist,list(as.numeric(e[wanted,sampleNames(pd[pd$Target==lev[i]])])))
				}
                                boxplot(mylist,xaxt="n",main=pl,col=rainbow(length(lev)))
				boxplot(mylist,cex.axis=1,names=lev,main=probe,col=rainbow(length(lev)))
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		}



}
dev.off()
}


if (technology == "affy_exon") {
png(file=paste(output_dir,study,"/",filename,".png",sep=""))
	par(las=2,cex.axis=0.5)
	print ("HERE I AM AFFY EXON!!!")
#file.create(paste(rootdir,study,"/transcript/norm/",filename,".txt",sep=""))
pl <- probeList$V1
colnames(A.data) <- pd$Name
rownames(A.data) <- A.data$groupName
e <- data.frame(A.data)
	if(pl %in% rownames(A.data)) {
	print ("I got to here exon_arrays")
				mylist<-list(as.numeric(e[pl,pd$Name[pd$Target==lev[1]]]))
				print ("this is MY LIST")
				print (mylist)
				for(i in 2:length(lev)){
				mylist<-c(mylist,list(as.numeric(e[pl,pd$Name[pd$Target==lev[i]]])))
					
					print (mylist)
				}
                                boxplot(mylist,xaxt="n",main=pl,col=rainbow(length(lev)))
				text(1:length(lev), par("usr")[3]-0.1, srt = 45, adj = 1,labels = lev, xpd = TRUE)
				count=1
		
		}
dev.off()
}



