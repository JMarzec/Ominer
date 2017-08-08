##########################################################################################
#
#File name:cluster_code_transcript.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#
#
#
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
##########################################################################################
##########################################################################################
#Description: function to generate cluster diagrams from the normalised matrix for transcripts
#targets - full path to target.txt
#n = number of clusters to include in the cluster plot
#project - name given to the analysis
################################################################################################
######

cluster_all_tr<-function(targets,n,project) {
	readdir <- paste("ominer_results",project,sep="/")
	A.data<-read.table(paste(readdir,"/transcript/norm/normalised.txt",sep=""),sep="\t",as.is=T,header=T,strip.white=TRUE)
    #A.data<- new("ExpressionSet",exprs = data.matrix(Normdata))
ncolumns <- ncol(A.data)
newdata <- A.data[6:ncolumns]

pd <- read.table(targets, sep="\t", as.is=TRUE,strip.white=TRUE,header=T);
source("code.R")
#setwd(paste(outdir,"/Cluster",sep=""))
print("cluster")
HCL=pd$Target
d.usa <- dist(t(newdata), "euc") 
h.usa <- hclust(d.usa, method="ward") 
png(paste("ominer_results/",project,"/","transcript/","cluster/","transcript_cluster.png",sep=""))
h.usa$labels=HCL
A2Rplot(h.usa,
k=2,
fact.sup=HCL, 
box=FALSE,
#criteria=hubertgamma,
show.labels=TRUE,
col.up = "gray",
main=" ")
dev.off()
}
