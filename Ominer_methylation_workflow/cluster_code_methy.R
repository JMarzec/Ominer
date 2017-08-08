############################################################################################################
#
#File name cluster_code_mrthy.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: script to cluster data from normalised methylation matrix
###inputs are - (1) limma target file - this is the name of teh target file
###k - this is the number of clusters
###project - name given to project


cluster_all <-function(limma_target_file,k,project) {
	readdir <- paste("ominer_results",project,"norm",sep="/")
	A.data<-read.table(paste(readdir,"/normalised.txt",sep=""),sep="\t",as.is=T,header=T,strip.white=TRUE,row.names=1)
    #A.data<- new("ExpressionSet",exprs = data.matrix(Normdata))

pd <- read.table(limma_target_file, sep="\t", as.is=TRUE,strip.white=TRUE,header=T);
source("code.R")
#setwd(paste(outdir,"/Cluster",sep=""))
print("cluster")
HCL=pd$Name
d.usa <- dist(t(A.data), "euc") 
h.usa <- hclust(d.usa, method="ward") 
png(paste("ominer_results/",project,"/","cluster/","cluster.png",sep=""))
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
