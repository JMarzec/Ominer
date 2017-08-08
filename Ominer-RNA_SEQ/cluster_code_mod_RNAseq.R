cluster_all <-function(target_file,k,project) {
	readdir <- paste("ominer_results",project,"norm",sep="/")
	A.data<-read.table(paste(readdir,"/filtered_data.txt",sep=""),sep="\t",as.is=T,header=T,strip.white=TRUE,row.names=1)
    #A.data<- new("ExpressionSet",exprs = data.matrix(Normdata))
require(fpc)
pd <- read.table(target_file, sep="\t", as.is=TRUE,strip.white=TRUE,header=T);
source("code.R")
#setwd(paste(outdir,"/Cluster",sep=""))
print("cluster")
HCL=pd$Target
d.usa <- dist(t(A.data), "euc") 
h.usa <- hclust(d.usa, method="ward") 
pdf(paste("ominer_results/",project,"/","cluster/","cluster.pdf",sep=""))
res = 900
h.usa$labels=HCL
A2Rplot(h.usa,
k,
fact.sup=HCL, 
box=FALSE,
#criteria=hubertgamma,
show.labels=TRUE,
col.up = "gray",
main=" ")
dev.off()

png(paste("ominer_results/",project,"/","cluster/","cluster.png",sep=""),width = 600, height = 600)
res = 900
h.usa$labels=HCL
A2Rplot(h.usa,
k,
fact.sup=HCL, 
box=FALSE,
#criteria=hubertgamma,
show.labels=TRUE,
col.up = "gray",
main=" ")
dev.off()


}
