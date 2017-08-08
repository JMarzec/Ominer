############################################################################################################
#
#generate_CNA_cluster.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to generate clusre plot (unsupervised) from filtered .txt file of gains and losses 
###relies on code.R 
getCluster <- function(Cancer,target,results) {
	#setwd(paste(resultdir,Cancer,"/Cluster/",sep=""))
	source("code.R")
	
	 targs <- read.table(target, sep = "\t", as.is = TRUE)
    targets = targs[1, ]
    pd <- read.table(targets, sep = "\t", as.is = TRUE, strip.white = TRUE, 
        header = T)
results.filtered <- read.table(results, sep = "\t", as.is = TRUE)

path_to_change <- paste("ominer_results/",Cancer,"/Cluster/",sep="")
	setwd(path_to_change)
	
	regions = results.filtered;
	cluster=regions[,4:ncol(regions)]
	#HCL=pd$Target
		fact = pd
		#rownames(fact)=pd[,3]
	#facteurs = fact[colnames(cluster),"Target"]
	facteurs <- fact$Group
	colnames(cluster) <- fact$Name
	require(fpc) 
	#require(A2R) 

	#remove the lines where all = zero
	ind = rowSums(cluster == 0) != ncol(cluster)
	cluster_filtered=cluster[ind,]
#cluster_filtered=cluster
	d.usa <- dist(t(cluster_filtered), "euc") 
	#d.usa <- dist(t(cluster), "euc") 
	h.usa <- hclust(d.usa, method="ward") 

	
	set.seed(1)

	hubertgamma <- rep(0,10)

	png("cluster.png")# Set plotting parameters.
	par(mar=rep(2,4)) 
	A2Rplot(h.usa,
			k=2,
			fact.sup= facteurs,
			box=FALSE,
#criteria=hubertgamma,
			show.labels=TRUE,
			col.up = "gray",
			main=" ")
	dev.off()
	pdf("cluster.pdf")
	par(mar=rep(2,4)) 
	A2Rplot(h.usa,
			k=2,
			fact.sup= facteurs,
			box=FALSE,
			#criteria=hubertgamma,
			show.labels=TRUE,
			col.up = "gray",
			main=" ")
	dev.off()
	setwd("../../../")
}
