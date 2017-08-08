##########################################################################################
#
#File name: bcctb.utils.R
#
#
#
#
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
##########################################################################################
##########################################################################################
#Description: function to match probe id from array platform to an Entrez gene id
#Arguments are:
#exp.data - normalised expression matrix
#platform - is the name of the array playform data is from e.g. hgu133plus2 for hgu133plus2.db
#type - datatype i.e raw CEL files "CEL" and normalised expression matrix - "normalised"
###################################################################################
######

library("genefilter")
library(org.Hs.eg.db)
library(annotate)

getEntrez <- function(exp.data,platform,type) {
        ## Now we need to switch from probe ids to entrez gene ids for each of the probes in the expression matrix
        #### Filter out records wih no entrezIds

        entrezIds <- mget(rownames(exp.data),envir=get(paste(platform,type,sep="")))
        haveEntrezId <- names(entrezIds)[sapply(entrezIds,function(x) !is.na(x))]            
        entrez.data <- exp.data[haveEntrezId,]

        ## map each probe set to exactly one Entrez gene id.  If multiple probes found select probe with largest IQR

        ## Removing rows with NAs here - we lost about 1300 rows
        entrez.data <- entrez.data[complete.cases(entrez.data),]

        data.var <- apply(entrez.data,1,var)
        uniqGenes <-findLargest(rownames(entrez.data),data.var,platform)
        entrez.data <- entrez.data[uniqGenes,]
                                
        entrezIds <- mget(rownames(entrez.data),envir=get(paste(platform,type,sep="")))
        entrezIds.matrix <- data.frame(sapply(entrezIds, "[",  1))
        colnames(entrezIds.matrix) <- c("EntrezId")
        all.data=merge(entrezIds.matrix,entrez.data,by.x=0,by.y=0)

        row.names(all.data) <- all.data[,"EntrezId"]
        exp.data <-all.data[,3:ncol(all.data)]

    return(exp.data)

}


getSymbolFromEntrezID <- function(exp.data) {

   maps <-lookUp(rownames(exp.data),'org.Hs.eg','SYMBOL')
   entrezIds.matrix <- data.frame(sapply(maps, "[",  1))
   colnames(entrezIds.matrix) <- c("id")
   all.data=merge(entrezIds.matrix,exp.data,by.x=0,by.y=0)
   exp.data <-all.data[,2:ncol(all.data)]
   
# Different ids may map to more than 1 symbol   
   exp.data.new=exp.data[!duplicated(exp.data[,"id"]),]
# Get rid of any NA symbols   
   exp.data.new.na=exp.data.new[!is.na(exp.data.new[,"id"]),]
   rownames(exp.data.new.na)=exp.data.new.na[,"id"]
   exp.data=exp.data.new.na[,2:length(exp.data.new.na)]  
   return(exp.data)
}



