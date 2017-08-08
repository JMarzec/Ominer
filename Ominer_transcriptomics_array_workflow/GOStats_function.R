##########################################################################################
#
#File name: GOStats_function.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#
#
#
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
##########################################################################################
##########################################################################################
#Description: function to perform GOstats analysis on the results of differential expression analysis 
#platform - platform from which data analysed is from e.g. for hgu133plus2 platform = hgu133plus2
#norm - "normalised.txt'
#project - name given to the analysis
#comp - full path to the comparisons file

################################################################################################KM plots - overall survival from a target file generates .txt file and KM plots for 5,10 and 15 years
######







#####GOstats function uses as an input the affyprobeids of the probes that pass the logfc and pvalue cutoffs
######takes as input platform
#####normalisation method that is used
#######project name
#####path to comp file - can get this from the master script i.e.

gostats_analysis <- function(platform,norm,project,comp) {
library("annotate")
library(paste(platform,".db",sep=""),character.only = TRUE)
#library("hgu133plus2.db")
library("GOstats")
library("xtable")
####Filter out records that do not have any annotation
#norm_file <- paste(norm,".exp",sep="")
#A.data <- read.table(norm_file,row.names = 1, header = T, sep = "\t") 
comps <- read.table(comp, sep = "\t", as.is = T, header = F, 
        strip.white = T)

A.data <- read.table(paste("ominer_results/",project,"/","norm","/","normalised.txt",sep=""), sep = "\t", as.is = T, header = T, 
        strip.white = TRUE, row.names = 1)
entrezIds <-mget(rownames(A.data),envir=get(paste(platform,"ENTREZID",sep="")))
                
                haveEntrezId<-names(entrezIds)[sapply(entrezIds,function(x) !is.na(x))]
                
                entrez.data <-A.data[haveEntrezId,]
                
### Filter out records with no GO annotation
                haveGo <- sapply(mget(rownames(entrez.data),get(paste(platform,"GO",sep="")))
                                                 , function(x) {
                                                 if(length(x) == 1 && is.na(x))
                                                 FALSE
                                                 else TRUE
                                                 })
                numNoGO <-sum(!haveGo)
                go.data <-entrez.data[haveGo,]


## map each probe set to exactly one Entrez gene id.  If multiple probes found select probe with largest IQR
                library("genefilter")
                data.iqr<-apply(go.data,1,IQR)
                uniqGenes<-findLargest(rownames(go.data),data.iqr,platform)
                go.data<-go.data[uniqGenes,]
                chipAffyUniverse<-rownames(go.data)
                chipEntrezUniverse<-mget(chipAffyUniverse,env=get(paste(platform,"ENTREZID",sep="")))
                chipEntrezUniverse<-unique(unlist(chipEntrezUniverse))
                        
        for (i in 1:length(comps$V1)) {
                print(comps$V1[i])
        fg <- strsplit(comps[,1],"=")
        ci <- fg[[i]][1]
        print ("This is the comparison")
        print(ci)
####need to do this for all comparisons
       filteredx <- read.table(paste("ominer_results/",project,"/","DifferentialExpression","/",ci,"annotated.txt",sep=""),row.names = 1, header = T, sep = "\t")
         #filteredx <- read.table(paste(ci,"_annotated.txt",sep=""),row.names = 1, header = T, sep = "\t")
probeid=rownames(filteredx)
#platform = "hgu133plus2"
####collect entrez IDs for all of the affy_probeids passing the FC and pval cutoffs
x_entrezIds <-mget(probeid,envir=get(paste(platform,"ENTREZID",sep="")))
x_haveEntrezId<-names(x_entrezIds)[sapply(x_entrezIds,function(x) !is.na(x))]                   
x_entrez.data <-x_entrezIds[x_haveEntrezId]
                                              
haveGo <- sapply(mget(names(x_entrez.data),get(paste(platform,"GO",sep="")))
                                                         , function(x) {
                                                         if(length(x) == 1 && is.na(x))
                                                         FALSE
                                                         else TRUE
                                                         })
                        
                        x_go <-x_entrez.data[haveGo]
                        affyUniverse<-names(x_go)
                        entrezUniverse<-unlist(mget(affyUniverse,envir=get(paste(platform,"ENTREZID",sep=""))))
                        entrezUniverse<-unique(unlist(entrezUniverse))
                        goTypes=c("BP","CC","MF")
                        
                        for(go in (goTypes))
                        {
                                print(go)
                                params<-new("GOHyperGParams",geneIds=entrezUniverse,universeGeneIds=chipEntrezUniverse,
                                annotation=paste(platform,".db",sep=""),ontology=go,pvalueCutoff=0.005,conditional=FALSE,testDirection="over")
                                hgOver<-hyperGTest(params)
                                ##hgOver
                                
                                ##file.create(paste("ominer_results/",project,"/GO/",ci,go,"_over.html",sep=""))
                                ##htmlReport(hgOver,file=paste("ominer_results/",project,"/GO/",ci,go,"_over.html",sep=""))
                                df <- summary(hgOver)
                                row <- nrow(df)
                                print (row)
                                if (row == 0) next
                                for (i in 1:length(df[,7])){
                                        df[i,7] <- paste('<a href="http://amigo.geneontology.org/amigo/term/',df[i,1],'" target="blank">',df[i,7],'</a>',sep='')
                                }
                                print(xtable(df, caption=paste("Gene to GO ",go," test for over representation"),align="llrrrrrl",digits=3,display=c("s","s","f","f","f","d","d","s")),type="html",file=paste("ominer_results/",project,"/GO/",ci,go,"_over.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)   
                                #print(xtable(df, caption=paste("Gene to GO ",go," test for over representation"),align="llrrrrrl",digits=3,display=c("s","s","f","f","f","d","d","s")),type="html",file=paste("../",ci,go,"_over.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)   

                                #htmlReport(hgOver,file="comp1_go_over.html")
                                params<-new("GOHyperGParams",geneIds=entrezUniverse,universeGeneIds=chipEntrezUniverse,
                                annotation=paste(platform,".db",sep=""),ontology=go,pvalueCutoff=0.005,conditional=FALSE,testDirection="under")
                                hgOver<-hyperGTest(params)
                                #hgOver
                                #file.create("GO_under.html")
                                ##file.create(paste("ominer_results/",project,"/GO/",ci,go,"_under.html",sep=""))
                                #htmlReport(hgOver,file="GO_under.html")
                                ##htmlReport(hgOver,file=paste("ominer_results/",project,"/GO/",ci,go,"_under.html",sep=""))
                                df <- summary(hgOver)
                                row <- nrow(df)
                                print (row)
                                if (row == 0) next
                                for (i in 1:length(df[,7])){
                                        df[i,7] <- paste('<a href="http://amigo.geneontology.org/amigo/term/',df[i,1],'" target="blank">',df[i,7],'</a>',sep='')
                                }
                               print(xtable(df, caption=paste("Gene to GO ",go," test for under representation"),align="llrrrrrl",digits=3,display=c("s","s","f","f","f","d","d","s")),type="html",file=paste("ominer_results/",project,"/GO/",ci,go,"_under.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)   
                                 #print(xtable(df, caption=paste("Gene to GO ",go," test for under representation"),align="llrrrrrl",digits=3,display=c("s","s","f","f","f","d","d","s")),type="html",file=paste("../",ci,go,"_under.html",sep=""),caption.placement="top",include.rownames=FALSE,sanitize.text.function = force)   
                                
                        }
                       }
                       }