############################################################################################################
#
#File filter_res.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to filter log2ratios according to the minimum number of SNPS 
#Input a text file of log2ratios - columns are sample names, probeids, chromosome and position are also present as column headers. NOTE a file of filtered log2ratios can also be used as input 
#Sample input file R_input.txt
#inputs to function are: (1) snp_number - minimum number of SNPs, (2) Cancer - name of project
filterBinary <- function (snp_number,Cancer) 
{
  path_to_change <- paste("ominer_results/",Cancer,"/output/",sep="")	
setwd(path_to_change)
  
    first.data <- read.table("R_input.txt", header = T, sep = "\t", as.is = T)
    results <- read.table("results.txt", header = T, sep = "\t", as.is = T)
SNP_nb = 15;
	print("filter binary")
	print(SNP_nb)
	results.filtered = results;
	a=colnames(results[,4:ncol(results)])
	for(j in a)
	{
        x=results[,j]
		y=first.data[,j]
		z.rle <- rle(x)
		
		ends <- cumsum(z.rle$lengths)
		starts <- ends - z.rle$lengths + 1
		
		indexes <- with(z.rle, data.frame(starts, ends, lengths, values))
		
#what you want to replace:
		s=starts[z.rle$lengths<snp_number &z.rle$values!=0]
		e=ends[z.rle$lengths<snp_number &z.rle$values!=0]
		
		for (i in 1:length(s))
		{
			x[s[i]:e[i]]=0
			y[s[i]:e[i]]=0
			
		}

		results.filtered[,j]=x
		first.data[,j]=y
	}
        write.table(first.data,"R_input_filtered.txt",sep="\t",row.names=FALSE,quote=FALSE)
        }