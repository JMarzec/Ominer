############################################################################################################
#
#File name get_regions_mod_annotated.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: R script to write regions gained and lost to separate .txt files according to the minimum number of SNPs 
###Arguments (1) snp_number - this is the minimum number of SNPs and (2) type - this is the type of data e.g. log2ratio, smoothed etc
getRegions <- function (snp_number, type) {
	path_to_change <- paste("ominer_results/",Cancer,"/output/",sep="")
	setwd(path_to_change)
    print("Calculating GLAD regions")
    if (type != "smoothed") {
    first.data <- read.table("R_input.txt", header = T, sep = "\t", 
        as.is = T)
        }
        
    regions <- read.table("binary_coded_filtered.txt", header = T, sep = "\t", as.is = T)
    results <- read.table("results.txt", header = T, sep = "\t", as.is = T)




a=colnames(regions[,4:ncol(results)])
	print(a)
	#setwd(paste(resultdir,Cancer,"/Regions/",sep=""))
	
	for(j in a)
	{
		print(j)
		x=regions[,j]
		if(type != "smoothed") {
			y=first.data[,j]
		}
		z.rle <- rle(x)
		
		ends <- cumsum(z.rle$lengths)
		starts <- ends - z.rle$lengths + 1
		
		indexes <- with(z.rle, data.frame(starts, ends, lengths, values))	
		
#gains
		s=starts[z.rle$lengths>=snp_number &z.rle$values==1]
		e=ends[z.rle$lengths>=snp_number &z.rle$values==1]
		
#
		gain=NULL;
		nbre=NULL;
		log_ratio=NULL;
		Chromosomes=NULL;
		probes_start=NULL;
		probes_end=NULL;
		pos_start=NULL;
		pos_end=NULL;
		patient=NULL;
		print(length(s))
		
		if(length(s) >0) {
			for (i in 1:length(s))
			{
				nbre[i]= e[i]-s[i]+1
				
				patient[i]= j;
				if(type == "smoothed") {
					log_ratio[i]="NA"
				} else {
				    log_ratio[i]= median(y[s[i]:e[i]])
				}
				Chromosomes[i] = results[s[i],c("Chromosome")]
				
				pos_start[i]=results[s[i],c("Position")]
				pos_end[i]=results[e[i],c("Position")]
			}
						 CNA = "gain"
                NOTE = "   "
                gain = data.frame(patient, nbre, Chromosomes,pos_start,pos_end,log_ratio,CNA,NOTE)
            write.table(gain, file = paste("../Regions/gains", j, ".txt", 
                sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

		}
#loss
			s=starts[z.rle$lengths>=snp_number &z.rle$values==-1]
			e=ends[z.rle$lengths>=snp_number &z.rle$values==-1]
			
			loss=NULL;
			nbre=NULL;
			log_ratio=NULL;
			Chromosomes=NULL;
			probes_start=NULL;
			probes_end=NULL;
			pos_start=NULL;
			pos_end=NULL;
			patient=NULL;
			if(length(s) >0) {
			for (i in 1:length(s))
			{
				nbre[i]= e[i]-s[i]+1
				patient[i]= j;
				if(type == "smoothed") {
					log_ratio[i]="NA"
				} else {
				    log_ratio[i]= median(y[s[i]:e[i]])
				}
				Chromosomes[i] = results[s[i],c("Chromosome")]
				pos_start[i]=results[s[i],c("Position")]
				pos_end[i]=results[e[i],c("Position")]
			}
			 CNA = "loss"
            NOTE = "   "
            loss = data.frame(patient, nbre, Chromosomes, pos_start,pos_end, log_ratio,CNA,NOTE)
            write.table(loss, file = paste("../Regions/losses", j, ".txt", 
                sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

			
			
		}
	}
	}
