############################################################################################################
#
#File: Exome_copy_number_mod.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to run an automated pipeline for the post-processing of exome sequencing data

for (e in commandArgs()) {
	ta = strsplit(e,"=",fixed=TRUE)
	ta
		temp = ta[[1]][2]
		assign(ta[[1]][1],temp)
		cat("assigned ",ta[[1]][1]," the value of |",temp,"|\n")
}

rootdir="/var/www/cgi-bin/onlinetool/version_2/"
setwd(rootdir)

target_path <- "../../../html/onlinetool/temp"


outputdir = paste("ominer_results",project,sep="/")
#Create directories for the results of analysis
dir.create(outputdir)
resultdir <- outputdir
dir.create(paste(resultdir,"/Regions/",sep=""))# when automating
dir.create(paste(resultdir,"/Regions/","segments",sep=""))# when automating
dir.create(paste(resultdir,"/Regions/","loh",sep=""))# when automating
dir.create(paste(resultdir,"/cghweb/",sep=""))
cghwebdir = paste(resultdir,"/cghweb/",sep="")
        dir.create(paste(resultdir,"/Cluster/",sep=""))
        dir.create(paste(resultdir,"/Heatmaps",sep=""))
        dir.create(paste(resultdir,"/FP",sep=""))
        dir.create(paste(resultdir,"/bed",sep=""))
        dir.create(paste(resultdir,"/MCR",sep=""))
        dir.create(paste(resultdir,"/threshold/",sep=""))# when automating
        dir.create(paste(resultdir,"/output/",sep=""))# when automating
        qcdir = paste(resultdir,"/QC",sep="")
                dir.create(qcdir)
target = paste(target_path,project,"targets.txt",sep="/")
###copy files to the correct directories to display the output correctly (.pl scripts in cgi-bin rely on these files being found in specific places)
file.copy(target,qcdir)
file.copy(target,resultdir)
annotation_file = paste(target_path,project,"annotations.out",sep="/")
file.copy(annotation_file,resultdir)
object_out_file = paste(target_path,project,"obj.out",sep="/")
file.copy(object_out_file,resultdir)
chip_file = paste(target_path,project,"chip.txt",sep="/")
###read in the target file 
pd <- read.table(chip_file, sep="\t", as.is=TRUE,strip.white=TRUE,header=T)
file.copy(chip_file,resultdir)


###assign variable to the file names of processed files 
input_files <- pd$Exome_CNV_files           
print (input_files)
rootdir="/var/www/cgi-bin/onlinetool/version_2/"
workdir <- rootdir
setwd(rootdir)

#####call R script to modify/manipulate output of varscan - (this has already been modified with process_varscan_output.pl)
source("exome_copynum/readvcf_modified.R")

####for each of the processed files 
for (i in 1:length(pd$Exome_CNV_files)) {
sample <- paste(pd$Exome_CNV_files[i],sep="")
outfile_name <- paste(pd$Name[i],sep="")


readvcf(folder,sample, outfile_name,project, outputdir)
}
  
source("exome_copynum/snp_paired_1_mod.R")
#get input samplenames from a target.txt that has been generated


vcf_files <- pd$Name
 for (i in 1:length(pd$Name)) {
sample <- paste(pd$Name[i],sep="")

snp_paired(sample,project,outputdir,workdir)
} 


source("ASCAT_output_parser_1.R")
for (i in 1:length(pd$Name)) {
sample <- paste(pd$Name[i],sep="")
ASCAT_segment_parser(sample,project)
}


#######change directory to where all the output files from ASCAT are found to create one file where they all are concatenated

chdir <- paste("ominer_results/",project,"/Regions",sep="")
setwd(chdir)
####put all of the segments files from each sample into one large file
awk_command <- "awk 'NR==1 {header=$_} FNR==1 && NR!=1 {$_ ~$header getline; }{print}' *segments_P*.txt > segments_final.out.txt"
system(awk_command)
setwd("../../../")
source("ASCAT_segments_plot_rev.R")
ASCAT_segments_plot(project)


resultdir = "/var/www/cgi-bin/onlinetool/version_2/ominer_results/"
setwd(resultdir)
zipcommand=paste("zip -r ",project,".zip ",project,sep="")
print(zipcommand)
system(zipcommand)
file.rename(paste(resultdir,project,".zip",sep=""),paste(resultdir,project,"/output/",project,".zip",sep=""))

command="/var/www/cgi-bin/onlinetool/version_2/mail.pl"
url="Your results are now available at : http://o-miner.org/cgi-bin/onlinetool/version_2/cghOutputView.pl?"
subject <- "O-miner Results"
mailcommand=paste(command," ",email," \"",subject,"\" \"",url,project,"\"",sep="")

system(mailcommand)







