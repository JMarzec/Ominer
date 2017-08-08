DESCRIPTION
Shell script runs automated analysis of methylation data from either raw idat files or unoormalised/normalised data
Data can be entered at any of the following points:
A: Raw IDAT files
B: unnormalised
C: normalised
PACKAGE DEPENDENCIES
This workflow is dependent on the following packages being installed:
cHamp
limma
affy
venn 
Gostats 

DESCRIPTION OF COMMAND LINE INPUTS 
The line below is the command line arguments that are required to execute the script:
/var/www/cgi-bin/onlinetool/version_2/usr/bin/R  --slave --no-restore --no-save --no-readline --silent  --args  email=a.sangaralingam@qmul.ac.uk platform=450k dataType=raw folder
=/var/www/html/onlinetool/temp/a.sangaralingam-METHYALTION_INT_TEST project=a.sangaralingam-METHYALTION_INT_TEST target=/var/www/html/onlinetool/temp/a.sangaralingam-METHYALTION_INT_
TEST/target.txt csv=/var/www/html/onlinetool/temp/a.sangaralingam-METHYALTION_INT_TEST/target.csv comp=/var/www/html/onlinetool/temp/a.sangaralingam-METHYALTION_INT_TEST/comp.txt nor
malisation=0 filter=sd filterval=40 adjust=BH pvalue=0.05 foldchange=2.0 diffmethod=limma limmamethod=separate analysis=unpaired replicates=no combat=0 champ=1 < /var/www/cgi-bin/onl
inetool/version_2/ominer_methylation_pipeline.R 
DESCRIPTION OF COMMAND LINE ARGUMENTS:
platform = 450k or 27k
dataType = raw, normalised/unnormalised
folder = full path to teh data uploaded rom the site
target = full path to target file uploaded from interface
comp = full path to comparisons file
normalisation = normalisationmethod
filetr = method used to flter datae.g. iqr, SD etc
filterval - 40, 10 etc
adjust - BH, BF etc
pvalue = value to sort results e.g 0.05
foldchange = value to sort results
diffmethod = limma
limmamethod = separate etc
analysis - paired/unpaired
replicates = yes/no
combat = 0/1 - 1 if it is used and 0 for no
champ = 1/0 = 1 if it is sued an 0 is it is not 

Functions:
these are called in the following order - function name followed by code
run_limma_methy.R = run_limma
cluster_all - cluster_code_methy.R
generate_venn - generate_venn.R
gostats_analysis - methylation_GO_analysis.R 
each function is annotated with inputs etc


OUTPUT OF WORKFLOW:
cluster - cluster_plot
DifferentialExpression - results of DE analysis
GO - statistically significant GO terms
norm - normalised/unnormalised data
QC - results of QC 

