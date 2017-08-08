DESCRIPTION
Shell script runs automated analysis of copy number data from genomics arrays using either the CBS pipleine or the ASCAT pipeline:
Data can be entered at any of the following points for the CBS pipeline:
A: Raw CEL files
B: log2ratios
C: Binary codes
D: segmented 
Examples of each of these files can be found on the Examples page of the O-miner website
For the ASCAT pipeline only raw CEL files can be entered

PACKAGE DEPENDENCIES
This workflow is dependent on the ASCAT code being installed (code can be downlaoded from the ASCAT website)
library biomaRt
library copynumber 

DESCRIPTION OF COMMAND LINE INPUTS 
The line below is the command line arguments that are required to execute the script:
/usr/local/bin/R --slave --no-restore --no-save --no-readline --silent  --args  Normaldir=a.sangaralingam-CBS_TEST_UPLOAD_Normal threshold= project=a.sangaral
ingam-CBS_TEST_UPLOAD analysisMethod=CBS folder=/var/www/html/onlinetool/temp/a.sangaralingam-CBS_TEST_UPLOAD baseline=user email=a.sangaralingam@qmul.ac.uk c
df_array1=GenomeWideSNP_6 patient=20 platform=six cdf_array2=NA snp_number=15 analysisType=unpaired hap= targetFile=chip.txt type=CEL mcr=none Normal=/var/www
/html/onlinetool/temp/a.sangaralingam-CBS_TEST_UPLOAD/normal.txt  < /var/www/cgi-bin/onlinetool/version_2/CNV_shell.R
perl /var/www/cgi-bin/onlinetool/version_2/generateAnnotation.pl a.sangaralingam-CBS_TEST_UPLOAD

DESCRIPTION OF COMMAND LINE ARGUMENTS:
Normaldir - this is the name given to the folder containing the data from normal samples
threshold - this is the threshold used for calling gains and losses in ASCAT 
project - this is the name given to the analysis
analysisMethod - name given to this particular type of analysis
folder - this is the full path to the folder containing all of the files for processing that have been downloaded from the web interface
baseline - this indicates the user want to use their own files for a comparison NOTE that the only option for this workflow for vaseline is user (it is left over from the genomics workflow, where it could be either hapmap or user defined baseline)
analysisMethod = ths depends on which pipeline is to be run either ASCAT or CBS
cdf_array1 = this is the corresponding cdf annotation file to be used with teh array 
patient = this is equal to the total number of patients present
platform = this is the array platform to be used i.e. six is equal to GenomeWideSNP_6 array 
cdf_array2 = this is teh corresponding cdf annotation file to be used with the array - only chips 500k and 100K have cdf_array1 and cdf_array2 
snp_number = this is the minimum number of consecutive SNPs that should be used to define a region 
analysisType = paired/unpaired analysis
hap = 0/1 indicating whether hapmap dataset is going to be used as the baseline
targetFile = this is the targetfile to be used e.g. chip.txt
type = this is the entrypoint to the pipeline i.e. CEL,log2ratio,smoothed or segmented
mcr= this is the type of minimum common region analysis to be used
Normal - this is the .txt file that acts as a target file for the normal samples

DESCRIPTION OF WORKFLOW OUTPUT:
CBS & ASCAT 
cghweb- contains subdirectories (a) input - contains separet .tx files for each sample containing log2ratios for each probe on the array, (b) Matrix - contains one .txt file with samples as column names, probeids as rownames, (c) output  and (d) profiles - this conatins a .tct file for each sample with teh smoothed log2ratios for each sample 
Cluster - contains the .png unsupervised clustering plot using the log2ratios
FP - frequency plots 
output - smoothed log2ratios for each of the probes on teh array with samplenames as columns
QC - contains a copy of the target file, locations -the names of the directories containing the normal and tumor samples
Regions - .txt files of gain/loss regions both annotated and unannotated
threshold - threshold to call gains/losses


FUNCTIONS - CBS:
function is followed by R cript that calls it:






