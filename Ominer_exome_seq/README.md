DESCRIPTION
Shell script runs automated analysis of the post-processing of exome-sequencing data to generate copy number estimates
A. Takes as input processed files from varcsan analysis and processed with process_varscan.pl 

PACKAGE DEPENDENCIES
This workflow is  dependent on the following R packages being installed:
1. biomaRt
2. copynumber 

DESCRIPTION OF COMMAND LINE INPUTS
The line below is the command line arguments that are required to execute the script:
R --slave --no-restore --no-save --no-readline --silent  --args  Normaldir=a.sangaralingam-EXOME_EXAMPLE_ANALYSIS_Normal threshold=0.6 project=a.sangaralingam-EXOME_EX
AMPLE_ANALYSIS analysisMethod=EXOME_CNV folder=/var/www/html/onlinetool/temp/a.sangaralingam-EXOME_EXAMPLE_ANALYSIS baseline=user email=a.sangaralingam@qmul.ac.uk cdf_array1=0 patien
t=20 platform=0 cdf_array2=NA snp_number=15 analysisType=unpaired hap= targetFile=chip.txt type=TXT mcr=none Normal=/var/www/html/onlinetool/temp/a.sangaralingam-EXOME_EXAMPLE_ANALYS
IS/normal.txt  < /var/www/cgi-bin/onlinetool/version_2/Exome_copy_number_mod.R

DESCRIPTION OF THE COMMAND LINE ARGUMENTS:
Normaldir - this is the name given to the folder containing the data from normal samples
threshold - this is the threshold used for calling gains and losses in ASCAT 
project - this is the name given to the analysis
analysisMethod - name given to this particular type of analysis
folder - this is the full path to the folder containing all of the files for processing that have been downloaded from the web interface
baseline - this indicates the user want to use their own files for a comparison NOTE that the only option for this workflow for vaseline is user (it is left over from the genomics workflow, where it could be either hapmap or user defined baseline)
cdf_array = this argument is needed only for an analysis involving genomics array data
patient = this is equal to the total number of patients present

DESCRIPTION OF WORKFLOW OUTPUT
Data is output to ominer_results directory where an output folder is automatically created with teh project name
subfolders created in this directory are:
QC
REGIONS - this contains the .txt and .xls files for each sample of gains and losses both annotated and unannotated
CLUSTER - contains aberration plots in both .png and .pdf format
FP - Contains all of the output plots from ASCAT - also contains frequency plots
NOTE: some of the folders here are empty but the same format as the genomics pieline has had to be followed as the same output .pl script was used and this output script relies on theat the same subfolders are present.

The running order of each of the functions required for this workflow are listed below - after each of the functions the relevant R script is listed.
FUNCTIONS:

1. readvcf - readvcf_modified.R
2. snp_paired - snp_paired_1_mod.R
3. ASCAT_segment_parser - ASCAT_output_parser_1.R
4. ASCAT_segments_plot - ASCAT_segments_plot_rev.R
