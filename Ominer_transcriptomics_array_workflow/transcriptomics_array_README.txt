{\rtf1\ansi\ansicpg1252\cocoartf1038\cocoasubrtf360
{\fonttbl\f0\froman\fcharset0 TimesNewRomanPSMT;\f1\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww29460\viewh17420\viewkind0
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\ql\qnatural\pardirnatural

\f0\b\fs33\fsmilli16920 \cf0 O-miner:Affymetrix Expression analysis workflow\

\fs30 Description: 
\b0 Shell script runs analysis of Affymetrix expression arrays - script can be used from different entry points A. Raw CEL files B. Normalised expression matrix\

\b Instructions: 
\b0 To run this pipeline you must have two separate text files, a target file and a comparisons file. Samples files of these can be found in this directory target.txt and comp.txt Code needed to run this pipeline and sample files are found on midplus /scratch/ominer/transcriptomics/Affymetrix\
Example.sh script to run an unpaired analysis from raw CEL files \'96 affy_exp_v2.sh ./affy_exp_v2.sh to execute pipeline Example.sh script to run an unpaired analysis from normalized data matrix \'96 ./affy_exp_normalised_v2.sh to execute pipeline\

\b Example target files:\

\b0 Format of the target files for different types of analysis are as follows: \
A. Unpaired analysis no replicates \'96 target.txt \
B. Unpaired analysis using COMBAT \'96 target_combat.txt \
C. Paired analysis no replicates \'96 target_paired.txt\
D. Unpaired analysis with replicates \'96 target_replicates.txt \
E. Survival analysis \'96 target_survival.txt (note column Surv_Period is in years)\

\b Package dependencies:\

\b0 This workflow is dependent on the following packages being installed: \
library("simpleaffy")\
 library("affy") \
library("arrayMvout")\
library("arrayqualitymetrics")\
 library("affyPLM") \
library("annotate") \
library("limma") \
library("sva") \
####This package is required if running combat library("hguplus133cdf") \
####This is the annotation package that is specific for the platform Affymetrix Human Genome U133 plus 2.0 array if you are analyzing data from another platform, please download the annotation file specific for that platform. \
library(\'93estimate\'94) ####This package is needed if running the ESTIMATE algorithm please follow the	instructions:
\f1 	\
 
\f0 http://bioinformatics.mdanderson.org/main/ESTIMATE:Overview#Installation\
library(\'93FC14.plotting.lib\'94) ####required for survival analysis library(\'93gplots\'94) ###required for survival analysis and heatmaps\

\f1 	\
1.
\f0\b Code dependencies:\

\b0 Venn.R Code.R bcctb.utils.R These can be copied from the Affymetrix directory\

\b NOTE: Before running the pipeline you will need to create a directory named \'93ominer_results\'94\
Description of command line inputs required:\

\b0 The line below is the command line arguments that are required to execute the script:\
R--file-/scratch/ominer/transcriptomics/Affymetrix/affy_expression.R --args platform=hgu133plus2 dataType=CEL folder=/scratch/ominer/transcriptomics/Affymetrix/CEL project=a.sangaralingam- AFFY_TEST target=/scratch/ominer/transcriptomics/Affymetrix/target.txt comp=/scratch/ominer/transcriptomics/Affymetrix/comp.txt normalization=rma filter=sd filterval=40 adjust=BH pvalue=0.05 foldchange=2 limmamethod=separate analysis=unpaired replicates=no combat=0 survival=0 estimate=0 aqm=0 gostats=1\

\b Description of command line arguments: platform 
\b0 - Affymetrix platform that is being analysed for example for AffymetrixU133plus2 platform - "hgu133plus2" 
\b dataType 
\b0 - Datatype of the data that is being used e.g CEL -- raw CEL files\
normalised - normalised expression matrix \

\b filtered 
\b0 - filtered expression matrix \

\b folder 
\b0 - Full path to to where the raw data is kept or normalised or filtered data. Note if using normalised or filtered data this is in the format of a .txt file \

\b project 
\b0 - name of your project; please note that in order for this workflow to work the target files and comparison files need to be placed in the same folder as your shell script command. The target file should be named - target.txt and the comparison file should be named comp.txt\

\b comp 
\b0 - full path to the comparison file \

\b normalisation 
\b0 - normalisation method that is to be used e.g gcrma, rma etc \

\b filter 
\b0 - method used to filter the normalised expression matrix e.g SD - standard deviation \

\b filterval 
\b0 - % to filter the normalised expression matrix - to retain the top 40% most variable genes filterval - 40. \

\b adjust 
\b0 -Method used to adjust the False Discovery (FDR) rate within limma BH - Benjamini_Hochberg method \

\b pvalue 
\b0 - Pvalue to use as cutoff when listing the differentially expressed genes \

\b foldchange 
\b0 - foldchange cutoff to use when listing the differentially expressed genes \

\b limmamethod 
\b0 - method to use within lima when comparing the groups\

\b analysis 
\b0 - whether a paired/unpaired analysis is executed \

\b replicates 
\b0 - yes if there are technical replicates within your experiment \

\b combat 
\b0 - This argument applies batch effect correction to your analysis 0 if batch effect correction is not to be applied and 1 if it is to be applied. Please note that if combat is to be applied to your analysis an extra column need to be added to the target file with the column header - "study" 
\b aqm 
\b0 - if arrayqualitymetrics is to be used to analyse your data 1- if it is to be used and 0 of it is not to be used.\

\b Description of output:\

\b0 Data is output to ominer_results directory a folder is automatically created with this directory with your chosen project name. Subfolders are created within this directory and these are QC, norm, DifferentialExpression QC - contains the results of running quality control metric on your data\
Norm - contains the normalised expression matrix and filtered expression matrices Differential Expression - contains the lists of differentially expressed genes, and Venn diagrams\

\b Functions: runQC function 
\b0 performs quality control metrics on Affymetrix data and prints out arrays that failed QC to outliers.txt usage:runQC(dataType,dat,target_file,project) 
\b Arguments: 
\b0 dataType - CEL,normalised dat - project - name of project target_file - Name of target file\

\b Normalisation function	
\b0 performs	normalisation on expression matrix 
\b usage
\b0 :Normalisation(normalisation,target_file,dat,project)\

\b Arguments:\

\b0 normalisation - normalisation method choose from gcrma,rma target_file - Name of target file dat - name of normalized matrix (this is normalized.txt) and is hardcoded into affy_expression.R\
project - name of project\

\b run_combat function 
\b0 runs batch effect correction on data from one or more studies 
\b usage
\b0 :run_combat(target_file,norm_matrix,project) 
\b Arguments: 
\b0 target_file - name of target file\
norm_matrix - name of the normalised expression matrix project - name of project\

\b Filtering: Filtering function 
\b0 extracts the most variable genes passing a certain percentage filter from the normalised expression matrix usage:filtering(filter,filterval,combat,project) 
\b Arguments: 
\b0 filter - Method of filtering to be used SD - standard deviation, filterval - % to filter the normalised expression matrix - to retain the top 40% most variable genes filterval - 40. combat - 0 if combat is not to be used and 1 if combat is to be used project - project name\

\b Differential Expression:\

\b0 Calculate the genes that are differentially expressed between groups e.g. two biological groups 
\b usage:run_limma
\b0 (target_file,analysis,data,comp_file,replicates,limmamethod,adjust,val ue,foldchange,platform,project) Arguments:\
target_file - name of target file analysis - whether analysis is paired/unpaired data - filtered normalised expression matrix comp_file - text file containing the comparison between biological groups that are to be made replicates - yes/no whether technical replicates limmamethod - method to be used within lima e.g separate, global, nestedF adjust - method to be used to adjust the false discovery rate (FDR) e.g BH - Benjamini- Hochberg value - value threshold genes below this value are listed as differentially expressed foldchange - log fold change value threshold genes with a fold change below this value are listed as differentially expressed platform - Affymetrix platform that is to be analysed e.g for Affymetrix HGU133plus2 argument to be used is hgu133plus2 project - name of project \

\b Clustering: 
\b0 Performs non-hierarchical clustering on normalised expression matrix 
\b usage
\b0 :cluster_all(target_file,n,project) 
\b Arguments: 
\b0 target_file - name of target file n - number of clusters to cluster the data into project - name of project\

\b Venn diagrams:\

\b0 Generates	venn	diagrams	from	lists	of	differentially expressed genes 
\b usage
\b0 :generate_venn("decideTestsSummary.txt",project) 
\b Arguments: 
\b0 decideTestsSummary.txt - output from running limma to identify differentially expressed genes project - name of project\

\b Gene Ontology\

\b0 Runs hypergeometric tests on Gene ontology terms 
\b Usage: 
\b0 gostats_analysis(platform,norm_matrix,project,comp_file) 
\b Arguments: 
\b0 platform- name of platform used norm_matrix - \'93normalised.txt\'94 project - name of project comp_file - path to a .txt file of comparisons\

\b Survival analysis\

\b0 Generates Kaplan-Meier survival plots (for 5,10 and 15 years) 
\b Usage: 
\b0 survival_os_plot(target_file,project) 
\b Arguments: 
\b0 target_file - path to target file (target file needs to be in a certain format for survival analysis (file survival_target is a template for this)\
project - name of project\

\b Estimation of tumor purity\

\b0 Runs ESTIMATE algorithm on Affymetrix expression data to generate tumour purity estimates 
\b Usage: 
\b0 run_estimate(project,platform,norm_matrix,target_file) 
\b Arguments:\

\b0 Project- name of project Platform - platform used norm_matrix - \'93filtered_data.txt\'94 target_file - path to target file\

\b Heatmap generation\

\b0 Generates a heatmap using the differentially expressed genes from a normalized expression matrix 
\b Usage:
\b0 heatmap_generate(targets,de_genes,data,project) 
\b Arguments:\

\b0 Targets - name/path to target file de_genes - name of the list of differentially expressed to use within a heatmap data- normalized matrix project - name of project\

\f1 	\
 \'a0 5	}