<b>DESCRIPTION</b> 
Shell script runs analysis of Illumina  expression arrays - script can be used from different entry points
A. Unnormalised expression matrix
B. Normalised expression matrix

<b>PLATFORMS</b> 
This workflow can be used for the following arrays: Illumina HumanHT-12 V3, Illumina HumanHT-12 V4 and Illumina MouseRef-8 v2.0

<b>PACKAGE DEPENDENCIES</b>
This workflow is dependent on the following packages being installed:
library("simpleaffy")
library("affy")
library("lumi")
library("arrayMvout")
library("arrayqualitymetrics")
library("affyPLM")
library("annotate")
library("limma")
library("sva")
library("GOstats")

<b>DESCRIPTION OF COMMAND LINE INPUTS</b>
The line below is the command line arguments that are required to execute the script:
R --file=/data/BCI-BioInformatics/ominer/illumina_exp.R  --args   platform=ht12v3 dataType=normalised folder=/var/www/html/onlinetool/temp/a.sangaralingam-illumina_expression_array project=a.sangaralingam-illumina_expression_array
target=/var/www/html/onlinetool/temp/a.sangaralingam-illumina_expression_array/target.txt comp=/var/www/html/onlinetool/temp/a.sangaralingam-illumina_expression_array/comp.txt normalisation=rsn filter=sd filterval=40 adjust=BH
pvalue=0.05 foldchange=0.5 limmamethod=separate analysis=unpaired replicates=no combat=0 aqm=0

<b>DESCRIPTION OF COMMAND LINE ARGUMENTS</b>
platform == Illumina platform that is being analysed for example for Illumina HumanHT-12 V3 platform == "ht12v3"
dataType == Datatype of the data that is being used e.g  normalised == normalised expression matrix, unnormalised == unnormalised expression matrix
folder = Full path to to where the unnormalised or normalised data is found. Normalised or filtered data is in the format of a .txt file
project = name of your project; please note that in order for this workflow to work the target files and comparison files need to be placed in a folder with the project name
target = Full path to the target file
comp = full path to the comparison file
normalisation - normalisation method that is to be used e.g gcrma, rma etc
filter = method used to filter the normalised expression matrix e.g SD = standard deviation
filterval = % to filter the normalised expression matrix - to retain the top 40% most variable genes filterval = 40.
adjust = Method used to adjust the False Discovery (FDR) rate within limma BH = Benjamini_Hochberg method
pvalue = Pvalue to use as cutoff when listing the differentially expressed genes
foldchange = foldchange cutoff to use when listing the differentially expressed genes
limmamethod = method to use within lima when comparing the groups
analysis = whether a paired/unpaired analysis is executed
replicates = yes if there are technical replicates within your experiment
combat = This argument applies batch efface correction to your analysis 0 if batch effect correction is not to be applied and 1 if it is to be applied.
Please note that if combat is to be applied to your analysis an extra column need to be added to the target file with the column header = "study"
aqm = if arrayqualitymetrics is to be used to analyse your data 1= if it is to be used and 0 of it is not to be used.

<b>DESCRIPTION OF WORKFLOW OUTPUT</b>
Data is output to ominer_results directory a folder is automatically created with this directory with your chosen project name.
Subfolders are created within this directory and these are QC, norm, DifferentialExpression 
QC folder contains the results of running quality control metric on your data
Norm - contains the normalised expression matrix and filtered expression matrices
Differential Expression - contains the lists of differentially expressed genes, and Venn diagrams 

<b>FUNCTIONS</b>
<b>1.Quality Control</b>
illumina_qc function performs quality control metrics on Illumina data and prints out arrays that failed QC to outliers.txt using rrayMvout & lumi this is run as default
usage:runQC(dataType,dat,target_file,project)
<b>Arguments</b>
dataset = name of .txt file containing either a normalised or a unnormalised matrix
project = name of project
target_file = full path to target file
analysis - unpaired.paired analysis
normalisation= normalisation method to be used i.e. rsn,vsn,ssn,quantile,loess

<b>2. Quality control with arrayQualityMetrics </b>
QC function using arrayqualitymetrics
usage:illumina_aqm(dataset,target_file,project)
<b>Arguments</b>
dataset = name of the .txt file of the unnormalised matrix
target_file = Name of target file
project = name of project

<b>3.Minimising batch effect using the COMBAT algorithm for a meta-analysis </b>
run_combat function runs batch effect correction on data from one or more studies
usage:run_combat(target_file,norm_matrix,project)
<b>Arguments</b>
target_file = name of target file
norm_matrix = name of the normalised expression matrix
project = name of project

<b>4. Filtering </b>
Filtering function extracts the most variable genes passing a certain percentage filter from the normalised expression matrix
usage:filtering(filter,filterval,combat,project)
<b>Arguments</b>
filter = Method of filtering to be used SD = standard deviation, 
filterval =  % to filter the normalised expression matrix - to retain the top 40% most variable genes filterval = 40.
combat = 0 if combat is not to be used and 1 if combat is to be used
project = project name

<b>5. Differential Expression </b>
DifferentialExpression:Calculate the genes that are differentially expressed between groups e.g. two biological groups
usage:run_limma(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,value,foldchange,platform,project)
<b>Arguments</b>
target_file = name of target file
analysis = whether analysis is paired/unpaired
data = filtered normalised expression matrix
comp_file = text file containing the comparison between biological groups that are to be made - up to six different comparisons can be made
replicates = yes/no whether technical replicates
limmamethod = method to be used within limma e.g separate, global, nestedF
adjust = method to be used to adjust the false discovery rate (FDR) e.g BH = Benjamini-Hochberg
value = value threshold genes below this value are listed as  differentially expressed
foldchange = log fold change value threshold genes with a fold change below this value are listed as differentially expressed
platform = Illumina platform that is to be analysed e.g. for Affymetrix HGU133plus2 argument to be used is hgu133plus2
project = name of project

<b>6. Clustering </b>
Clustering:Performs non-hierarchical clustering on normalised expression matrix 
usage:cluster_all(target_file,k,project)
<b>Arguments</b>
target_file = name of target file
k = number of clusters to cluster the data
project = name of project
 

<b>7. Venn diagram generation </b>
Venn diagrams:Generates venn diagrams from lists of differentially expressed genes
usage:generate_venn("decideTestsSummary.txt",project)
<b>Arguments</b>
decideTestsSummary.txt = output from running limma to identify differentially expressed genes
project = name of project

<b>8. Heatmap generation </b>
Heatmaps: A heatmap is created for each comparison where there is a list of statistically significant probes
usage:heatmap_generate(target,project,comp_file)
<b>Arguments</b>
target - full path to target file
project - name of project
comp_file - full path to comparisons file

<b>9. Gene Ontology analysis </b>
Gene Ontology analysis of statistically important Gene Ontologies
illum_GO: Analysis of statistically importnat Gene Ontology terms
usage:illum_GO(project,comp)
<b>Arguments</b>
project - name of project
comp- full path to the comp.txt file