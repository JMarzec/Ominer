<b>DESCRIPTION</b> 
Shell script runs analysis of Affymetrix  expression arrays - script can be used from different entry points
A. Raw CEL files
B. Normalised expression matrix

<b>PLATFORMS</b> 
This workflow can be used for the following arrays: Affymetrix GeneChip Human Genome Array U133 Plus 2.0, GeneChip Human Genome Array U133,
GeneChip Human Genome array U95 and GeneChip Mouse Genome 430 2.0

<b>PACKAGE DEPENDENCIES</b>
This workflow is dependent on the following packages being installed:
library("simpleaffy")
library("affy")
library("arrayMvout")
library("arrayqualitymetrics")
library("affyPLM")
library("annotate")
library("limma")
library("sva")
library("hgu133plus2.db") (Annotation database relevant to the data from the array platform that is being analysed)
library("GOstats")

<b>DESCRIPTION OF COMMAND LINE INPUTS</b>
The line below is the command line arguments that are required to execute the script:
 R --file=/scratch/ominer/transcriptomics/Affymetrix/affy_expression.R --args  platform=hgu133plus2 dataType=CEL folder=/scratch/ominer/transcriptomics/Affymetrix/CEL project=a.sangaralingam-affy_expression_test target=/scratch/ominer/transcriptomics/Affymetrix/target.txt comp=/scratch/ominer/transcriptomics/Affymetrix/comp.txt  normalisation=rma filter=sd filterval=40 adjust=BH pvalue=0.05 foldchange=0.5 limmamethod=separate analysis=unpaired replicates=no combat=0 aqm=0 

<b>DESCRIPTION OF COMMAND LINE ARGUMENTS</b>
platform == Affymetrix platform that is being analysed for example for AffymtetrixU133plus2 platform == "hgu133plus2"
dataType == Datatype of the data that is being used e.g CEL == raw CEL files, normalised == normalised expression matrix, filtered == filtered expression matrix
folder = Full path to to where the raw data is kept or normalised or filtered data. Note if using normalised or filtered data this is in the format of a .txt file
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
runQC function performs quality control metrics on Affymetrix data and prints out arrays that failed QC to outliers.txt 
usage:runQC(dataType,dat,target_file,project)
<b>Arguments</b>
dataType = CEL,normalised
dat = Affymetrix object created from raw CEL files
project = name of project
target_file = Name of target file

<b>2. Normalisation </b>
Normalisation function performs normalisation on expression matrix 
usage:Normalisation(normalisation,target_file,dat,project)
<b>Arguments</b>
normalisation = normalisation method choose from gcrma,rma
target_file = Name of target file
dat = Affymetrix object created from raw CEL files
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
platform = Affymetrix platform that is to be analysed e.g. for Affymetrix HGU133plus2 argument to be used is hgu133plus2
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