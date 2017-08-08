<b>DESCRIPTION</b> 
Shell script runs analysis of Affymetrix exon arraya - script can be used from 
A. Raw CEL files


<b>PLATFORMS</b> 
This workflow can be used for the following arrays: 

<b>PACKAGE DEPENDENCIES</b>
This workflow is dependent on the following packages being installed:
library("aroma.affymetrix")
library("affy")
library("biomaRt")
library("limma")
library("survival")
library("sva")
library("gostats")




<b>DESCRIPTION OF COMMAND LINE INPUTS</b>
The line below is the command line arguments that are required to execute the script:

R--file-/scratch/ominer/transcriptomics/Affymetrix/EA_run.R platform=st1 dataType=CEL folder=/var/www/html/onlinetool/temp/a.sangaralingam-TEST_EXON_ARRAY project=a.sangaralingam-TEST_EXON_ARRAY target=/var/www/html/onlinetool/temp/a.sangaralingam-TEST_EXON_ARRAY/target.txt
comp=/var/www/html/onlinetool/temp/a.sangaralingam-TEST_EXON_ARRAY/comp.txt normalisation=0 filter=sd filterval=40 adjust=BH pvalue=0.05 foldchange=0.5 limmamethod=separate analysis=unpaired replicates=no combat=0 

<b>DESCRIPTION OF COMMAND LINE ARGUMENTS</b>
platform == Affymetrix platform that is being analysed for example for Illumina HumanHT-12 V3 platform == "ht12v3"
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


<b>DESCRIPTION OF WORKFLOW OUTPUT</b>
Data is output to ominer_results directory a folder is automatically created with this directory with your chosen project name.
Subfolders are created within this directory and these are: tarnscript, exon and splicing
Within each of these folders the following are created as subfolders:
cluster folder contains the results of unsupervised hierarchicial clustering on the normalised matrix either exon, transcript or splicing
Norm - contains the normalised expression matrix and filtered expression matrices
Differential Expression - contains the lists of differentially expressed genes, and Venn diagrams 

<b>FUNCTIONS</b>
<b>1.Reading in CEL files</b>
Read in Affymetrix CEL files and create an affybatch object
usage:exon_qc(platform,project,aromadir,mydir)
<b>Arguments</b>
platform = platform to be analysed i.e. st1
project = name given to project
aromadir = full path to aroma.affymetrix directory
mydir = full path to directory where master R script is found



<b>2.Minimising batch effect using the COMBAT algorithm for a meta-analysis </b>
run_combat function runs batch effect correction on data from one or more studies
usage:run_combat(target_file,norm_matrix,project)
<b>Arguments</b>
target_file = name of target file
norm_matrix = name of the normalised expression matrix
project = name of project

<b>3. Differential Expression for exons </b>
DifferentialExpression:Calculate the exons that are differentially expressed between groups e.g. two biological groups
usage:run_limma_exon(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,value,foldchange)
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

<b>4. Differential Expression for transcripts </b>
DifferentialExpression:Calculate the transcripts that are differentially expressed between groups e.g. two biological groups
usage:run_limma_exon(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,value,foldchange)
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

<b> 5. Differential splicing score i.e. FIRMA scores to identify splice variants </b>
DifferentialExpression:Calculate the transcripts that are differentially expressed between groups e.g. two biological groups
usage:run_limma_firma(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,value,foldchange)
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

<b>6. Clustering for exons </b>
Clustering:Performs non-hierarchical clustering on normalised expression matrix of exons
usage:cluster_all_Exon(target_file,k,project)
<b>Arguments</b>
target_file = name of target file
k = number of clusters to cluster the data
project = name of project
 
<b> 7. Clustering for transcripts </b>
Clustering:Performs non-hierarchical clustering on normalised expression matrix of transcripts
usage:cluster_all_tr(target_file,k,project)
<b>Arguments</b>
target_file = name of target file
k = number of clusters to cluster the data
project = name of project

<b>8. Venn diagram generation </b>
Venn diagrams:Generates venn diagrams from lists of differentially expressed transcripts
usage:generate_venn("decideTestsSummary.txt",project)
<b>Arguments</b>
decideTestsSummary.txt = output from running limma to identify differentially expressed transcripts
project = name of project


<b>9. Gene Ontology analysis </b>
Gene Ontology analysis of statistically important Gene Ontologies from results of differentially expressed exons
Gene Ontology analysis : Analysis of statistically important Gene Ontology terms
usage:Exon_GO(project,comp)
<b>Arguments</b>
project - name of project
comp- full path to the comp.txt file