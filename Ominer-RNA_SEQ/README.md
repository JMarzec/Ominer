<b>DESCRIPTION</b>
Shell script runs analysis of read counts data from RNA-Seq experiments - script can be used from different entry points
A. Normalised matrix of read counts data
B. Unnormalised matrix of read counts data 

<b>TARGET FILE FORMAT</b>
For the format of the target files, please follow the same format as those target files prepared/present in the Transcriptomics_array_workflow. e.g target_combat.txt, target_paired.txt,
target_replicates.txt and target_survival.txt 
<b>PACKAGE DEPENDENCIES</b>
This workflow is dependent on the following packages being installed:
library("edgeR")
library("goseq")
library("xtable")
library("affy")
library("limma")
library("sva")
library("biomart") 
library("GO.db")

<b>DESCRIPTION OF COMMAND LINE INPUTS:</b>
The line below is the command line arguments that are required to execute the script:
 /usr/bin/R --slave --no-restore --no-save --no-readline --silent  --args  folder=/data/BCI-BioIn
formatics/ominer/RNA_seq/post_processing  project=test target=target.txt count_table=normalised.txt combat=0 dataType=
normalised diffmethod=limma analysis=unpaired replicates=no filterval=40 filter=iqr foldchange=0.5 limmamethod=separate com
p=RNA_seq_comp.txt adjust=BH pvalue=0.05 gostats=0 normalisation=rma  survival=0  < run_RNA_SEQ.R 


<b>DESCRIPTION OF COMMAND LINE ARGUMENTS</b>
folder = Full path to to where the raw data is kept or normalised or filtered data. Note if using normalised or filtered data this is in the format of a .txt file
project = name of your project; please note that in order for this workflow to work the target files and comparison files need to be placed in a folder with the project name
target = Full path to the target file
count_table = name of the normalised or unnormalised .txt file of read counts data
combat = This argument applies batch efface correction to your analysis 0 if batch effect correction is not to be applied and 1 if it is to be applied.
dataType = normalised (for normalised read counts data) and unnormalised for unnormalised read counts data
diffmethod = limma (to use the LIMMA R package), edge (to use the edgeR package)
analysis = whether a paired/unpaired analysis is executed
replicates = yes if there are technical replicates within your experiment
filterval = % to filter the normalised expression matrix - to retain the top 40% most variable genes filterval = 40.
filter = method used to filter the normalised expression matrix e.g SD = standard deviation
foldchange = foldchange cutoff to use when listing the differentially expressed genes (this is typically set to log 2.0)
limmamethod = method to use within lima when comparing the groups, e.g, separate
comp = full path to the comparison file
adjust = Method used to adjust the False Discovery (FDR) rate within limma BH = Benjamini_Hochberg method
pvalue = Pvalue to use as cutoff when listing the differentially expressed genes
gostats = 0 if users do not want to perform Gene Ontology analysis and 1 if they do
normalisation - normalisation method that is to be used e.g gcrma, rma etc
survival - 0 if user does not wish to run survival analysis and 1 if they do
Please note that if combat is to be applied to your analysis an extra column need to be added to the target file with the column header = "study"

<b>DESCRIPTION OF OUTPUT</b>
Data is output to ominer_results directory a folder is automatically created with this directory with your chosen project name.
Subfolders are created within this directory and these are norm, DifferentialExpression 
Norm - contains the normalised expression matrix and filtered expression matrices
Differential Expression - contains the lists of differentially expressed genes, and Venn diagrams 
cluster - contains results of unsupervised hierarchicial clustering analysis
GO - this folder contains the results of statistical analysis of Gene Ontology terms

<b>FUNCTIONS</b>
<b>1.Combat</b>
run_combat function performs applies the COMBAT algorithm to read counts data from different experiments (it is advised to use this option when conducting a meta-analysis)
usage: run_combat(target_file,norm_matrix,project)
<b>Arguments:</b>
target_file = Name of target file
norm_matrix = this is the name of the .txt file of read counts data (this can be either normalised or unnormalised)
project = name of project

<b>2.Filtering:</b>
Filtering function extracts the most variable genes passing a certain percentage filter from the normalised expression matrix
usage:filtering(filter,filterval,combat,project,data)
<b>Arguments:</b>
filter = Method of filtering to be used SD = standard deviation, 
filterval =  % to filter the normalised expression matrix - to retain the top 40% most variable genes filterval = 40.
combat = 0 if combat is not to be used and 1 if combat is to be used
project = project name
data = this is the name of the normalised/unnormalised .txt file of read counts data

<b>3. DifferentialExpression:</b>
Calculate the genes that are differentially expressed between groups e.g. two biological groups using the LIMMA R package
usage:run_limma(target_file,analysis,data,comp_file,replicates,limmamethod,adjust,pvalue,foldchange,platform,project,dataType)
<b>Arguments:</b>
target_file = name of target file
analysis = whether analysis is paired/unpaired
data = filtered normalised expression matrix
comp_file = text file containing the comparison between biological groups that are to be made - up to six different comparisons can be made
replicates = yes/no whether technical replicates are present
limmamethod = method to be used within limma e.g separate, global, nestedF
adjust = method to be used to adjust the false discovery rate (FDR) e.g BH = Benjamini-Hochberg
pvalue = value threshold genes below this value are listed as  differentially expressed
foldchange = log fold change value threshold genes with a fold change below this value are listed as differentially expressed
platform = Affymetrix platform that is to be analysed e.g for Affymetrix HGU133plus2 argument to be used is hgu133plus2
project = name of project
dataType = normalised or unnormalised

Calculate the genes that are differentially expressed between groups e.g. between two biological groups tumor vs normal (using the edgeR package)
usage:run_edgeR(data,target,analysis,comp_file,replicates,adjust,pvalue,project)
<b>Arguments:</b>
data = filtered normalised expression matrix
target = name of target file
analysis = whether analysis is paired/unpaired
comp_file = text file containing the comparison between biological groups that are to be made - up to six different comparisons can be made
replicates = yes/no whether technical replicates are present
adjust = method to be used to adjust the false discovery rate (FDR) e.g BH = Benjamini-Hochberg
pvalue = value threshold genes below this value are listed as  differentially expressed
project = name of project

<b>4. Annotation</b>
Annotate the list(s) of genes found to be differentially expressed from edgeR
usage:annotated_edgeR(comp_file,project)
<b>Arguments:</b>
comp_file = text file containing the comparison between biological groups that are to be made - up to six different comparisons can be made
project = name of project

<b>5. Gene Ontology analysis</b>
Perform Gene Ontology statistical analysis on the list of differentially expressed genes
usage:run_rna_seq_GO(project,comp_file,diffmethod,data)
<b>Arguments:</b>
project = name of project
comp_file = text file containing the comparison between biological groups that are to be made - up to six different comparisons can be made
diffmethod = limma (to use the LIMMA R package), edge to use the edgeR package
data = this is the name of the .txt file containing the normalised or unnormalised read counts data

<b>6. Clustering</b>
Performs non-hierarchical clustering on normalised expression matrix 
usage:cluster_all(target_file,k,project)\
<b>Arguments:</b>
target_file = name of target file
k = number of clusters to cluster the data
project - name of project
 

<b>7. Venn diagrams</b>
Generates venn diagrams from lists of differentially expressed genes
usage:generate_venn("decideTestsSummary.txt",project)
<b>Arguments:</b>
decideTestsSummary.txt = output from running limma to identify differentially expressed genes
project = name of project

<b>8. Survival</b>
Generation of Kaplan-Meier survival plots
<b>Arguments:</b>
usage:survival_os_plot(target_file,project)
target_file = name of target file
project - name of project

<b>9. Heatmap generation </b>
Generation of heatmaps
usage:heatmap_generate(target,project,comp_file)
<b>Arguments:</b>
target - full path to target.txt
project - name of project
comp_file - full path to comparisons.txt

=======
## Ominer-RNA-Seq

O-miner code for RNA-Seq post-processing


