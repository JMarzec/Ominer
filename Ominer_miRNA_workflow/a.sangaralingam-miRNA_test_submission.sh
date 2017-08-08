R --file=/scratch/ominer/transcriptomics/Affymetrix/miRNA_commands.R --args  email=a.sangaralingam@qmul.ac.uk platform=mirna2 dataType=CEL folder=/var/www/html/onlinetool/temp/a.sangaralingam-miRNA_test_submission
project=a.sangaralingam-miRNA_test_submission target=/var/www/html/onlinetool/temp/a.sangaralingam-miRNA_test_submission/target.txt comp=/var/www/html/onlinetool/temp/a.sangaralingam-miRNA_test_submission/comp.txt normalisation=RMA
filter=sd filterval=40 adjust=BH pvalue=0.05 foldchange=2.0 limmamethod=separate analysis=unpaired replicates=no combat=0 aqm=0 
