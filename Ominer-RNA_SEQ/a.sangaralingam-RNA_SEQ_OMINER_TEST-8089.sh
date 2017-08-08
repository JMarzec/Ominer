#!/bin/sh
R --file=/data/BCI-BioInformatics/ominer/transcriptomics/RNA-Seq/run_RNA_SEQ.R --args dataType=unnormalised folder=/var/www/html/onlinetool/temp/a.sangaralingam-RNA_SEQ_OMINER_TEST-8089 project=a.sangaralingam-RNA_SEQ_OMINER_TEST-8089
target=/var/www/html/onlinetool/temp/a.sangaralingam-RNA_SEQ_OMINER_TEST-8089/target.txt comp=/var/www/html/onlinetool/temp/a.sangaralingam-RNA_SEQ_OMINER_TEST-8089/comp.txt normalisation=0 filter=sd filterval=40 adjust=BH pvalue=0.05 foldchange=2.0 diffmethod=edge
limmamethod= analysis=unpaired replicates=no combat=0 survival=0 

