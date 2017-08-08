#!/bin/sh
 /usr/bin/R --slave --no-restore --no-save --no-readline --silent  --args email=r.j.cutts@qmul.ac.uk folder=/data/BCI-BioInformatics/ominer/RNA_seq/post_processing  project=test target=target.txt mf=1 count_table=normalised.txt combat=0 diffmethod=edge analysis=unpaired replicates=no filterval=40 filter=iqr foldchange=0.5 limmamethod=limma comp=RNA_seq_comp.txt adjust=BH pvalue=0.05 gostats=0 normalisation=rma  survival=0  < run_RNA_SEQ.R 

