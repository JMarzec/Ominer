R --file=/scratch/ominer/transcriptomics/Affymetrix/affy_expression.R  --args  platform=hgu133plus2 dataType=CEL folder=/scratch/ominer/transcriptomics/Affymetrix/CEL project=a.sangaralingam-AFFY_TEST target=/scratch/ominer/transcriptomics/Affymetrix/target.txt comp=/scratch/ominer/transcriptomics/Affymetrix/comp.txt normalisation=rma filter=sd filterval=40 adjust=BH pvalue=0.05 foldchange=0.5 limmamethod=separate analysis=unpaired replicates=no combat=0 survival=0 estimate=1 aqm=0 gostats=1 