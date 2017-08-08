##########################################################################################
#
#File name: KM_surv_plots.R
#
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#
#
#
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
##########################################################################################
##########################################################################################
#Description: function to perform survival analysis generates KM plots for 5,10 and 15 years
#Arguments are:
#target_file = path to the target file (please note: target file for a survival analysis MUST contain four columns headed sample,survtime (in months),survstat (this is survival status 0 for not alive and 1 for alive), and also riskgroups this is the same as the target column). We take the riskgroups automatically from the target file (and assign them as either riskgroup 1 or riskgroup 2).
#project is the given project name
################################################################################################KM plots - overall survival from a target file generates .txt file and KM plots for 5,10 and 15 years
######
survival_os_plot <- function(target_file,project) {
	
library("FC14.plotting.lib");
library("gplots");
my.round <- function(number) {
        rounded.value <- round(number, digits = 3)
        return(rounded.value)
}
cox.file <- paste( project, "_univariate_modeling.txt", sep = "" );
x.label <- "Survival time";
y.label <- "Survival probability";
res <- 400;
ann.data <- as.matrix( read.table( file = target_file, row.names = 1, header = T, sep = "\t" ) );
ann.data <- ann.data[ complete.cases(ann.data[ , "Surv_Period" ]) , ];


surv.time <- as.numeric( ann.data[ , "Surv_Period"] );
surv.stat <- as.numeric( ann.data[ , "Surv_Status"] );

samples.rg1 <- vector();
samples.rg2 <- vector();

rg <- rep( 0, length( 1:nrow(ann.data) ) );
names(rg) <- rownames(ann.data);

# identify riskgroup names
#unique.rgs <- unique( ann.data[ , "Riskgroups"] );
unique.rgs <- unique(ann.data[ , "Target"])
samples.rg1 <- names( which( ann.data[ , "Target"] == unique.rgs[1] ) );
samples.rg2 <- names( which( ann.data[ , "Target"] == unique.rgs[2] ) );    

cat("\nSamples in riskgroup1 = ", length(samples.rg1),"\tSamples in riskgroup2 ", length(samples.rg2), sep = "");

rg[samples.rg1] <- 1;

# apply univariate model to the data
cox.fit <- summary( coxph( Surv(surv.time, surv.stat) ~  rg  ) ); 

# vector of results to be written to outfile
cox.results <- c(
 "HR" = my.round(cox.fit$conf.int[1,1]), "CI95L" = my.round(cox.fit$conf.int[1,3]), "CI95U" = my.round(cox.fit$conf.int[1,4]), "WaldP" = cox.fit$coef[1,5] 
);
write.table( t(cox.results), file = paste("ominer_results/",project,"/","survival","/coxfile",sep=""), row.names = "", col.names = NA, sep = "\t" );

# KM PLOTS
# generate KM plots for 5, 10 and 15 year truncation points
for( trunc.surv in c(5, 10, 15) ) {

		filename = paste( "ominer_results/",project,"/","survival/", "_survival_analysis_",  trunc.surv, ".tiff", sep = "" );

		# check
		cat("\ncreating KM file: ", filename);

		make.KM.curve(
				riskgroup = rg,
				survtime = surv.time,
				survstat = surv.stat,
				file.name = filename, 
				truncate.survival = trunc.surv,
				cex.lab = 1.1,
				cex.axis = 1.0,
				main.title = paste(project, " (n=" , nrow(ann.data), ")", sep = "" ),
				xaxis.label = x.label,
				yaxis.label = y.label,
				resolution = res
		);
}

####Now convert images so they can be correctly displayed on the website
change_image_dir <- paste( "ominer_results/",project,"/","survival/",sep="")
setwd(change_image_dir)
system("convert _survival_analysis_5.tiff -resize 25% _survival_analysis_5.png")
system("convert _survival_analysis_10.tiff -resize 25% _survival_analysis_10.png")
system("convert _survival_analysis_15.tiff -resize 25% _survival_analysis_15.png")
setwd("../../../")
}
	
	
	
	
