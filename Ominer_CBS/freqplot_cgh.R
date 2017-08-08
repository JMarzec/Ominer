############################################################################################################
#
#File freqplot_cgh.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: Script to generate frequency plots from binary_code_data and binary coded filtered data
#Input a text file of log2ratios - columns are sample names, probeids, chromosome and position are also present as column headers.
#Sample input file R_input.txt
#inputs to function are: (1) results = .txt file of binary coded data, (2) target = .txt of target file, (3) Cancer - name of project, (4) type = datatype i.e. CEL 





freq_plot_cgh <- function(target,results,type,Cancer) {
    results <- read.table(results, header = T, sep = "\t", as.is = T)
      targs <- read.table(target, sep = "\t", as.is = TRUE)
    targets = targs[1, ]
    pd <- read.table(targets, header = T, sep = "\t", as.is = T)
    chromosomes = c(1:22)
    levels = unique(pd$Group)
    par(mar = c(8, 4, 4, 4))
    for (chr in 1:length(chromosomes)) {
        for (lev in 1:length(levels)) {
            if (type == "CEL") {
            if(platform=="500" || platform=="100") {
                x = (pd[pd$Group == levels[lev], 3])
            }
            else {
             x = (pd[pd$Group == levels[lev], 2])
            }
            }
            else {
                x = (pd[pd$Group == levels[lev], 1])
            }
            results$Loss = rep(0, length(nrow(results)))
            results$Gain = rep(0, length(nrow(results)))
            data.plot = results[, c("ProbeID", "Chromosome", 
                "Position", x, "Gain", "Loss")]
            data.plot = subset(data.plot, Chromosome == chromosomes[chr])
            data.plot.subset = data.plot[order(as.numeric(data.plot$Chromosome), 
                as.numeric(data.plot$Position)), ]
            if (length(x) > 1) {
                data.plot.subset[, "Gain"] = rowSums(data.plot.subset[, 
                  x] == 1)
                data.plot.subset[, "Loss"] = rowSums(data.plot.subset[, 
                  x] == -1)
            }
            if (length(x) == 1) {
                data.plot.subset[, "Gain"] = data.plot.subset[, 
                  x]
                data.plot.subset[, "Loss"] = data.plot.subset[, 
                  x]
                data.plot.subset[, "Gain"] = replace(data.plot.subset[, 
                  "Gain"], data.plot.subset[, "Gain"] == -1, 
                  0)
                data.plot.subset[, "Loss"] = replace(data.plot.subset[, 
                  "Loss"], data.plot.subset[, "Loss"] == 1, 0)
                data.plot.subset[, "Loss"] = replace(data.plot.subset[, 
                  "Loss"], data.plot.subset[, "Loss"] == -1, 
                  1)
            }
            png(paste("ominer_results/",Cancer,"/FP/",levels[lev],"FrequencyPlots","chr",chr, ".png", sep = ""), width = 1024, height = 300)
            n = length(x)
            y1 = 100 * (data.plot.subset$Gain/n)
            y2 = paste("-", 100 * (data.plot.subset$Loss/n), 
                sep = "")
            plot(y1, type = "h", xaxt = "n", yaxt = "n", col = "green", 
                main = paste(levels[lev], n, sep = " "), ylim = range(-100, 
                  100), xlab = chromosomes[chr], ylab = "Fraction of Patient for Gain or Loss", 
                xaxs = "i", yaxs = "i")
            points(y2, type = "h", col = "red")
            x = c(-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 
                100)
            y = chromosomes[chr]
            axis(2, at = c(x), labels = x, tick = TRUE, las = 1, 
                col = "black", lty = "dotted", tck = 1)
            dev.off()
        }
    }
    dev.off()
    
}