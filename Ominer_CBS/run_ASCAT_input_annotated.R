############################################################################################################
#
#File name run_ASCAT_input.R
#Authors: Ajanthah Sangaralingam (a.sangaralingam@qmul.ac.uk)
#Barts Cancer Institute
#Queen Mary, University of London
#Charterhouse Square, London EC1M 6BQ
############################################################################################################
##################################################################################################################
#Description: R scripts for calculating log2ratios from raw CEL files using aroma.affymetrix works for all Affymetrix SNP platforms ##########
#writes out LRR and BAF values for each sample to separate .txt files #########################################

#####Requires the following R packages to be installed (1) aroma.affymetrix, (2) sfit, (3) aroma.cn and (4) calmate
#####Input arguments to this function are (1) locations - this is a .txt file with the full paths to the directory containing the rawdata for the normal sa
#######mples and the tumor samples. Each directory is on a separate line, (2) targets - .txt file containing the names of two target files .txt one for normals and one 
#######for tumor samples, (3) platform - this is the name of the Affymetrix platform type (4) dataSet - this relates to the .cdf file and depends o which platform is being analysed
########this is automatically calculated in CNV_shell.R, (5) anakysisType = paired/unpaired analysis.



run_ASCAT_input <- function (locations, targets, platform, dataSet,analysisType) {
    library("aroma.affymetrix")
    library("sfit")
    library("aroma.cn")
    library("calmate")
    verbose <- Arguments$getVerbose(-10, timestamp = TRUE)
    setOption(aromaSettings, "memory/ram", 200)
    locs <- read.table(locations, sep = "\t", as.is = T)
    Cancer = locs[1, ]
    targs <- read.table(targets, sep = "\t", as.is = TRUE)
    targets = targs[1, ]
    pd <- read.table(targets, sep = "\t", as.is = TRUE, strip.white = TRUE, 
        header = T)
         if (platform == "100") {
        cdf_arrays = c("Mapping50K_Hind240", "Mapping50K_Xba240")
        plm = "RmaCnPlm"
        if (analysisType == "paired") {
            chip_pairs = targs[2, ]
        }
    }
if (platform == "50xba") {
        cdf_arrays = "Mapping50K_Xba240"
        plm = "RmaCnPlm"
        if (analysisType == "paired") {
            chip_pairs = targs[2, ]
        }
    }
    if (platform == "50hind") {
        cdf_arrays = "Mapping50K_Hind240"
        plm = "RmaCnPlm"
        if (analysisType == "paired") {
            chip_pairs = targs[2, ]
        }
    }
    if (platform == "500") {
        cdf_arrays = c("Mapping250K_Sty", "Mapping250K_Nsp")
        plm = "RmaCnPlm"
        if (analysisType == "paired") {
            chip_pairs = targs[2, ]
        }
    }
    if (platform == "250sty") {
        cdf_arrays = "Mapping250K_Sty"
        plm = "RmaCnPlm"
        if (analysisType == "paired") {
            chip_pairs = targs[2, ]
        }
    }
    if (platform == "250nsp") {
        cdf_arrays = "Mapping250K_Nsp"
        plm = "RmaCnPlm"
        if (analysisType == "paired") {
            chip_pairs = targs[2, ]
        }
    }
    if (platform == "six") {
        cdf_arrays = "GenomeWideSNP_6,Full"
        plm = "AvgCnPlm"
        if (analysisType == "paired") {
            chip_pairs = targs[2, ]
        }
    }
    if (platform == "five") {
        cdf_arrays = "GenomeWideSNP_5,Full,r2"
        plm = "AvgCnPlm"
        if (analysisType == "paired") {
            chip_pairs = targs[2, ]
        }
    }
    if (platform == "131") {
        cdf_arrays = "Mapping10K_Xba131"
        plm = "RmaCnPlm"
        if (analysisType == "paired") {
            chip_pairs = targs[2, ]
        }
    }
    if (platform == "142") {
        cdf_arrays = "Mapping10K_Xba142"
        plm = "RmaCnPlm"
        if (analysisType == "paired") {
            chip_pairs = targs[2, ]
        }
    }
for (j in 1:length(cdf_arrays)) {
        chipType = cdf_arrays[j]
    cdfs <- lapply(chipType, FUN = function(chipType) {
        AffymetrixCdfFile$byChipType(chipType)
    })
    gis <- lapply(cdfs, getGenomeInformation)
    dataf <- lapply(gis, FUN = function(chipType) {
        readDataFrame(chipType)
    })
    sis <- lapply(cdfs, getUnitNames)
    data = do.call("rbind", dataf)
    data$unitName = unlist(sis)
    csR <- AffymetrixCelSet$byName(Cancer, chipType = chipType)
    print(csR)
    dsList <- doASCRMAv2(csR, plm = "RmaCnPlm", verbose = verbose)
    print(dsList)
    cmt <- CalMaTeCalibration(dsList)
    print(cmt)
    dsCList <- process(cmt, verbose = verbose)
    print(dsCList)
    extractSignals <- function(dsList, sampleName, reference = c("none", 
        "median"), refIdxs = NULL, ..., verbose = FALSE) {
        reference <- match.arg(reference)
        idx <- indexOf(dsList$total, sampleName)
        dfT <- getFile(dsList$total, idx)
        dfB <- getFile(dsList$fracB, idx)
        tcn <- extractRawCopyNumbers(dfT, logBase = NULL, ..., 
            verbose = verbose)
        baf <- extractRawAlleleBFractions(dfB, ..., verbose = verbose)
        if (reference == "median") {
            if (!is.null(refIdxs)) {
                dsR <- extract(dsList$total, refIdxs)
            }
            else {
                dsR <- dsList$total
            }
            dfTR <- getAverageFile(dsR, verbose = verbose)
            tcn <- divideBy(tcn, tcnR)
            setSignals(tcn, 2 * getSignals(tcn))
        }
        list(tcn = tcn, baf = baf, data = data)
    }
    rootPath <- "totalAndFracBData"
    dsfracb <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType = "*", 
        paths = rootPath)
    dstotalcn <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType = "*", 
        paths = rootPath)
            
    for (j in 1:length(pd$Name)) {
        sampleName <- pd$Name[j]
        dsCNList <- list(total = dstotalcn, fracB = dsfracb)
        chromo_vec <- vector()
        position_vec <- vector()
        log2_vec <- vector()
        betat_vec <- vector()
    
        
        for (chromo in 1:22) {
            dataTC_T <- extractSignals(dsCNList, sampleName = sampleName, 
                chromosome = chromo, verbose = verbose)
            log2_val <- log2(dataTC_T$tcn$cn/2)
            chromosome = dataTC_T$tcn$chromosome
            x = dataTC_T$tcn$x
            betaT = dataTC_T$baf$y
            chromo_vec <- c(chromo_vec, chromosome)
            position_vec <- c(position_vec, x)
            log2_vec <- c(log2_vec, log2_val)
            betat_vec <- c(betat_vec, betaT)
          
        }
        if (platform == "100" || platform == "500") {
        #if (platform == "500") {
            data1 <- data.frame(Chr = chromo_vec, Position = position_vec, 
            LRR = log2_vec)
        data2 <- data.frame(Chr = chromo_vec, Position = position_vec, 
            BAF = betat_vec)
        write.table(data1, paste("LRR_", sampleName, "_", chipType,".txt", sep = ""), 
            sep = "\t",quote = FALSE, row.names = TRUE)
        write.table(data2, paste("BAF_", sampleName, "_", chipType,".txt", sep = ""), 
            sep = "\t",quote = FALSE, row.names = TRUE)
        #}
        
        
        }
      if (platform == "50xba" || platform == "50hind" || platform == "250sty" || platform == "250nsp" || platform == "six"|| platform == "five") {
        data1 <- data.frame(Chr = chromo_vec, Position = position_vec, 
            LRR = log2_vec)
        data2 <- data.frame(Chr = chromo_vec, Position = position_vec, 
            BAF = betat_vec)
        write.table(data1, paste("LRR_", sampleName, ".txt", sep = ""), 
            sep = "\t", quote = FALSE, row.names = TRUE)
        write.table(data2, paste("BAF_", sampleName, ".txt", sep = ""), 
            sep = "\t", quote = FALSE, row.names = TRUE)
            }
    
    }
    
    
    
    
    }
    
}
