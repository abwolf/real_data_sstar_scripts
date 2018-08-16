library(qvalue)

## SETUP ##
## grab command line arguments ## 
ARGV <- commandArgs(trailingOnly=TRUE)
tag <- ARGV[1]
samplistfile <- ARGV[2]
outputdir <- ARGV[3]
input.req.window.snps <- ARGV[4]



## set filtering thresholds ##
req.window.snps <- as.numeric(input.req.window.snps)
req.N <- 100  # minimum number of matching YRI haplotypes used for match p-val calculation
samps <- scan(samplistfile, what="")
sstar_pval_thresholds <- c(seq(0.05, 0.01, -0.01), 0.001)
match_qval_thresholds <- c(0.05, 0.01)

outfile <- paste0(outputdir, "/", tag, "_match_pvalue_thresholds.txt")
pvals.allsamps <- list()

ReadAndAnnotateSampleData <- function(infile) {
    dat <- read.table(infile, header=TRUE, na.strings=c("NA", "None"), as.is=TRUE)
    dat[which(dat=="None", arr.ind=TRUE)] <- NA

    ## window masking ##
    ## filter windows without S* value or that failed filters used for S* significance calculation
    ## filter windows based on number of S* SNPs
    
    dat$filter_sampspecific <- NA   ## new column
    dat$filter_sampspecific[dat$num_s_star_snps < req.window.snps] <- TRUE
   
    ## haplotype-specific masking ##
    ## based on N ##
    ## and on whether haplotype underlies S* score ##
    dat$filter_hap1 <- NA  ## new column
    dat$filter_hap2 <- NA  ## new column
    dat$hap_1_window_pval_table[dat$hap_1_window_match_N_table<req.N] <- NA
    dat$hap_2_window_pval_table[dat$hap_2_window_match_N_table<req.N] <- NA
    dat$hap_1_window_pval_table[dat$s_star_hap1 != TRUE] <- NA
    dat$hap_2_window_pval_table[dat$s_star_hap2 != TRUE] <- NA

    ## create and update new filter columns based on filters applied above
    dat$filter_sampspecific <- NA   ## new column
    dat$filter_sampspecific[dat$num_s_star_snps < req.window.snps] <- TRUE
   dat$filter_hap1 <- NA  ## new column
    dat$filter_hap2 <- NA  ## new column
    dat$filter_hap1[dat$hap_1_window_match_N_table<req.N] <- TRUE
    dat$filter_hap2[dat$hap_2_window_match_N_table<req.N] <- TRUE
    dat$filter_hap1[dat$s_star_hap_1 != TRUE] <- TRUE
    dat$filter_hap2[dat$s_star_hap_2 != TRUE] <- TRUE

    
    ## update existing filter column using new criteria (but this won't include haplotype-level filtering)
    dat$filter[dat$num_s_star_snps < req.window.snps] <- TRUE
    return(dat)

}

FindSstarSigWindows <- function(adat, threshold_sstarpval) {
    windows <- list()
    for(i in c(1,2)) {
        windows[[i]] <- subset(adat, (adat$sstarpval <= threshold_sstarpval & adat[,paste0("s_star_hap_", i)] & !adat$filter))
    }
    return(windows)
}

FindNotSigWindows <- function(adat, threshold_sstarpval) {
    windows <- list()
 for(i in c(1,2)) {
        windows[[i]] <- subset(adat, (adat$sstarpval > threshold_sstarpval & adat[,paste0("s_star_hap_", i)] & !adat$filter))
    }
    return(windows)
}


PrintWindows <- function(windowslist, prefix) {
    for(i in c(1,2)) {
            toprint <- windowslist[[i]]
            toprint <- toprint[,1:3]
            toprint$chrom <- paste0("chr", toprint$chrom)
            write.table(toprint, file=paste0(prefix, "_hap_", i, ".bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
        }
}
        

### PART 1 ###
### Use all samples to calculate match q-values ###
### and corresponding match p-value thresholds ### 

for(sampname in samps) {
    print(paste("Calculating q-values: processing sample", sampname))

    ## read in data #
    ## and annotate filter
    ## then write back out
    wcfile <- paste0(outputdir, "/SstarSigFiles/", tag, "_", sampname, "_null_glm.sstar_sig_out")
    dat  <- ReadAndAnnotateSampleData(wcfile)
    dat <- subset(dat, dat$filter==FALSE)

    
    ## save data from relevant haplotypes to all-sample object
    pvals.samp <- rbind(
        data.frame(matchpval=dat$hap_1_window_pval_table,sstarpval=dat$sstarpval), # grab hap1 values
        data.frame(matchpval=dat$hap_2_window_pval_table,sstarpval=dat$sstarpval)) # grab hap2 values

    pvals.samp <- subset(pvals.samp, !is.na(pvals.samp$matchpval))   # keep only windows with a match p-value
    pvals.allsamps[[length(pvals.allsamps)+1]] <- pvals.samp  
}

pvals.allsamps <- Reduce(rbind, pvals.allsamps) # convert from list of data frames to a single big data frame


## for each S* p-value threshold,
## and each Neanderthal match FDR value,
## output Neanderthal match p-value threshold

match_pval_thresholds <- matrix(nrow=length(sstar_pval_thresholds), ncol=length(match_qval_thresholds),
                                dimnames=list(paste0("p", sstar_pval_thresholds), paste0("p",match_qval_thresholds)))

for(xi in 1:length(sstar_pval_thresholds)) {
    x <- sstar_pval_thresholds[xi]
    matchpvals_forq<- subset(pvals.allsamps$matchpval, pvals.allsamps$sstarpval <= x)
    npvals <- length(matchpvals_forq)
    if(npvals>0) {
        print(paste0(npvals, " p-values being used to calculate q-values at S* sig value ", x ))
        qobj <- qvalue(matchpvals_forq, pi0.method="bootstrap")
    }else{
        print(paste0("no haplotypes passed thresholds to calculate q-values at S* sig value ", x ))
        qobj <- NA 
    }

    df <- data.frame(qvals=qobj$qvalues, pvals=qobj$pvalues)
    df <- df[order(df$qvals, df$pvals),]
    for(yi in 1:length(match_qval_thresholds)) {
        y <- match_qval_thresholds[yi]
        thresh <- df$pvals[tail(which(df$qvals <= y),1)]  # this will have value numeric(0) if no windows met the qvalue threshold...
        match_pval_thresholds[xi,yi] <- ifelse(length(thresh) > 0, thresh, -1) # in which case we will set the p-value threshold to -1

    }
}

## FINISH PART 1 ##
write.table(match_pval_thresholds, file=outfile, quote=FALSE, sep="\t")
#save.image(file="testing_qvalue_threshold_calculation.rdata")


### PART 2 ###
### Add column indicating Neanderthal match status at all match p-value thresholds ### 

for(sampname in samps) {
    print(paste("Determining match significance: processing sample", sampname))

    # read in data #
    # and filter 
    wcfile <- paste0(outputdir, "/SstarSigFiles/", tag, "_", sampname, "_null_glm.sstar_sig_out")
    dat  <- ReadAndAnnotateSampleData(wcfile)
    dat$chrom <- paste0("chr", dat$chrom) # to use bed file naming convention 

    for(i in c(1,2)) {
            colname <- paste0("hap_", i, "_match_sstarsig0.05_FDR0.05")
            threshold_matchpval <- match_pval_thresholds[paste0("p", 0.05), paste0("p", 0.05)]
            dat[,colname] <- (dat[,paste0("hap_", i,"_window_pval_table")] <=  threshold_matchpval & dat$sstarpval <= 0.05 & dat[,paste0("s_star_hap_", i)] & !dat$filter)

            colname <- paste0("hap_", i, "_match_sstarsig0.01_FDR0.05")
            threshold_matchpval <- match_pval_thresholds[paste0("p", 0.01), paste0("p", 0.05)]
            dat[,colname] <- (dat[,paste0("hap_", i,"_window_pval_table")] <=  threshold_matchpval & dat$sstarpval <= 0.05 & dat[,paste0("s_star_hap_", i)] & !dat$filter)
        }

    write.table(dat, file=paste0(wcfile, "_withmatchstatus"), row.names=FALSE, quote=FALSE, sep="\t")
}




####### FOR TEST
FindSigMatchingWindows <- function(adat, threshold_sstarpval, threshold_matchqval) {
    threshold_matchpval <- match_pval_thresholds[paste0("p", threshold_sstarpval), paste0("p", threshold_matchqval)]
    windows <- list()
    for(i in c(1,2)) {
        windows[[i]] <- subset(adat, (adat[,paste0("hap_", i,"_window_pval_table")] <= threshold_matchpval & adat$sstarpval <= threshold_sstarpval & adat[,paste0("s_star_hap_", i)] & !adat$filter))
    }
    return(windows)
}

if(0) {
    sigwindows <- FindSstarSigWindows(dat, 0.05)
    notsigwindows <- FindNotSigWindows(dat, 0.05)
    matchingwindows <- FindSigMatchingWindows(dat, 0.05, 0.05)

    PrintWindows(sigwindows, "tmp_sigwindows")
    PrintWindows(notsigwindows, "tmp_notsigwindows")
    PrintWindows(matchingwindows, "tmp_matchingwindows")

    print('PrintWindows')
}
