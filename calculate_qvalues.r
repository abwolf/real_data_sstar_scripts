library(qvalue)
library(GenomicRanges)


## SETUP ##
## grab command line arguments ## 
ARGV <- commandArgs(trailingOnly=TRUE)
tag <- ARGV[1]
samplistfile <- ARGV[2]
outputdir <- ARGV[3]
input.req.window.snps <- ARGV[4]
scriptdir <- ARGV[5]
source(paste0(scriptdir, "/calculate_qvalues_functions.r"))

## set filtering thresholds ##
## I have the last two hardwired for now, but 
## can change the required number of matching YRI haplotypes used for match p-val calculation
## and the  fraction of S* SNPs that are required to have variant allele haplotype
req.window.snps <- as.numeric(input.req.window.snps)
req.N <- 100
req.snp.frac <- 0.5

## more setup ##
samps <- scan(samplistfile, what="")
sstar_sig_thresholds <- c("0.95", "0.99", "0.999", "0.9995", "0.9999")
FDR_thresholds <- c("prefilter", "FDR30pc", "FDR20pc", "FDR10pc", "FDR5pc", "FDR1pc")

## these correspond to FDR_thresholds vector
## FDR equal to 1 is a hacky way of getting all windows that passed S* filter
FDR_thresholds_numeric <- c(1, 0.30, 0.20, 0.10, 0.05, 0.01) 

outfileprefix <- paste0(tag, "_match_pvalue_hists_min", req.window.snps, "SNPs")
sigcolor <- rgb(1,0,0,1) 
#nonsigcolor <- rgb(99,184,255, 150, maxColorValue=255) # transparent color
nonsigcolor <- rgb(99,184,255, 255, maxColorValue=255) # non-transparent color
qobjs <- list()
regions_matchfilter <- list()


## MAIN ##
pdf(paste0(outfileprefix, ".pdf"), width=10, height=7)
par(mfrow=c(2,3), oma=c(0,0,1,0), cex.main=0.8)

for(sampname in samps) {
    print(paste("processing sample", sampname))
    regions_sig <- GRangesList()
    
    ## read in data # 
    wcfile <- paste0(outputdir, "/SstarSigFiles/", tag, "_", sampname, "_null_glm", ".sstar_sig_out")
    dat <- read.table(wcfile, header=TRUE, na.strings=c("NA", "None"))
    dat[which(dat=="None", arr.ind=TRUE)] <- NA
    datprefilter <- dat # for debug

    ## filter windows #
    ## filter windows without S* value or that failed filters used for S* significance calculation
    dat <- subset(dat, dat$filter==FALSE) # this column comes from S* sig calc script

    ## filter windows based on number of S* SNPs
    dat <- subset(dat, dat$num_s_star_snps >= req.window.snps)
    
    ## haplotype-specific masking based on N
    dat$hap_1_window_pval_table[dat$hap_1_window_match_N_table<req.N] <- NA
    dat$hap_2_window_pval_table[dat$hap_2_window_match_N_table<req.N] <- NA

    ## define columns for fraction of S* snps on each haplotype
    ## since some sites may be homozyous for S* variant, these may sum to greater than 1
    ## then do haplotype-specific masking based on fraction of S* SNPs on haplotype
    dat$prop_s_star_snps_hap1 <- dat$n_s_star_snps_hap1 / dat$num_s_star_snps
    dat$prop_s_star_snps_hap2 <- dat$n_s_star_snps_hap2 / dat$num_s_star_snps

    dat$haps_sstar <- apply(dat, 1, function(x) {
                                if(x["prop_s_star_snps_hap1"] > req.snp.frac) {
                                    if(x["prop_s_star_snps_hap2"] > req.snp.frac) { return("both") }
                                    else{return("hap1") }
                                }
                                if(x["prop_s_star_snps_hap2"] > req.snp.frac) { return("hap2") }
                            })
    
    dat$hap_1_window_pval_table[dat$prop_s_star_snps_hap1 < req.snp.frac] <- NA
    dat$hap_2_window_pval_table[dat$prop_s_star_snps_hap2 < req.snp.frac] <- NA

    ## for each S* significance threshold, #
    ## plot match p-value histograms and calculate q-values
    qobjs[[sampname]] <- list()
    regions_matchfilter[[sampname]] <- list()
    plot.new()
    
    for(sigval in c("all", sstar_sig_thresholds)) {
        ## plot histogram
        if(sigval!= "all") {
        nplotted <- histfunc(dat, asigval=sigval, axlab="match p-value\nfrom both haplotypes",
                             abreaks=100, aylim=c(0,0.3))
        mtext(sampname, side=3, outer=TRUE, line=-1, font=2)
        Drawlegend(nplotted, sigval)
    }
        ## calculate qvalues
        if(sigval=="all") {
            rowstokeep <- 1:nrow(dat) 
        }else {
            rowstokeep <- which(dat[,paste0("sig", sigval)]==TRUE)  # keep only windows significant at current threshold
        }
        
        sigdat1 <- dat[rowstokeep, c("chrom", "winstart", "winend", "prop_s_star_snps_hap1",
                                     "num_s_star_snps", "hap_1_window_pval_table", "haps_sstar")] # info and pvals for hap 1
        sigdat2 <- dat[rowstokeep, c("chrom", "winstart", "winend", "prop_s_star_snps_hap2",
                                     "num_s_star_snps", "hap_2_window_pval_table", "haps_sstar")] # info and pvals for hap 2
        sigdat1$hapnum <- 1 # add label column to keep phase information 
        sigdat2$hapnum <- 2
        names(sigdat1) <- names(sigdat2) <- c("chrom", "winstart", "winend", "prop_s_star_snps_this_hap",
                                              "num_s_star_snps", "match_pval", "haps_sstar", "hapnum")
        sigdat <- rbind(sigdat1,sigdat2)
        rm(list=c("sigdat1", "sigdat2"))

        sigdat <- subset(sigdat, !is.na(sigdat$match_pval))
        pvals <- as.vector(unlist(sigdat$match_pval))
        npvals <- length(pvals)
        if(npvals >0) {
        print(paste0(length(pvals), " p-values being used to calculate q-values at S* sig value ", sigval))
        qobjs[[sampname]][[sigval]] <- qvalue(pvals, pi0.method="bootstrap")
        qobjs[[sampname]][[sigval]]$wininfo <- sigdat[, c("chrom", "winstart", "winend","prop_s_star_snps_this_hap", "num_s_star_snps", "hapnum")]
        sigdat$qvalues <- qobjs[[sampname]][[sigval]]$qvalues

        regions_sig[[sigval]] <-  Make_reduced_gr_and_preserve_hapinfo(sigdat)

        
    } else {
        print(paste0("no haplotypes passed thresholds to calculate q-values at S* sig value ", sigval))
        qobjs[[sampname]][[sigval]] <- NA
        qobjs[[sampname]][[sigval]]$wininfo <- NA
        sigdat$qvalues <- NA
    }
       
        ## determine which windows pass at several FDRs
        regions_matchfilter[[sampname]][[sigval]] <-
            lapply(1:length(FDR_thresholds), function(x) {
                       Make_granges(sampname, sigval, afdr=FDR_thresholds_numeric[x])
                   })
        names(regions_matchfilter[[sampname]][[sigval]]) <- FDR_thresholds
    } # end of loop for given sigval


    ## generate mutually exclusive regions
    ## with annotation of status
    ## and print to bed file
    for(sigval in c("0.95", "0.99")) {
   #for(sigval in c("0.99")) {
        regions_forprint.gr <- Make_granges_allstatus(dat, sstarsiglevel=sigval)
        toprint.df <- Convert_granges_to_dataframe(regions_forprint.gr)
        write.table(toprint.df, file=paste0(outputdir, "/", tag, "_", sampname, "_sstarsig", sigval, "_matchFDR5pc.bed"),
                    row.names=FALSE, quote=FALSE, sep="\t")
    }
}


## SUMMARIZE USING QVALUES ## 
## make data frame of total Mb recovered at different FDRs for each individual
Totmb <- lapply(samps, function(sampname) {
                    sapply(c("all", sstar_sig_thresholds), function(thresh1) {
                               sapply(FDR_thresholds, function(thresh2) {
                                          sum(as.numeric(width(regions_matchfilter[[sampname]][[thresh1]][[thresh2]])))/1e6
                                      })})})
names(Totmb) <- samps

Totmb_sig0.95_matchFDR5pc <- sapply(Totmb, function(x) { x["FDR5pc", "0.95"] })  
write.table(data.frame(samp=names(Totmb_sig0.95_matchFDR5pc), tot_sig_mb=Totmb_sig0.95_matchFDR5pc),
                 file=paste0(tag, "_totmb_sig0.95_matchFDR5pc.txt"), row.names=FALSE, quote=FALSE, sep="\t")

Totmb_matchFDR5pc <- list()
for(yy in c("0.95", "0.99") ) {
    Totmb_matchFDR5pc[[yy]] <- sapply(Totmb, function(x) { x["FDR5pc", yy] })  
    write.table(data.frame(samp=names(Totmb_matchFDR5pc[[yy]]), tot_sig_mb=Totmb_matchFDR5pc[[yy]]),
                file=paste0(tag, "_totmb_sig", yy, "_matchFDR5pc.txt"), row.names=FALSE, quote=FALSE, sep="\t")
}

## FINISH ## 
save.image(paste0(outfileprefix, ".rdata"))



####################################################################
####################################################################

#if(0) {
### plot Mb passed filter for combinations of S* and FDR
#plotfunc <- function(afdr) {
#    acolor <- ifelse(afdr=="prefilter", "SlateBlue4", "DeepPink2")
#    aylab <- ifelse(afdr=="prefilter", "MB passed filters", "Mb passed filters")
#    atitle <- ifelse(afdr=="prefilter", "Before Neanderthal match filter", afdr)
#    ymax <- max(sapply(samps, function(x) {max(Totmb[[x]][afdr,]) }))
#
#    plot(1,1, xlim=c(1,5), ylim=c(0,ymax), type="n", xaxt="n", xlab="S* significance threshold", ylab=aylab, main=atitle);
#    for(x in samps) { points(Totmb[[x]][afdr,], type="o", col=acolor) };
#    axis(c("all", sstar_sig_thresholds), at=seq(1,6,1), side=1)
#
#}
#
#par(mfrow=c(3,3), oma=c(1,1,5,1)) # leave enough oma for filter info
#for (f in FDR_thresholds) { plotfunc(f) }
#filter.info <- paste(
#    paste0("Minimum ", req.window.snps, " S* SNPs in window,"),
#    paste0("Minimum ", req.snp.frac, " of S* SNPs on haplotype,"),
#    paste0("Minimum ", req.N, " matching YRI haplotypes"),
#    paste0("plus filters from S* significance calculation"),
#    sep="\n"
#    )
#mtext(filter.info, side=3, outer=TRUE, line=-1, font=1)
#dev.off()
#
### plot results across FDRs for a given S* threshold
### in terms of fraction of significant windows
#
#frac.sig <- lapply(1:length(qobjs), function(asampnum) {
#                       amat <- matrix(nrow=length(FDR_thresholds_numeric), ncol=length(sstar_sig_thresholds),
#                                      dimnames=list(fdr=FDR_thresholds, sstar_sig=sstar_sig_thresholds))
#                       for(y in 1:length(sstar_sig_thresholds)) {
#                           qvals.thisind <- qobjs[[asampnum]][[sstar_sig_thresholds[y]]][["qvalues"]]
#                           for(x in 1:length(FDR_thresholds)) {
#                                   amat[x,y] <- length(which(qvals.thisind <= FDR_thresholds_numeric[x]))
#                           }
#                       }
#                       return(amat)
#                   })
#
#
#plot(1,xlim=c(1,length(FDR_thresholds)), ylim=c(0,1), type="n", xaxt="n",
#     ylab="fraction of windows passed FDR threshold", xlab="FDR threshold")
#axis(side=1, at=1:length(FDR_thresholds), FDR_thresholds)
#
#asstar.sig <- "0.99"
#for(asamp in 1:length(frac.sig)) {
#        n.sig <- frac.sig[[asamp]][,asstar.sig]
#        frac.sig.thisind <- n.sig / n.sig[1]
#        lines(1:length(frac.sig.thisind), frac.sig.thisind, type="b")
#    }
#}
#
#
#
### WORKING ## 
#if(0) {
#
## plot histogram for just one sample for one sig value for S* significant windows only,
#    # and highlight FDR
#    sampname <- samps[14]
#    wcfile <- paste0("sstar5_mpval_", sampname, ".sstar_sig_out")
#    dat <- read.table(wcfile, header=TRUE, na.strings=c("NA", "None"))
#    dat[which(dat=="None", arr.ind=TRUE)] <- NA
#    datprefilter <- dat # for debug
#
#    # filter windows #
#    # filter windows without S* value or that failed filters used for S* significance calculation
#    dat <- subset(dat, dat$filter==FALSE) # this column comes from S* sig calc script
#    
#    # haplotype-specific masking based on N
#    dat$hap_1_window_pval_table[dat$hap_1_window_match_N_table<100] <- NA
#    dat$hap_2_window_pval_table[dat$hap_2_window_match_N_table<100] <- NA
#    sigval <- sstar_sig_thresholds[5]
#    nplotted <- histfunc(dat, asigval=sigval, axlab="Neanderthal match p-value", abreaks=100, sigonly=TRUE, aylim=c(0,0.3))
#    mtext(sampname, side=3, outer=TRUE, line=-1.5, font=2)
#    Drawlegend(nplotted, sigval, sigonly=TRUE)
#    mtext("Distribution of Neanderthal match p-values from S*-significant windows", side=3, line=1.5, font=2)
#
#
#    #### looking at haplotypes from window separately ####
#    pdf("plot_haps_separately.pdf")
#    par(mfrow=c(2,1))
#    plot(hmin, col="grey", ylim=c(0,0.25), xlab="p-value", main="Distribution of lower match p-values");
#    plot(hmax, col="grey", ylim=c(0,0.25), xlab="p-value", main="Distribution of higher match p-values")
#    mtext("from 1406 significant windows that had match p-values for both haplotypes", side=1, outer=TRUE, line=-1)
#    dev.off()
#
#
#    pdf("simulated_null.pdf");
#    par(mfrow=c(2,1));
#    plot(smin.h, col="grey");
#    plot(smax.h, col="grey");
#    mtext("simulated 1406 pairs of p-values from uniform distribution", side=1, line=-1, outer=TRUE);
#    dev.off()
#
#    # test GenomicRanges width function
#
#    df1 <- data.frame(seqnames=rep(1,10), start=c(1,11, 21,31,41,51,61,71,81,91), end=c(15,25,35,45,55,65,75,85,95,105))
#    gr1 <- GRanges(seqnames=df1$seqnames, ranges=IRanges(start=df1$start, end=df1$end))
#    gr1.red <- reduce(gr1)
#                         
#}
