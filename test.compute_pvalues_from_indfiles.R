## requires data.table
library(data.table)
library(dplyr)
library(stringr)
require(bit64)

## SETUP ##
ARGV <- commandArgs(trailingOnly=TRUE)
tag <- "YRI_1"
neand.callable.file <- "~/AKEYLab//Sstar_files/real_data_files.sv//windows.50kb.10kb.bed"
recomb_rates_file <- "~/AKEYLab//Sstar_files/real_data_files.sv//recombination_rates_per_50Kb_window_allchr.txt"
outputdir <- "~/DATALab/SimulatedDemographic/Sstar/chr1_variable_ref//Output_YRI_1"
nullmode <- "glm"
sstar_null_model<- "~/AKEYLab//Sstar_files/real_data_files.sv//gravel_asn_scale_60k.model.Rdata"
region_starts <- 9000001
region_chrs <- 1
refname <- "NA18522"
n_req_window_snps <- 3
spval <- 0.01
matchpval <- 0.05



try(if(!any(c(nullmode=="ecdf", nullmode=="glm"))) {stop("nullmode must be either 'ecdf' or 'glm'" ) } )
req.snp.frac <- 0.8  # hardwiring this for now, but can make a variable
# colstowrite <- c("chrom", "winstart", "winend", "callable_bases_neand", "n_snps", "n_ind_snps", "n_region_ind_snps", "ind_id",
#                  "pop", "s_star", "num_s_star_snps", "s_star_snps",
#                  "hap_1_window_pval_table", "hap_2_window_pval_table",
#                  "hap_1_window_match_pct_table", "hap_2_window_match_pct_table",
#                  "hap_1_window_match_N_table", "hap_2_window_match_N_table",
#                  "hap_1_window_match_len_table", "hap_2_window_match_len_table",
#                  "hap_1_window_match_mapped_table", "hap_2_window_match_mapped_table",
#                  "hap_1_window_match_mh_table", "hap_2_window_match_mh_table",
#                  "hap_1_s_start", "hap_1_s_end", "hap_2_s_start", "hap_2_s_end", "s_start", "s_end",
#                  "n_s_star_snps_hap1", "n_s_star_snps_hap2", "s_star_haps",
#                  "recomb", "filter", "filter_min_callable", "s_star_hap_1", "s_star_hap_2", "sstarpval")


colstowrite <- c("chrom", "winstart", "winend", "callable_bases_neand", "n_snps", "n_ind_snps", "n_region_ind_snps", "ind_id",
                  "pop", "s_star", "num_s_star_snps", "s_star_snps",
                  "hap_1_window_pval_table", "hap_2_window_pval_table",
                  "hap_1_s_start", "hap_1_s_end", "hap_2_s_start", "hap_2_s_end", "s_start", "s_end",
                  "n_s_star_snps_hap1", "n_s_star_snps_hap2", "s_star_haps",
                  "recomb", "filter", "filter_min_callable", "s_star_hap_1", "s_star_hap_2", "sstarpval")

outfile <- paste0(outputdir, "/SstarSigFiles/", tag, "_chr_", region_chrs, "_start_", region_starts,".",refname,".sstar_sig_out.gz")
outbed <- paste0(outputdir, "/SstarSigFiles/", tag, "_chr_", region_chrs, "_start_", region_starts,".",refname,".sstar_sig.bed.gz")
## load the object G2s if using glm model
## or object F_null  if using ecdf
print('Loading precomputed null model...')
load(sstar_null_model)


## load callable sequence
print('Loading callable sequence...')
neand.callable = fread(neand.callable.file,
                        col.names=c('chrom', 'winstart', 'winend', 'callable_bases_neand'))
neand.callable <- neand.callable %>% filter(chrom==region_chrs) %>%
                mutate(chrom=as.numeric(chrom),
                        winstart=as.numeric(winstart),
                        winend=as.numeric(winend),
                        callable_bases_neand=as.numeric(callable_bases_neand)) %>% as.data.table()
setkey(neand.callable, chrom, winstart, winend)
# Test neand.callable loaded correctly
head(neand.callable)



## FUNCTIONS ##
TestThreshold.fn <- function(anullmode, apval, adat) {
    totest <- which(adat$filter==FALSE)

    # figure out which windows met the pvalue threshold
    if(anullmode=="glm") {
        passed <- adat[totest, s_star] > predict(G2s, adat[totest, list(snps = n_region_ind_snps, lr = log(recomb), q=(1-apval))])
    }
    if(anullmode=="ecdf") {
        passed <- F_null(adat[totest, s_star]) > 1-apval
    }

    # if any windows passed the pvalue threshold, store the pvalue in 'sstarpval' column
    toreplace <- which(passed)
    if(length(toreplace) > 0) {
        adat$sstarpval[totest][toreplace] <- apval
    }
    return(adat)
}

write.filtered.bed.fn <- function(dt, outbed, spval, matchpval){
 print(' Writing filtered .bed file')
 dt_1 <- dt %>%
     filter(filter==FALSE) %>%
     filter(sstarpval<=spval) %>%
     filter(hap_1_window_pval_table<=matchpval | hap_2_window_pval_table<=matchpval) %>%
     filter(s_star_hap_1==TRUE) %>%
     filter(num_s_star_snps>=n_req_window_snps) %>%
     select(chrom, ind_id, winstart, winend) %>%
     #select(chrom, ind_id, hap_1_s_start, hap_1_s_end) %>%
     mutate(msp_ID=paste0(ind_id,':',1,'_',chrom)) %>%
     select(msp_ID, winstart, winend) %>%
     #select(msp_ID, hap_1_s_start, hap_1_s_end) %>%
     setnames(c('msp_ID','start','end')) %>%
     as.data.table()

 dt_2 <- dt %>%
    filter(filter==FALSE) %>%
    filter(sstarpval<=spval) %>%
    filter(hap_1_window_pval_table<=matchpval | hap_2_window_pval_table<=matchpval) %>%
    filter(s_star_hap_2==TRUE) %>%
    filter(num_s_star_snps>=n_req_window_snps) %>%
    select(chrom, ind_id, winstart, winend) %>%
    #select(chrom, ind_id, hap_1_s_start, hap_1_s_end) %>%
    mutate(msp_ID=paste0(ind_id,':',2,'_',chrom)) %>%
    select(msp_ID, winstart, winend) %>%
    #select(msp_ID, hap_1_s_start, hap_1_s_end) %>%
    setnames(c('msp_ID','start','end')) %>%
    as.data.table()

 dt.bed <- rbind(dt_1,dt_2)

 print(paste0('nrow dt_1: ',nrow(dt_1)))
 print(paste0('nrow dt_2: ',nrow(dt_2)))
 print(paste0('nrow dt.bed: ',nrow(dt.bed)))

 options(scipen=10)
 dt.bed %<>% mutate(start=str_trim(as.character(start),side = "both")) %>% mutate(end=str_trim(as.character(end), side = "both"))
 write.table(x = dt.bed,
           file = gzfile(outbed),
           quote = FALSE,
           sep = '\t',
           row.names = FALSE,
           col.names = TRUE)
 options(scipen=0)
}



## MAIN ##
wroteheader <- FALSE

  ## load the S* output as dat
print('Loading S* output...')
windowcalcfile <- paste0(outputdir, "/RegionFiles/", tag, "_chr_", region_chrs, "_start_", region_starts,".",refname,".windowcalc_out.gz")
dat.allchr = fread(input = paste0('gunzip -c ', windowcalcfile), select = colstowrite, header=TRUE)

## determine which haplotype(s) contributed to S* score
dat.allchr$s_star_hap_1 <- (dat.allchr$n_s_star_snps_hap1 / dat.allchr$num_s_star_snps) > req.snp.frac
dat.allchr$s_star_hap_2 <- (dat.allchr$n_s_star_snps_hap2 / dat.allchr$num_s_star_snps) > req.snp.frac

## load the recombination rates
print('Loading recombination rates...')
recomb_per_interval <- fread(recomb_rates_file, header=TRUE)
colnames(recomb_per_interval) <- c("chrom", "winstart", "recomb")


for (achr in unique(sort(dat.allchr$chrom))) { # real data has 22 chromosomes, sim data probably has only 1
## merge dt with recombination rates
print(paste0("chromosome ", achr))
print('Merging S* output and recombination rates...')
dat <- subset(dat.allchr, dat.allchr$chrom==achr)
dat <- merge(dat, recomb_per_interval, by=c("chrom", "winstart"), all.x=TRUE)
dat <- data.table(dat)
setkey(dat, chrom, winstart, winend)

##  merge dt with callable sequence
print('Merging S* output and callable sequence...')
dat = left_join(dat, neand.callable, by=c('chrom','winstart','winend')) %>% as.data.table()
dat[, min_callable := callable_bases_neand]


## filter windows
print('Filtering...')
dat[, filter := (recomb == 0 | min_callable < 25000 | s_star <= 0)]
dat[, filter_min_callable := min_callable < 25000]
dat[,"sstarpval" := 0.00]  # make new column of double type
dat[,"sstarpval" := NA]  # replace values with NA
dat[filter==FALSE, "sstarpval" := 1]  ## initiate with value 1 if window passed filters

print('Determining maximum significant p-value...')
pvalstotest <- c(0.05, 0.01, 0.001)
for(pp in pvalstotest) {
    dat <- TestThreshold.fn(nullmode, pp, dat)
}

print(paste0('nrow dat: ', nrow(dat)))

## output --
print('Writing to file...')
# if(!wroteheader) {
#
#     write.table(colstowrite, file=outfile, ncol=ncol(dat))
#     wroteheader <- TRUE
# }


    ## write expanded output file including significant windows only
# if(0) {
#     forfilter <- dat[,paste0('sig', qt), with=FALSE]
#     write.table(dat[as.logical(unlist(forfilter)),], file=outfile,
#                 append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
# }

write.filtered.bed.fn(dt=dat,
                    outbed="~/DATALab/SimulatedDemographic/Sstar/chr1_variable_ref/Output_YRI_1/SstarSigFiles/test.bed.gz",
                    spval=spval,
                    matchpval=matchpval)

    ## or write expanded output file including all windows
write.table(dat[,colstowrite, with=FALSE], file=gzfile(outfile),
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

}
