## FUNCTIONS ##
histfunc <- function(adf, asigval="0.99", acolor.sig=sigcolor, acolor.nonsig=nonsigcolor,
                     abreaks=100, aylim=c(0,0.2), axlab="p-value", sigonly=FALSE) {

    Get.freq.hist <- function(avec, N, b=abreaks) {
        h <- hist(avec, breaks=b, plot=FALSE)
        h$counts <- h$counts/N
        return(h)
    }

    # get just the two pval columns
    pvaldat <- adf[,c("hap_1_window_pval_table", "hap_2_window_pval_table")]  
    
    toplot <- list()
    toplot[["sig"]] <- subset(pvaldat, adf[,paste0("sig", asigval)]==TRUE)
    toplot[["nonsig"]] <- subset(pvaldat, adf[,paste0("sig", asigval)]==FALSE)
    toplot[["all"]] <- rbind(toplot[["sig"]], toplot[["nonsig"]])
                             
    plotvec <- lapply(toplot, function(x) { as.numeric(as.vector(unlist(x))) } ) 
    nplotted <- sapply(plotvec, function(x) { length(which(!is.na(x)))} )
    N <- nplotted[["all"]]

    if(sigonly) {
        hobj <- Get.freq.hist(plotvec[["sig"]], nplotted[["sig"]])
        plot(hobj, col=acolor.sig, add=FALSE, xlab=axlab, ylim=aylim,  main="")
    } else {
    # plot sig (actually plotting all, but will overlay with nonsig so sig will appear 'stacked' 
    hobj <- Get.freq.hist(plotvec[["all"]], N)
    plot(hobj, col=acolor.sig, add=FALSE, xlab=axlab, ylim=aylim,  main="")

    # add nonsig
    if(nplotted[["nonsig"]]>0) { # to account for situation where there are no windows in this category
        hobj <- Get.freq.hist(plotvec[["nonsig"]], N)
        plot(hobj, col=acolor.nonsig, add=TRUE)
    }
}
return(nplotted)
}


pN <- function(x) {
   prettyNum(x, big.mark=",")
}


Drawlegend <- function(anplotted, asigval, asigcolor=sigcolor, anonsigcolor=nonsigcolor, sigonly=FALSE) {

    nonsigtext <- paste0(pN(anplotted[["nonsig"]]), " windows - ", "S* not significant")
    sigtext <- paste0(pN(anplotted[["sig"]]), " windows - ", "S* significant (p<=" , round(1- as.numeric(asigval), 4),  ")" )
    if(sigonly) {
        textvec <- c(sigtext)
        colvec <- c(asigcolor)
    } else {
        textvec <- c(sigtext, nonsigtext)
        colvec <- c(asigcolor, anonsigcolor)
    }
   
    legend(textvec, col=colvec, 
           x="topright", bty="n", pch=15)
}

Make_granges <- function(asampname, asigval, afdr=NA) {
    ## assume that input coordinates are given in bed convention,
    ## that is the start coordinate is included in the region
    ## and the end coordinate is not included in the region
    ## The granges convention is both are included in the region,
    ## so adjust end coordinate accordingly

    if(nrow(qobjs[[asampname]][[asigval]]$wininfo)==0) { # no windows pass S* threshold
         print(paste0("no significant windows at S* threshold ", asigval))
        return(GRanges())
    }
    
    if(is.na(afdr)) {
        df <- qobjs[[asampname]][[asigval]]$wininfo # use all windows
    }
    else{ # use only windows that pass some FDR threshold
        df <- subset(qobjs[[asampname]][[asigval]]$wininfo, qobjs[[asampname]][[asigval]]$qvalues <= afdr) 
    }
    if(nrow(df) == 0) { # no windows significant at given FDR
        print(paste0("no matching windows at sigval ", asigval, " and FDR ", afdr))
        return(GRanges())
    }
    
    #gr <- GRanges(seqnames=df$chrom, ranges=IRanges(start=df$winstart, end=df$winend - 1))
    #gr.red <- reduce(gr) # the GenomicRanges reduce function
    #return(gr.red)
    
    # divide df by which haplotype the q-value came from
    dfbyhap <- lapply(1:2, function(x) { subset(df, df$hapnum==x) })
    grbyhap <- GRangesList(lapply(1:2, function(y) {
                          x <- dfbyhap[[y]]
                          tmp <- GRanges(seqnames=x$chrom, ranges=IRanges(start=x$winstart, end=x$winend - 1))
                          tmp <- reduce(tmp)
                          if(length(tmp) >0) { tmp$id <- y }
                          return(tmp)
                                     }))
    ## assign haplotype info to each region
    ## indicating whether hap1, hap2, or both met the match FDR threshold
    onlyhap1 <- setdiff(grbyhap[[1]], grbyhap[[2]])
    onlyhap2 <- setdiff(grbyhap[[2]], grbyhap[[1]])
    bothhaps <- intersect(grbyhap[[1]], grbyhap[[2]])

    if(length(onlyhap1) > 0)  { onlyhap1$id <- "hap1" }
    if(length(onlyhap2) > 0)  { onlyhap2$id <- "hap2" }
    if(length(bothhaps) > 0)  { bothhaps$id <- "both"}

    combined <- c(onlyhap1, onlyhap2, bothhaps)

    return(combined)
}     

 Make_granges_allstatus <- function(datobj, sstarsiglevel="0.95") {
     ## pull out windows for which at least one haplotype was considered
     ## and make reduced GRanges object
     dat.fil <- subset(datobj, (!is.na(datobj$hap_1_window_pval_table) | !is.na(datobj$hap_2_window_pval_table)))
     regions_considered <- Make_reduced_gr_and_preserve_hapinfo(dat.fil)

    ## make mutually exclusive sets of ranges: 
    ## regions passed both Sstar and archaic match thresholds   -- i.e. "sig0.95_matchFDR5pc"
    ## regions passed Sstar threshold but failed archaic match threshold    -- i.e. "sig0.95"
    ## regions failed Sstar threshold  -- "not_sstar_sig"
     ## I have to do it using tmps because GRangesList objects must all have the same structure BEFORE they're added to list :<
     ## so I have to add the id column 
     regions_forprint.gr <- GRangesList()
     regions_forprint.gr[[paste0("sstar_sig",sstarsiglevel,"_matchFDR5pc")]] <- regions_matchfilter[[sampname]][[sstarsiglevel]][["FDR5pc"]]

     ## and have to do setdiff within a hap category so I can reassign id, which is lost during setdiff, in bulk
     #tmp.sstarsig <- GRangesList()
     #tmp.notsig <- GRangesList()

      Setdiff_and_preserve_id <- function(ahapcat, gr1, gr2) {
         tmp.thishap <- subset(gr1, gr1$id==ahapcat)
         tmp.setdiff <- setdiff(tmp.thishap, gr2)
          if(length(tmp.setdiff) >0) { tmp.setdiff$id <- ahapcat } 
         return(tmp.setdiff)
     }
     
     tmp.sstarsig <- GRangesList(lapply(c("hap1", "hap2", "both"), Setdiff_and_preserve_id,
                                        regions_sig[[sstarsiglevel]],
                                        regions_forprint.gr[[paste0("sstar_sig", sstarsiglevel, "_matchFDR5pc")]]))
     tmp.notsig <- GRangesList(lapply(c("hap1", "hap2", "both"), Setdiff_and_preserve_id,
                                      regions_considered,
                                      regions_sig[[sstarsiglevel]]))
     names(tmp.sstarsig) <- names(tmp.notsig) <- c("hap1", "hap2", "both")

     # now there may be some overlapping windows between hap subsets
     # any overlaps between the "hap1" and "hap2" sets will be marked as NA, and 
     # anything in "both" will be removed from
     # the "hap1",  "hap2", and "NA" sets
     

     Setdiff_between_hapcats <- function(grlist) {
         tosetNA <- intersect(grlist[["hap1"]], grlist[["hap2"]]) # overlapping intervals annotated to opposite chromosomes
         
         newhap1 <- setdiff(grlist[["hap1"]], grlist[["both"]])
         newhap2 <- setdiff(grlist[["hap2"]], grlist[["both"]])
         tosetNA <- setdiff(tosetNA, grlist[["both"]])
         
         newhap1 <- setdiff(newhap1, tosetNA)
         newhap2 <- setdiff(newhap2, tosetNA)

         if(length(newhap1) >0) { newhap1$id <- "hap1" }
         if(length(newhap2) >0) { newhap2$id <- "hap2" }
         if(length(tosetNA  ) >0) { tosetNA$id <- "NA" }

         return(GRangesList(newhap1, newhap2, tosetNA, grlist[["both"]]))
     }

     tmp.sstarsig.nooverlap <- Setdiff_between_hapcats(tmp.sstarsig)
     tmp.notsig.nooverlap <- Setdiff_between_hapcats(tmp.notsig)
     
     regions_forprint.gr[[paste0("sstar_sig", sstarsiglevel)]] <- unlist(tmp.sstarsig.nooverlap)
     regions_forprint.gr[[paste0("not_sstar_sig", sstarsiglevel)]] <- unlist(tmp.notsig.nooverlap)
     
     return(regions_forprint.gr)
 }


Convert_granges_to_dataframe <- function(grobj) {
    regions_forprint.df <- lapply(names(grobj), function(astatus) { # make a list of data frames
                                      obj <- grobj[[astatus]]
                                      if(length(obj) ==0) { # return dummy  data frame
                                          return(data.frame(chrom=NA, chromStart=NA, chromEnd=NA, status=NA, hap=NA, width=NA))
                                      }  else {
                                          df <- data.frame(chrom=seqnames(obj), chromStart=start(obj), chromEnd=end(obj)+1,
                                                           status=astatus, hap=obj$id, width=width(obj))
                                          return(df)
                                      }})
    toprint <- Reduce(rbind, regions_forprint.df) # unlist the data frames
    toprint <- toprint[!is.na(toprint$chrom),] # remove any NA lines that were inserted if a status category was empty
    toprint <- toprint[order(as.numeric(as.vector(toprint$chrom)), as.numeric(as.vector(toprint$chromStart))), ] # order by start position
    return(toprint)
}


Make_reduced_gr_and_preserve_hapinfo <- function(adf, hapcolumn="haps_sstar") {
    ## divide df by which haplotype(s) meet threshold for proportion of S* SNPs
    ## then reduce each and add id column 
    dfbyhap <- lapply(c("hap1", "hap2", "both"), function(x) { subset(adf, adf[,hapcolumn]==x) })
    names(dfbyhap) <- c("hap1", "hap2", "both")
    grbyhap <- GRangesList(lapply(c("hap1", "hap2", "both"), function(y) {
                          x <- dfbyhap[[y]]
                          tmp <- GRanges(seqnames=x$chrom, ranges=IRanges(start=x$winstart, end=x$winend - 1))
                          tmp <- reduce(tmp)
                          if(length(tmp) >0) { tmp$id <- y }
                          return(tmp)
                                     }))
    
    ## use unlist to combine GRangesList elements into one GRanges object
    return(unlist(grbyhap))
}

Make_reduced_gr_and_preserve_hapinfo_noboth <- function(adf, hapcolumn="haps_sstar") {
    ## divide df by which haplotype(s) meet threshold for proportion of S* SNPs
    ## then reduce each and add id column 
    dfbyhap <- lapply(c("hap1", "hap2", "both"), function(x) { subset(adf, adf[,hapcolumn]==x) })
    names(dfbyhap) <- c("hap1", "hap2", "both")

    # copy "both" haplotypes into dfs for hap1 or hap2, so those will be complete 
    dfbyhap[["hap1"]] <- rbind(dfbyhap[["hap1"]], dfbyhap[["both"]])
    dfbyhap[["hap2"]] <- rbind(dfbyhap[["hap2"]], dfbyhap[["both"]])
                               
    grbyhap <- GRangesList(lapply(c("hap1", "hap2"), function(y) {
                          x <- dfbyhap[[y]]
                          tmp <- GRanges(seqnames=x$chrom, ranges=IRanges(start=x$winstart, end=x$winend - 1))
                          tmp <- reduce(tmp)
                          if(length(tmp) >0) { tmp$id <- y }
                          return(tmp)
                                     }))
    
    ## use unlist to combine GRangesList elements into one GRanges object
    return(unlist(grbyhap))
}


