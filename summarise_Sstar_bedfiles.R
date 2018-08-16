library(data.table)
library(ggplot)
library(dplyr)
library(tidyr)


NA18522 <- fread(input = 'gunzip -c ~/DATALab/SimulatedDemographic/Sstar/chr1_variable_ref/Output_YRI_1/SstarSigFiles/YRI_1_chr_1_start_ALL.NA18522.sstar_sig.bed.gz', col.names=c('msp_ID','start','end'))


NA19223 <- fread(input = 'gunzip -c ~/DATALab/SimulatedDemographic/Sstar/chr1_variable_ref/Output_YRI_1/SstarSigFiles/YRI_1_chr_1_start_ALL.NA19223.sstar_sig.bed.gz', col.names=c('msp_ID','start','end'))

NA19239 <- fread(input = 'gunzip -c ~/DATALab/SimulatedDemographic/Sstar/chr1_variable_ref/Output_YRI_1/SstarSigFiles/YRI_1_chr_1_start_ALL.NA19239.sstar_sig.bed.gz', col.names=c('msp_ID','start','end'))


NA18522 <- NA18522[,kb:=(end-start)/1000] %>% group_by(msp_ID) %>% summarise(YRI_kb=sum(kb)) %>% as.data.table()
NA19223 <- NA19223[,kb:=(end-start)/1000] %>% group_by(msp_ID) %>% summarise(YRI_kb=sum(kb)) %>% as.data.table()
NA19239 <- NA19239[,kb:=(end-start)/1000] %>% group_by(msp_ID) %>% summarise(YRI_kb=sum(kb)) %>% as.data.table()

dt <- left_join(NA18522, NA19223, by=c('msp_ID')) %>% left_join(NA19239,) %>% as.data.table() 

as.data.table(t(apply(X = dt, MARGIN = 1, FUN = summary.fn )))



summary.fn = function(x){
  vector <- c( as.numeric(x[["YRI_kb.x"]]), as.numeric(x[["YRI_kb.y"]]), as.numeric(x[["YRI_kb"]]) ) 
  sd <- sd( vector )
  length <- length(vector)
  x[["mean"]] <- round(mean(vector), digits=2)
  x[["se"]] <- round(sd/sqrt(length),digits = 2)
  x[["ci"]] <- paste0(round(t.test(vector)$conf.int[[1]], digits = 0),':', round(t.test(vector)$conf.int[[2]], digits = 0))
  return(x[c("msp_ID","YRI_kb.x", "YRI_kb.y", "YRI_kb","mean","se","ci")])
}






