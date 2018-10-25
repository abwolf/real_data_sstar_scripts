library(data.table)
library(ggplot)
library(dplyr)
library(tidyr)


NA18876 <- fread(input = 'gunzip -c ~/DATALab/SimulatedDemographic/Sstar/chr1_variable_ref/Output_YRI_1/SstarSigFiles/NA18876/YRI_1_chr_1_start_ALL.NA18876.sstar_sig.bed.gz', col.names=c('msp_ID','start','end'))


NA19137 <- fread(input = 'gunzip -c ~/DATALab/SimulatedDemographic/Sstar/chr1_variable_ref/Output_YRI_1/SstarSigFiles/NA19137/YRI_1_chr_1_start_ALL.NA19137.sstar_sig.bed.gz', col.names=c('msp_ID','start','end'))

NA19247 <- fread(input = 'gunzip -c ~/DATALab/SimulatedDemographic/Sstar/chr1_variable_ref/Output_YRI_1/SstarSigFiles/NA19247/YRI_1_chr_1_start_ALL.NA19247.sstar_sig.bed.gz', col.names=c('msp_ID','start','end'))


NA18876 <- NA18876[,kb:=(end-start)/1000] %>% group_by(msp_ID) %>% summarise(YRI_kb=sum(kb)) %>% as.data.table()
NA19137 <- NA19137[,kb:=(end-start)/1000] %>% group_by(msp_ID) %>% summarise(YRI_kb=sum(kb)) %>% as.data.table()
NA19247 <- NA19247[,kb:=(end-start)/1000] %>% group_by(msp_ID) %>% summarise(YRI_kb=sum(kb)) %>% as.data.table()

dt <- left_join(NA18876, NA19137, by=c('msp_ID')) %>% left_join(NA19247) %>% as.data.table() 

as.data.table(t(apply(X = dt, MARGIN = 1, FUN = summary_for_ind.fn )))



summary_for_ind.fn = function(x){
  vector <- c( as.numeric(x[["YRI_kb.x"]]), as.numeric(x[["YRI_kb.y"]]), as.numeric(x[["YRI_kb"]]) ) 
  sd <- sd( vector )
  length <- length(vector)
  x[["mean"]] <- round(mean(vector), digits=2)
  x[["se"]] <- round(sd/sqrt(length),digits = 2)
  x[["ci"]] <- paste0(round(t.test(vector)$conf.int[[1]], digits = 0),':', round(t.test(vector)$conf.int[[2]], digits = 0))
  return(x[c("msp_ID","YRI_kb.x", "YRI_kb.y", "YRI_kb","mean","se","ci")])
}






