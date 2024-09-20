# Copyright © 2024 fabio-affaticati
# tools.r

diff_abundance <- function(count_data, meta_data, form){

  require(multidiffabundance)
  require(tidyr)
  
  D <- mda.create(count_data, meta_data, form)
  
  df1 <- mda.corncob(D)
  df2 <- mda.limma(D)
  df3 <- mda.lmclr(D)
  #df4 <- mda.maaslin2(D)
  df5 <- mda.deseq2(D)
  #df6 <- mda.aldex2(D)
  #df7 <- mda.ancom(D)

  
  return(rbind(df1$res, df2$res, df3$res, df5$res)) #df4$res))
}
