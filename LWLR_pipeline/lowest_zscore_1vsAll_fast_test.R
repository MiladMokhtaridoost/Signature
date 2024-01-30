# choose lowest absolute value of z-score
##########################################################################
# Developer: Jordan Chalmers
##########################################################################


cell = c("Macrophage_patientF_Zhang",
         "Mesoderm_cell_day02_Zhang")


for (i in 1:length(cell)) {
  
  cat(sprintf("Index %s: %s",i,cell[i]),"\n")
  tmp <- sprintf("Z:/Signature/results/output_new/diploid/trans-1vsAll/100KB/%s/delete/tmp.trans.%s.trans-1vsAll.100KB.%s.100000_zscore.all.txt",cell[i],cell[i],cell[i])
  
  print("load file")  
  interactions<-read.table(tmp)
  colnames(interactions)=c("st_anchor" ,"st_target" ,"ID", "chrA", "stA", "endA", "chrB", "stB", "endB", "freq", "target_chr", "estimation", "sd", "zscore", "pvalue","qvalue") 
  
  print("order all zscores")
  #interactions$interID <- do.call(paste0, interactions[c("chrA", "stA", "endA", "chrB", "stB", "endB")])
  #sDat <- interactions[order(interactions$interID, abs(interactions$zscore) ), ]
  sDat <- interactions[order(interactions$ID,interactions$stA,interactions$stB,abs(interactions$zscore)), ]
  
  print("remove absolute highest zscore")
  #uDat <- sDat[ !duplicated(sDat$interID), ]
  index <- seq(from = 1, to = nrow(interactions), by = 2)
  uDat <- sDat[index,]
  
  ##################save file##########################
  print("write output file")
  write.table(uDat,gsub(" ", "", paste(tmp,".tmp.zscores.trans.lowest.score.txt")), sep="\t", quote=FALSE, row.name=FALSE,col.name=FALSE)
  
}
