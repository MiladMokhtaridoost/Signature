options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
tmp <- args[1]
#########################################################################

#read in data

interactions<-read.table(tmp)
##pairwise
colnames(interactions)=c("index","id","chr1","st1","end1","chr2","st2","end2","dist","freq","mean","sd","zscore","pvalue")
interactions$interID <- do.call(paste0, interactions[c("chr1","st1","end1", "chr2", "st2", "end2")])
sDat <- interactions[order(interactions$interID, abs(interactions$zscore) ), ]
uDat <- sDat[ !duplicated(sDat$interID), ]

##################save file##########################
write.table(uDat,gsub(" ", "", paste(tmp,".tmp.zscores.trans.lowest.score.txt")), sep="\t", quote=FALSE, row.name=FALSE,col.name=FALSE)
