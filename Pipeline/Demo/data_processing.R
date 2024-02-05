library(dplyr)
library(tidyr)

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
resolution <- 1000000
cells <- c("Astrocyte_Spine",
           "H9hESC_day00_Zhang")

CELL <- cells[as.numeric(args[[1]])]

int_table <- read.delim(sprintf("./%s/%s.pairs.res1000000.cool.txt",CELL,CELL), header = F)

cat(sprintf("Cell type = %s", CELL), sep="\n")
colnames(int_table) <- c("chrA", "stA", "endA","chrB", "stB", "endB", "rawread", "freq")

int_table <- int_table[!is.na(int_table$freq),]
int_table <- int_table[,c("chrA", "stA", "chrB", "stB", "freq")]

########deleting rows correspond to chrM
int_table <- int_table[!(int_table$chrA == "chrM" | int_table$chrB == "chrM"), ]

int_table$chrA <- sub("chr", "", as.character(int_table$chrA))
int_table$chrB <- sub("chr", "", as.character(int_table$chrB))

resolution <- as.numeric(resolution)

int_table$stA <- as.numeric(int_table$stA)/resolution 
int_table$stB <- as.numeric(int_table$stB)/resolution

int_table$ID_chrA <- paste(int_table$chrA, int_table$stA, sep = "_")
int_table$ID_chrB <- paste(int_table$chrB, int_table$stB, sep = "_")

int_table$freq <- as.numeric(int_table$freq)

int_table <- int_table[,c("ID_chrA", "ID_chrB", "freq")]

#########################################################
write.table(int_table, sprintf("%s/%s_network.txt", CELL, CELL), row.names = F)
