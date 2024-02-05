library(dplyr)
library(tidyr)

PATHWAY <- "../Demo"
resolution <- 1000000

cells <- c("Astrocyte_Spine",
           "H9hESC_day00_Zhang")

for 
CELL <- cells[as.numeric(args[[1]])]
CELL <- cells[2]
###DATA_PATH <- "Z:/3D-flow/normalized_data_4DNuc_pipeline/human/CHLA9_Maass"

#int_table <- read.table(sprintf("%s/%s/trans.100000_iced.sorted.txt", DATA_PATH, CELL))
#int_table <- read.table(sprintf("%s/%s/trans.100000_iced.sorted.txt", DATA_PATH, CELL))
int_table <- read.delim(sprintf("%s/%s/%s.pairs.res1000000.cool.txt",DATA_PATH,CELL,CELL), header = F)

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
#write.csv(int_table, sprintf("%s/%s_network.csv", RESULT_PATH, CELL))
write.table(int_table, sprintf("%s/%s_network.txt", RESULT_PATH, CELL), row.names = F)
