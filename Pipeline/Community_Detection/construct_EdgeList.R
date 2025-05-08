library(dplyr)
library(tidyr)

comms <- read.table("./data/average_netwrok_comms_46.txt", sep = ",", header = T)
comms <- comms %>% separate(ID.name, sep = "\\_", into = c("chr","st"), remove = FALSE)
#comms$chr <- as.numeric(comms$chr)
comms$st <- as.numeric(comms$st)

sorted_comms <- comms[order(comms$chr, comms$st),]

chr_list <- c(1:22, "X", "Y")
edge_list_final <- data.frame()
for(i in 1:length(chr_list)){
  chr_name <- chr_list[i]
  data_chr <- sorted_comms[which(sorted_comms$chr == chr_name),]
 if(nrow(data_chr) > 1){
  edge_list <- matrix(NA, ncol = 2, nrow = nrow(data_chr)-1)
  for(j in 1:(nrow(data_chr)-1)){
  edge_list[j,1] <- data_chr$ID.name[j]
  edge_list[j,2] <- data_chr$ID.name[j+1]
  }
  edge_list_final <- rbind(edge_list_final, edge_list)
  }
}
colnames(edge_list_final) <- c("Source", "Target")

edge_list_final <- edge_list_final %>% separate(Target, sep = "\\_", into = c("chr","st"), remove = FALSE)
edge_list_final <- edge_list_final[,-4]

write.table(edge_list_final, "./data/edgelist_average_netwrok.csv", sep = ",", row.names = F)

