library(dplyr) 
library(purrr)

PATHWAY <- "../Demo"

cells <- c("Astrocyte_Spine",
           "H9hESC_day00_Zhang")

read_function <- function(i){
  
  read.table(sprintf("%s/%s/%s_network.txt",PATHWAY,cells[i],cells[i]), header = T)
  
}

print("reading in files ...")
l_dfs <- lapply(1:length(cells), read_function)

# merge all batches together
print("merging ...")
mdat <- reduce(l_dfs, full_join, by = c("ID_chrA","ID_chrB"))

for(i in 3:ncol(mdat)){
  mdat[,i] <- as.numeric(as.character(mdat[,i])) 
}

colnames(mdat) <- c("ID_chrA", "ID_chrB", cells)

mdat$average <- rowMeans(mdat[, -(1:2)], na.rm = TRUE)

avrg <- mdat[,c("ID_chrA", "ID_chrB", "average")]

write.table(avrg, file = sprintf("%s/average_1MB_network.txt",PATHWAY), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

