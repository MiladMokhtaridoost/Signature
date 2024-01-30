#______Read in arguments________________________________________________________

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

data <- args[1]
outpath <- args[2] 

#______Get the data imported and ready_________________________________________

row_sums <- rowSums(data[, -1], na.rm = TRUE) 
data$Total <- row_sums
print(data)
data_subset <- data[, c(1, ncol(data))] 
print(data_subset)

write.table(data_subset, file = sprintf("%s/diploid_merged.trans.1000000_iced.sorted_total.txt", outpath), quote=FALSE, row.names = FALSE, col.names = TRUE)
