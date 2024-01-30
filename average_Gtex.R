gtex <- read.table("Z:/Jordan/signature/CD_gtex/output/gtex_filtered_trans1MB_all_bins.txt", header = T)

# Split the column by ":"
split_values <- strsplit(gtex$bin_ID, ":", fixed = TRUE)

# Create two new columns for the split values
gtex$chr <- sapply(split_values, function(x) x[1])
gtex$StEnd <- sapply(split_values, function(x) x[2])

# Split the column by "-"
split_values2 <- strsplit(gtex$StEnd, "-", fixed = TRUE)

# Create two new columns for the split values
gtex$st <- sapply(split_values2, function(x) x[1])
gtex$end <- sapply(split_values2, function(x) x[2])

###
gtex$st <- as.numeric(gtex$st)/1000000

gtex$chr <- sub("chr", "", gtex$chr)

gtex$ID <- paste(gtex$chr, gtex$st, sep = "_")


gtex$gtex_avg[which(is.na(gtex$gtex_avg)==TRUE)] = 0
#######
gtex <- gtex[,c("ID", "gtex_avg")]
write.csv(gtex, "Z:/Milad/Community_detection/Gephi/average_gtex_raw.csv", row.names = F)

#####
gtex$gtex_avg <- log2(gtex$gtex_avg+1)

#######
max <- max(gtex$gtex_avg)
min <- min(gtex$gtex_avg)

gtex$normalized_gtex <- (gtex$gtex_avg - min) / (max - min)

gtex <- gtex[,c("ID", "normalized_gtex")]

write.csv(gtex, "Z:/Milad/Community_detection/Gephi/average_gtex_filtered_log2.csv", row.names = F)
