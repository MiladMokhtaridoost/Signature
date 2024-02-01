#_____Load required packages_____________________________________________________

library(tidyr)
library(dplyr)
library(data.table)

#_____Read in arguments_________________________________________________________

args = commandArgs(trailingOnly = TRUE)

data_path <- args[1]
outpath <- args[2]

#____split column ID____________________________________________________________
data <- read.table(sprintf("%s/signature_trans1vsAll_1MB_merged_zscores.txt"), header= T)
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
data$ID <- sub("B", "\\.B", as.character(data$ID))
data <- data %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
head(data)
data<- data[,-c(1,4,7)]
data

#make 62 files for each data set
for (col in names(data)[5:ncol(data)]) {


  if (file.exists(paste0(cells_pathway, col, ".txt"))) {
    next
  }

  new_table <- data[c( "chrA", "st1", "chrB", "st2", col)]

  file_name <- paste0(col, ".txt")
  write.table(new_table, file = sprintf("%s/%s", outpath, file_name), sep = "\t", row.names = FALSE)
  
  
}

#______________________remove rows with NA in the last column____________________

na.omit_last <- function(df) {
  last_col <- ncol(df)
  na_rows <- is.na(df[, last_col])
  df[!na_rows, ]
}

for (filename in list.files(path = outpath ,pattern = "\\.txt$")) {
  
  df <- read.table(filename, header = TRUE, sep = "\t")
  
  df <- na.omit_last(df)
  
  write.table(df, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
}



