################################################################################
# NOTE: THIS SCRIPT IS A PRE-PROCESSING STEP AND ONLY SAVES AN RDS OBJECT
################################################################################

#______Read in arguments________________________________________________________
# running with a scheduler #

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

zdat_file <- args[1]
out <- args[2]




#______Load required packages___________________________________________________

library(tidyr)
library(dplyr)
library(pracma)
library(ggplot2)
library(ggpubr)

options(scipen = 999) #turning off scientific notation




#______Get the data imported and ready__________________________________________

print("## read in file")
dat <- read.table(zdat_file, header = TRUE)

print("## tidy data")
# gather/verticalize data
zdat = gather(dat, key = "cell", value = "zscore", 2:length(colnames(dat)))

# split up interaction ID information into new columns
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
zdat$ID <- sub("B", "\\.B", as.character(zdat$ID))
zdat <- zdat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)

# remove A and B from chrom names
zdat$chrA <- gsub("A", "", zdat$chrA)
zdat$chrB <- gsub("B", "", zdat$chrB)
zscore_df <- zdat

# grab only columns we want
zscore_df = subset(zscore_df, select = -c(1,4,7))

# change columns to numeric
zscore_df$st1 <- as.numeric(zscore_df$st1)
zscore_df$st2 <- as.numeric(zscore_df$st2)

# scale genomic positions by 1Mb
zscore_df$st1 <- zscore_df$st1/1000000
zscore_df$st2 <- zscore_df$st2/1000000

# add column that contains both chromosome numbers of interest
zscore_df$chrs <- paste0(zscore_df$chrA,zscore_df$chrB)

head(zscore_df,4L)

#save for the cluster later
saveRDS(zscore_df, paste0(out,"/zscore_df.rds"))

