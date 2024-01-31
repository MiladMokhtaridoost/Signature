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




#______Get the sig data ready___________________________________________________

print("## read in file")
dat <- read.table(zdat_file, header = TRUE)

# sort data and determine which bins are 'significant'
zscore_df_sig <- c()
zscore_df_sig$ID <- dat$ID
zscore_df_sig <- as.data.frame(zscore_df_sig)
zscore_df_sig$count <- NA
zscore_df_sig$NAs <- NA
zscore_df_sig$percent <- NA

total <- ncol(dat)-1

for (i in 1:nrow(zscore_df_sig)){

  count=0
  for (j in 2:ncol(dat)) {
    if(is.na(dat[i,j]) == FALSE){
      if (dat[i,j] >= 1.959) { count=count+1 }
    }
  }

  zscore_df_sig$count[i] <- count
  zscore_df_sig$NAs[i] <- sum(is.na(dat[i,]))
  zscore_df_sig$percent[i] <- ((zscore_df_sig$count[i]) / (total - zscore_df_sig$NAs[i])) * 100
  zscore_df_sig$percent[i] <- as.numeric(zscore_df_sig$percent[i])
  cat(sprintf("index=%s",i), sep = "\n")

}


# split up interaction ID information into new columns
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
zscore_df_sig$ID <- sub("B", "\\.B", as.character(zscore_df_sig$ID))
zscore_df_sig <- zscore_df_sig %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)

# remove A and B from chrom names
zscore_df_sig$chrA <- gsub("A", "", zscore_df_sig$chrA)
zscore_df_sig$chrB <- gsub("B", "", zscore_df_sig$chrB)

# grab only columns we want
zscore_df_sig = subset(zscore_df_sig, select = -c(1,4,7))

# change columns to numeric
zscore_df_sig$st1 <- as.numeric(zscore_df_sig$st1)
zscore_df_sig$st2 <- as.numeric(zscore_df_sig$st2)

# scale genomic positions by 1Mb
zscore_df_sig$st1 <- zscore_df_sig$st1/1000000
zscore_df_sig$st2 <- zscore_df_sig$st2/1000000

# add column that contains both chromosome numbers of interest
zscore_df_sig$chrs <- paste0(zscore_df_sig$chrA,zscore_df_sig$chrB)

head(zscore_df_sig,4L)


#save for the cluster later
saveRDS(zscore_df_sig, paste0(out,"/zscore_df_sig.rds"))
