################################################################################
# NOTE: THIS SCRIPT IS A PRE-PROCESSING STEP AND ONLY SAVES AN RDS OBJECT
################################################################################

#______Read in arguments________________________________________________________
# running with a scheduler #

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

zdat_file <- args[1]
out <- args[2]
f_in <- args[3]
m_in <- args[4]

print("## read in files")
m_list <- read.table(m_in)
f_list <- read.table(f_in)
dat <- read.table(zdat_file, header = TRUE)
dat <- dat[grep("chrX", dat$ID), ]




#______Load required packages___________________________________________________

library(tidyr)
library(dplyr)
library(pracma)
library(ggplot2)
library(harrypotter)
options(scipen = 999) #turning off scientific notation




#_______________________________________________________________________________
# sort data and determine which bins are 'significant' for XX samples

# set up null dataframe
zscore_df_sig_X <- c()
zscore_df_sig_X$ID <- dat$ID
zscore_df_sig_X <- as.data.frame(zscore_df_sig_X)
zscore_df_sig_X$count <- NA
zscore_df_sig_X$NAs <- NA
zscore_df_sig_X$percent <- NA

# subset XX female samples
f_list <- c("ID",f_list$V1)
fdat <- dat %>% select(contains(c(f_list))

# get total number of samples included
total <- ncol(fdat)-1

# calculate percent significant
for (i in 1:nrow(zscore_df_sig_X)){
  count=0
  for (j in 2:ncol(fdat)) {
    if(is.na(fdat[i,j]) == FALSE){
      if (fdat[i,j] >= 1.959) { count=count+1 }
    }
  }

  zscore_df_sig_X$count[i] <- count
  zscore_df_sig_X$NAs[i] <- sum(is.na(fdat[i,]))
  zscore_df_sig_X$percent[i] <- ((zscore_df_sig_X$count[i]) / (total - zscore_df_sig_X$NAs[i])) * 100
  zscore_df_sig_X$percent[i] <- as.numeric(zscore_df_sig_X$percent[i])
  cat(sprintf("index=%s",i), sep = "\n")

}

# split up interaction ID information into new columns
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
zscore_df_sig_X$ID <- sub("B", "\\.B", as.character(zscore_df_sig_X$ID))
zscore_df_sig_X <- zscore_df_sig_X %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)

# remove A and B from chrom names
zscore_df_sig_X$chrA <- gsub("A", "", zscore_df_sig_X$chrA)
zscore_df_sig_X$chrB <- gsub("B", "", zscore_df_sig_X$chrB)

# grab only columns we want
zscore_df_sig_X = subset(zscore_df_sig_X, select = -c(1,4,7))

# change columns to numeric
zscore_df_sig_X$st1 <- as.numeric(zscore_df_sig_X$st1)
zscore_df_sig_X$st2 <- as.numeric(zscore_df_sig_X$st2)

# scale genomic positions by 1Mb
zscore_df_sig_X$st1 <- zscore_df_sig_X$st1/1000000
zscore_df_sig_X$st2 <- zscore_df_sig_X$st2/1000000

# add column that contains both chromosome numbers of interest
zscore_df_sig_X$chrs <- paste0(zscore_df_sig_X$chrA,zscore_df_sig_X$chrB)

saveRDS(zscore_df_sig_X, paste0(out,"/zscore_df_sig_X.rds"))
print("zscore_df_sig_X.rds ... SAVED")


#_______________________________________________________________________________
# sort data and determine which bins are 'significant' for XY samples

# set up null dataframe
zscore_df_sig_Y <- c()
zscore_df_sig_Y$ID <- dat$ID
zscore_df_sig_Y <- as.data.frame(zscore_df_sig_Y)
zscore_df_sig_Y$count <- NA
zscore_df_sig_Y$NAs <- NA
zscore_df_sig_Y$percent <- NA

# subset XX female samples
m_list <- c("ID",m_list$V1)
mdat <- dat %>% select(contains(m_list))


# get total number of samples included
total <- ncol(mdat)-1

# calculate percent significant
for (i in 1:nrow(zscore_df_sig_Y)){
  count=0
  for (j in 2:ncol(mdat)) {
    if(is.na(mdat[i,j]) == FALSE){
      if (mdat[i,j] >= 1.959) { count=count+1 }
    }
  }

  zscore_df_sig_Y$count[i] <- count
  zscore_df_sig_Y$NAs[i] <- sum(is.na(mdat[i,]))
  zscore_df_sig_Y$percent[i] <- ((zscore_df_sig_Y$count[i]) / (total - zscore_df_sig_Y$NAs[i])) * 100
  zscore_df_sig_Y$percent[i] <- as.numeric(zscore_df_sig_Y$percent[i])
  cat(sprintf("index=%s",i), sep = "\n")

}

# split up interaction ID information into new columns
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
zscore_df_sig_Y$ID <- sub("B", "\\.B", as.character(zscore_df_sig_Y$ID))
zscore_df_sig_Y <- zscore_df_sig_Y %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)

# remove A and B from chrom names
zscore_df_sig_Y$chrA <- gsub("A", "", zscore_df_sig_Y$chrA)
zscore_df_sig_Y$chrB <- gsub("B", "", zscore_df_sig_Y$chrB)

# grab only columns we want
zscore_df_sig_Y = subset(zscore_df_sig_Y, select = -c(1,4,7))

# change columns to numeric
zscore_df_sig_Y$st1 <- as.numeric(zscore_df_sig_Y$st1)
zscore_df_sig_Y$st2 <- as.numeric(zscore_df_sig_Y$st2)

# scale genomic positions by 1Mb
zscore_df_sig_Y$st1 <- zscore_df_sig_Y$st1/1000000
zscore_df_sig_Y$st2 <- zscore_df_sig_Y$st2/1000000

# add column that contains both chromosome numbers of interest
zscore_df_sig_Y$chrs <- paste0(zscore_df_sig_Y$chrA,zscore_df_sig_Y$chrB)

saveRDS(zscore_df_sig_Y, paste0(out,"/zscore_df_sig_Y.rds"))
print("zscore_df_sig_Y.rds ... SAVED")

