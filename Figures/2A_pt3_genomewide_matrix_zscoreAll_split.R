################################################################################
# NOTE: THIS SCRIPT IS A PRE-PROCESSING STEP AND ONLY SAVES AN RDS OBJECT
################################################################################

#______Read in arguments________________________________________________________
# running with a scheduler #

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

zscr_all <- args[1]
out <- args[2]
f_in <- args[3]
m_in <- args[4]



#______Load required packages___________________________________________________

library(tidyr)
library(dplyr)
library(pracma)
library(ggplot2)
library(ggpubr)
library(utils)
options(scipen = 999)




#______Get the data imported and ready__________________________________________

print("## read in file zscore_df")
zscore_df <- readRDS(zscr_all)
head(zscore_df,4L)


#female XX
f_list <- read.table(f_in)
zscore_df_X <- zscore_df[grep(paste0("(", paste(f_list$V1, collapse=")|("),")"), zscore_df$cell), ]

saveRDS(zscore_df_X, paste0(out,"/zscore_df_X.rds"))
print("zscore_df_X.rds ... SAVED")


#male XY
m_list <- read.table(m_in)
zscore_df_Y <- zscore_df[grep(paste0("(", paste(m_list$V1, collapse=")|("),")"), zscore_df$cell), ]

saveRDS(zscore_df_Y, paste0(out,"/zscore_df_Y.rds"))
print("zscore_df_Y.rds ... SAVED")

