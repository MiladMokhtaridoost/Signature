################################################################################
# NOTE: THIS SCRIPT IS THE PRE-PROCESSING STEP AND ONLY SAVES AN RDS OBJECT
################################################################################


#_____Read in arguments_________________________________________________________
args = commandArgs(trailingOnly = TRUE)

zdat_file <- args[1]
pdat_file <- args[2]
output <- args[3]
type <- args[4]


#_____Load required packages____________________________________________________

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))


#_____Get the data imported and ready___________________________________________

print("#read in files")
zdat <- read.table(zdat_file, header = TRUE)
pdat <- read.table(pdat_file, header = TRUE)

print("summary of ALL zscores per cell type")
summary(zdat)
print("summary of ALL pvalues per cell type")
summary(pdat)

#combine zscore and pvalue into one df
zdatL <- gather(zdat, key = "cell", value = "zscore", 2:length(colnames(zdat)))
pdatL <- gather(pdat, key = "cell", value = "pvalue", 2:length(colnames(pdat)))
dat <- merge(zdatL, pdatL, by=c("ID","cell"))


#split up interaction ID information into new columns
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
dat$ID <- sub("B", "\\.B", as.character(dat$ID))
dat2 <- dat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
dat2$chrA <- gsub("A", "", dat2$chrA)
dat2$chrB <- gsub("B", "", dat2$chrB)
hm_zscore_df <- dat2
head(hm_zscore_df)

#change number columns to numeric
hm_zscore_df$st1 <- as.numeric(hm_zscore_df$st1)
hm_zscore_df$end1 <- as.numeric(hm_zscore_df$end1)
hm_zscore_df$st2 <- as.numeric(hm_zscore_df$st2)
hm_zscore_df$end2 <- as.numeric(hm_zscore_df$end2)

#add column that contains both chromosome numbers
hm_zscore_df$chrs <- paste0(hm_zscore_df$chrA,hm_zscore_df$chrB)
head(hm_zscore_df)



# export RDS object 
saveRDS(hm_zscore_df, file=paste0(output,"/prepared_data_",type,".rds"))

