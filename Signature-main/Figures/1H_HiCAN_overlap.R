#_____Load required packages_____________________________________________________
options(scipen=999)
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(ggVennDiagram)
  library(VennDiagram)
})


#_____Read in arguments_________________________________________________________
args = commandArgs(trailingOnly = TRUE)

hican_nhub <- read.table(args[1], header = T)
hican_shub <- read.table(args[2], header = T)
sign <- read.table(args[3], header = T)
out <- args[4]
ann_file <- args[5]

#_____Prepare data______________________________________________________________

# filter to appropriate cell types
hican_nhub <- hican_nhub %>% select(-c("GM23248", "WI38_raf"))
hican_shub <- hican_shub %>% select(-c("GM23248", "WI38_raf"))

# split up interaction ID information into new columns
sign$ID <- sub("B", "\\.B", as.character(sign$ID))
sign <- sign %>% separate(ID, sep = "\\.", into = c("chrA", "st1", "end1","chrB","st2","end2"), remove = FALSE)
sign$chrA <- gsub("A", "", sign$chrA)
sign$chrB <- gsub("B", "", sign$chrB)



#_____Convert HiCAN res to Signature res________________________________________

# convert hican_shub #
df <- hican_shub

for (ct in colnames(hican_shub)){
  
  df <- df %>% separate(ct, sep = "\\.", into = c("chrA", "st1", "end1"), remove = FALSE)
  
  # convert
  for (r in 1:nrow(df)) {
    
    if (grepl("\\.",as.numeric(df$st1[r])/1000000)){
      
      df$st1[r] <- floor(as.numeric(df$st1[r])/1000000)*1000000
      
    } else if (grepl("\\.",as.numeric(df$end1[r])/1000000)){
      
      df$end1[r] <- ceiling(as.numeric(df$end1[r])/1000000)*1000000
      
    }
    
  }
  
  # replace 500kb loci with new 1mb loci
  coln <- grep(ct, colnames(df))
  df[coln] <- paste(df$chrA, df$st1, df$end1, sep = ".")
  
  # remove these columns for next cell type
  df <- df %>% select(-c("chrA", "st1", "end1"))

}

hican_shub <- df
rm(df)


# convert hican_nhub #
df <- hican_nhub

for (ct in colnames(hican_nhub)){
  
  df <- df %>% separate(ct, sep = "\\.", into = c("chrA", "st1", "end1"), remove = FALSE)
  
  # convert
  for (r in 1:nrow(df)) {
    
    if (grepl("\\.",as.numeric(df$st1[r])/1000000)){
      
      df$st1[r] <- floor(as.numeric(df$st1[r])/1000000)*1000000
      
    } else if (grepl("\\.",as.numeric(df$end1[r])/1000000)){
      
      df$end1[r] <- ceiling(as.numeric(df$end1[r])/1000000)*1000000
      
    }
    
  }
  
  # replace 500kb loci with new 1mb loci
  coln <- grep(ct, colnames(df))
  df[coln] <- paste(df$chrA, df$st1, df$end1, sep = ".")
  
  # remove these columns for next cell type
  df <- df %>% select(-c("chrA", "st1", "end1"))
  
}

hican_nhub <- df
rm(df)



#_____Create interaction loci lists_____________________________________________

# HiCAN speckle
hican_shub_list <- unique(c(hican_shub[,1],hican_shub[,2],hican_shub[,3],hican_shub[,4],hican_shub[,5],hican_shub[,6]))

# HiCAN nucleolus
hican_nhub_list <- unique(c(hican_nhub[,1],hican_nhub[,2],hican_nhub[,3],hican_nhub[,4],hican_nhub[,5],hican_nhub[,6]))

# Signature
sign <- sign[!rowSums(is.na(sign)) == length(colnames(sign))-7, ]   # remove interactions with only NAs (first 7 columns are related to ID)
sign$bin1 <- paste(sign$chrA, sign$st1, sign$end1, sep = ".")
sign$bin2 <- paste(sign$chrB, sign$st2, sign$end2, sep = ".")
sign_list <- c(sign$bin1,sign$bin2)
sign_list <- unique(sign_list)



#_____Produce Venn Diagram______________________________________________________

gglist <- list("signature" = sign_list,
               "hican_speckle" = hican_shub_list,
               "hican_nucleolus" = hican_nhub_list)


vd <- ggVennDiagram(gglist, set_size = 3, label_alpha = 0, label_size = 4, label = "count") +
              labs(title = "Signature (q<0.05) vs HiCAN", fill = "1MB loci") +
              scale_fill_distiller(palette = "Reds", direction = 1) +
              scale_color_brewer(palette = "Greys") +
              theme(plot.title = element_text(face="bold", colour="black", size=14, hjust=0.5, vjust=3),
                    plot.subtitle = element_text(face="plain", colour="black", size=12, hjust=0.5),
                    legend.title = element_text(face="plain", colour="black", size=10, hjust=0.5),
                    legend.text = element_text(face="plain", colour="black", size=8),
                    legend.position = NULL)

ggsave(paste0(out,"/HiCAN_VD.pdf"), width = 8, height = 6)


# get list of unique bins to extract genes later
all_overlaps <- calculate.overlap(gglist) 

total <- lengths(all_overlaps)

ggg <- c(rep("all",total[1]),
         rep("sp+si",total[2]),
         rep("nu+si",total[3]),
         rep("sp+nu",total[4]),
         rep("si",total[5]),
         rep("sp",total[6]),
         rep("nu",total[7]))

bbb <- c(all_overlaps[[1]],
         all_overlaps[[2]],
         all_overlaps[[3]],
         all_overlaps[[4]],
         all_overlaps[[5]],
         all_overlaps[[6]],
         all_overlaps[[7]])

genes_df <- data.frame("bins" = bbb,
                       "group" = ggg)



#_____Annotate with genes_______________________________________________________

dat <- read.table(ann_file, sep = "\t", header = TRUE)

dat$bin_ID <- gsub(":",".",dat$bin_ID)
dat$bin_ID <- gsub("-",".",dat$bin_ID)

genes_df$gene_IDs <- "NA"

for (b in 1:nrow(genes_df)){
  
  gbi <- which(dat$bin_ID == genes_df$bins[b])
  
  if(length(gbi) > 0){genes_df$gene_IDs[b] <- dat$gene_IDs[gbi]}

}

# export
write.table(genes_df, paste0(out,"/overlap_list_annotated_geneIDs.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




################################################################################
#     Repeat analysis for q<0.0001
################################################################################
rm(list = ls()) # clear environment

#_____Prepare data______________________________________________________________

# load data
hican_nhub <- read.table(args[1], header = T)
hican_shub <- read.table(args[2], header = T)
sign <- read.table(args[3], header = T)

# filter to appropriate cell types
hican_nhub <- hican_nhub %>% select(-c("GM23248", "WI38_raf"))
hican_shub <- hican_shub %>% select(-c("GM23248", "WI38_raf"))

# filter to q < 0.0001
for (c in 2:ncol(sign)){
  
  index <- which(sign[,c] > 0.0001)
  sign[index,c] <- NA
  
}
sign <- sign[!rowSums(is.na(sign)) == length(colnames(sign))-1, ]   # remove interactions with only NAs (first column is ID)


# split up interaction ID information into new columns
sign$ID <- sub("B", "\\.B", as.character(sign$ID))
sign <- sign %>% separate(ID, sep = "\\.", into = c("chrA", "st1", "end1","chrB","st2","end2"), remove = FALSE)
sign$chrA <- gsub("A", "", sign$chrA)
sign$chrB <- gsub("B", "", sign$chrB)



#_____Convert HiCAN res to Signature res________________________________________

# convert hican_shub #
df <- hican_shub

for (ct in colnames(hican_shub)){
  
  df <- df %>% separate(ct, sep = "\\.", into = c("chrA", "st1", "end1"), remove = FALSE)
  
  # convert
  for (r in 1:nrow(df)) {
    
    if (grepl("\\.",as.numeric(df$st1[r])/1000000)){
      
      df$st1[r] <- floor(as.numeric(df$st1[r])/1000000)*1000000
      
    } else if (grepl("\\.",as.numeric(df$end1[r])/1000000)){
      
      df$end1[r] <- ceiling(as.numeric(df$end1[r])/1000000)*1000000
      
    }
    
  }
  
  # replace 500kb loci with new 1mb loci
  coln <- grep(ct, colnames(df))
  df[coln] <- paste(df$chrA, df$st1, df$end1, sep = ".")
  
  # remove these columns for next cell type
  df <- df %>% select(-c("chrA", "st1", "end1"))
  
}

hican_shub <- df
rm(df)


# convert hican_nhub #
df <- hican_nhub

for (ct in colnames(hican_nhub)){
  
  df <- df %>% separate(ct, sep = "\\.", into = c("chrA", "st1", "end1"), remove = FALSE)
  
  # convert
  for (r in 1:nrow(df)) {
    
    if (grepl("\\.",as.numeric(df$st1[r])/1000000)){
      
      df$st1[r] <- floor(as.numeric(df$st1[r])/1000000)*1000000
      
    } else if (grepl("\\.",as.numeric(df$end1[r])/1000000)){
      
      df$end1[r] <- ceiling(as.numeric(df$end1[r])/1000000)*1000000
      
    }
    
  }
  
  # replace 500kb loci with new 1mb loci
  coln <- grep(ct, colnames(df))
  df[coln] <- paste(df$chrA, df$st1, df$end1, sep = ".")
  
  # remove these columns for next cell type
  df <- df %>% select(-c("chrA", "st1", "end1"))
  
}

hican_nhub <- df
rm(df)



#_____Create interaction loci lists_____________________________________________

# HiCAN speckle
hican_shub_list <- unique(c(hican_shub[,1],hican_shub[,2],hican_shub[,3],hican_shub[,4],hican_shub[,5],hican_shub[,6]))

# HiCAN nucleolus
hican_nhub_list <- unique(c(hican_nhub[,1],hican_nhub[,2],hican_nhub[,3],hican_nhub[,4],hican_nhub[,5],hican_nhub[,6]))

# Signature
sign$bin1 <- paste(sign$chrA, sign$st1, sign$end1, sep = ".")
sign$bin2 <- paste(sign$chrB, sign$st2, sign$end2, sep = ".")
sign_list <- c(sign$bin1,sign$bin2)
sign_list <- unique(sign_list)



#_____Produce Venn Diagram______________________________________________________

gglist <- list("signature" = sign_list,
               "hican_speckle" = hican_shub_list,
               "hican_nucleolus" = hican_nhub_list)


vd <- ggVennDiagram(gglist, set_size = 3, label_alpha = 0, label_size = 4, label = "count") +
  labs(title = "Signature (q<0.0001) vs HiCAN", fill = "1MB loci") +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_color_brewer(palette = "Greys") +
  theme(plot.title = element_text(face="bold", colour="black", size=14, hjust=0.5, vjust=3),
        plot.subtitle = element_text(face="plain", colour="black", size=12, hjust=0.5),
        legend.title = element_text(face="plain", colour="black", size=10, hjust=0.5),
        legend.text = element_text(face="plain", colour="black", size=8),
        legend.position = NULL)

ggsave(paste0(out,"/HiCAN_VD.0001.pdf"), width = 8, height = 6)


# get list of unique bins to extract genes later
all_overlaps <- calculate.overlap(gglist) 

total <- lengths(all_overlaps)

ggg <- c(rep("all",total[1]),
         rep("sp+si",total[2]),
         rep("nu+si",total[3]),
         rep("sp+nu",total[4]),
         rep("si",total[5]),
         rep("sp",total[6]),
         rep("nu",total[7]))

bbb <- c(all_overlaps[[1]],
         all_overlaps[[2]],
         all_overlaps[[3]],
         all_overlaps[[4]],
         all_overlaps[[5]],
         all_overlaps[[6]],
         all_overlaps[[7]])

genes_df <- data.frame("bins" = bbb,
                       "group" = ggg)


#_____Annotate with genes_______________________________________________________

dat <- read.table(ann_file, sep = "\t", header = TRUE)

dat$bin_ID <- gsub(":",".",dat$bin_ID)
dat$bin_ID <- gsub("-",".",dat$bin_ID)

genes_df$gene_IDs <- "NA"

for (b in 1:nrow(genes_df)){
  
  gbi <- which(dat$bin_ID == genes_df$bins[b])
  
  if(length(gbi) > 0){genes_df$gene_IDs[b] <- dat$gene_IDs[gbi]}
  
}

# export
write.table(genes_df, paste0(out,"/overlap_list_annotated_geneIDs.0001.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
