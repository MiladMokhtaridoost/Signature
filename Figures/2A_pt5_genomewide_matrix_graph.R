################################################################################
# Generate genome-wide heatmap of 1vsAll
################################################################################

#______Read in arguments________________________________________________________
# running with a scheduler #

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

rds1 <- args[1]
rds2 <- args[2]
rds3 <- args[3]
rds4 <- args[4]
rds5 <- args[5]
rds6 <- args[6]
out <- args[7]



#______Load required packages___________________________________________________
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(pracma))
suppressMessages(library(ggplot2))
suppressMessages(library(harrypotter))
suppressMessages(library(ggpubr))
suppressMessages(library(utils))
options(scipen = 999) #turning off scientific notation



#______Set ggplot graph theme___________________________________________________

theme_set(theme_bw() + theme(panel.background = element_rect(fill = "white", colour = NA),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.spacing = unit(0, "lines"),
                             plot.title = element_text(hjust = 0.5, size = 16),
                             plot.margin = unit(c(0,0,0,0), "cm"),
                             #legend.key = element_blank(),
                             legend.background=element_blank(),
                             #legend background
                             legend.key = element_rect(fill = NA),
                             legend.title = element_text(vjust = 0.8),
                             legend.position="none",
                             #remove x-axis
                             axis.title.x = element_blank(),
                             axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             #remove y-axis
                             axis.title.y = element_blank(),
                             axis.text.y = element_blank(),
                             axis.ticks.y = element_blank()
)
)

sigcolour <- c("#FFFFFF","#d60202")



#______Get the data imported and ready__________________________________________
## ran this section in the cluster and saved as RDS object 

print("Read in RDS object: zscore_df")
zscore_df <- readRDS(rds1)
head(zscore_df,4L)



#______split 'all' dataframe into male and female_______________________________
## ran this section in the cluster and saved as RDS objects

print("Read in RDS object: zscore_df_X")
zscore_df_X <- readRDS(rds2)
head(zscore_df_X,4L)

print("Read in RDS object: zscore_df_Y")
zscore_df_Y <- readRDS(rds3)
head(zscore_df_Y,4L)




#______Get the sig data ready___________________________________________________
## ran these sections in the cluster and saved as RDS objects 

print("Read in RDS object: zscore_df_sig")
zscore_df_sig <- readRDS(rds4)
head(zscore_df_sig,4L)


print("Read in RDS object: zscore_df_sig_X")
zscore_df_sig_X <- readRDS(rds5)
head(zscore_df_sig_X,4L)


print("Read in RDS object: zscore_df_sig_Y")
zscore_df_sig_Y <- readRDS(rds6)
head(zscore_df_sig_Y,4L)




#______get chromosomal indexes__________________________________________________

chrpairs <- unique(zscore_df$chrs)

chrpairs_x <- chrpairs[grepl(glob2rx("*chrX*"), chrpairs) == TRUE]
chrpairs_all <- chrpairs[grepl(glob2rx("*chrX*"), chrpairs) == FALSE]



#______generate trans plots for chr1-22 and chrY________________________________
print("Generate individual plots for autosomes and chrY")

for (i in 1:length(chrpairs_all)){
  
  cat(sprintf("pair %s/%s: %s", i, length(chrpairs_all),chrpairs_all[i]), "\n")
  
  # all data #
  # select for only your pairwise interaction
  pwdata <- zscore_df %>% filter(chrs == chrpairs_all[i])
  # get mean zscore for graph
  mzscore_df_sum <- pwdata %>% group_by(st1,st2) %>% dplyr::summarize(mzscore = mean(zscore, na.rm = TRUE))
  
  # sig data #
  # select for only your pairwise interaction
  pwdata_sig <- zscore_df_sig %>% filter(chrs == chrpairs_all[i])
  
  
  # heat maps #
  x=unique(gsub("chr","",pwdata$chrA))
  y=unique(gsub("chr","",pwdata$chrB))
  if (x == "X") { x = 23 }
  if (y == "X") { y = 23 }
  if (y == "Y") { y = 24 }
  if (x == "Y") { x = 24 }
  x=as.numeric(x)
  y=as.numeric(y)
  
  if (x < y){
    
    # all
    hm1 <- (ggplot(mzscore_df_sum, aes(x=st1, y=st2, fill = mzscore))
            + geom_tile(aes(fill = mzscore))
            + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score", na.value = "white")
            + labs(x = paste0("Chromosome ", x),y = paste0("Chromosome ", y))
            + scale_y_continuous(expand = c(0, 0))
            + scale_x_continuous(expand = c(0, 0))
    )
    assign(paste0("hm_",chrpairs_all[i],"_all"), hm1) ; rm(hm1)
    
    # sig  
    hm2 <- (ggplot(pwdata_sig, aes(y=st1, x=st2, fill = percent))
            + geom_tile(aes(fill = percent))
            + scale_fill_gradientn(colors = sigcolour, name = "Mean z-score", na.value = "white")
            + labs(x = paste0("Chromosome ", y),y = paste0("Chromosome ", x))
            + scale_x_continuous(expand = c(0, 0))
            + scale_y_reverse(expand = c(0, 0))
    )
    assign(paste0("hm_",chrpairs_all[i],"_sig"), hm2) ; rm(hm2)
    
    
  } else {
    
    
    # all
    hm1 <- (ggplot(mzscore_df_sum, aes(y=st1, x=st2, fill = mzscore))
            + geom_tile(aes(fill = mzscore))
            + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score", na.value = "white")
            + labs(x = paste0("Chromosome ", y),y = paste0("Chromosome ", x))
            + scale_y_continuous(expand = c(0, 0))
            + scale_x_continuous(expand = c(0, 0))
    )
    assign(paste0("hm_",chrpairs_all[i],"_all"), hm1) ; rm(hm1)
    
    # sig  
    hm2 <- (ggplot(pwdata_sig, aes(x=st1, y=st2, fill = percent))
            + geom_tile(aes(fill = percent))
            + scale_fill_gradientn(colors = sigcolour, name = "Mean z-score", na.value = "white")
            + labs(x = paste0("Chromosome ", x),y = paste0("Chromosome ", y))
            + scale_x_continuous(expand = c(0, 0))
            + scale_y_reverse(expand = c(0, 0))
    )
    assign(paste0("hm_",chrpairs_all[i],"_sig"), hm2) ; rm(hm2)  
    
  }
}




#______generate trans plots for chrX (female XX samples)________________________
print("Generate individual plots for chrX (XX datasets)")

for (i in 1:length(chrpairs_x)){
  
  cat(sprintf("pair %s/%s: %s", i, length(chrpairs_x),chrpairs_x[i]), "\n")
  
  # all data #
  # select for only your pairwise interaction
  pwdata <- zscore_df_X %>% filter(chrs == chrpairs_x[i])
  # get mean zscore for graph
  mzscore_df_sum <- pwdata %>% group_by(st1,st2) %>% dplyr::summarize(mzscore = mean(zscore, na.rm = TRUE))
  
  # sig data #
  # select for only your pairwise interaction
  pwdata_sig <- zscore_df_sig_X %>% filter(chrs == chrpairs_x[i])
  
  
  # heat maps #
  x=unique(gsub("chr","",pwdata$chrA))
  y=unique(gsub("chr","",pwdata$chrB))
  if (x == "X") { x = 23 }
  if (y == "X") { y = 23 }
  if (y == "Y") { y = 24 }
  if (x == "Y") { x = 24 }
  x=as.numeric(x)
  y=as.numeric(y)
  
  if (x < y){
    
    # all
    hm1 <- (ggplot(mzscore_df_sum, aes(x=st1, y=st2, fill = mzscore))
            + geom_tile(aes(fill = mzscore))
            + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score", na.value = "white")
            + labs(x = paste0("Chromosome ", x),y = paste0("Chromosome ", y))
            + scale_y_continuous(expand = c(0, 0))
            + scale_x_continuous(expand = c(0, 0))
    )
    assign(paste0("hm_",chrpairs_x[i],"_all_X"), hm1) ; rm(hm1)
    
    # sig  
    hm2 <- (ggplot(pwdata_sig, aes(y=st1, x=st2, fill = percent))
            + geom_tile(aes(fill = percent))
            + scale_fill_gradientn(colors = sigcolour, name = "Mean z-score", na.value = "white")
            + labs(x = paste0("Chromosome ", y),y = paste0("Chromosome ", x))
            + scale_x_continuous(expand = c(0, 0))
            + scale_y_reverse(expand = c(0, 0))
    )
    assign(paste0("hm_",chrpairs_x[i],"_sig_X"), hm2) ; rm(hm2)
    
    
  } else {
    
    
    # all
    hm1 <- (ggplot(mzscore_df_sum, aes(y=st1, x=st2, fill = mzscore))
            + geom_tile(aes(fill = mzscore))
            + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score", na.value = "white")
            + labs(x = paste0("Chromosome ", y),y = paste0("Chromosome ", x))
            + scale_y_continuous(expand = c(0, 0))
            + scale_x_continuous(expand = c(0, 0))
    )
    assign(paste0("hm_",chrpairs_x[i],"_all_X"), hm1) ; rm(hm1)
    
    # sig  
    hm2 <- (ggplot(pwdata_sig, aes(x=st1, y=st2, fill = percent))
            + geom_tile(aes(fill = percent))
            + scale_fill_gradientn(colors = sigcolour, name = "Mean z-score", na.value = "white")
            + labs(x = paste0("Chromosome ", x),y = paste0("Chromosome ", y))
            + scale_x_continuous(expand = c(0, 0))
            + scale_y_reverse(expand = c(0, 0))
    )
    assign(paste0("hm_",chrpairs_x[i],"_sig_X"), hm2) ; rm(hm2)  
    
  }
}



#______generate trans plots for chrX (male XY samples)__________________________
print("Generate individual plots for chrX (XY datasets)")

for (i in 1:length(chrpairs_x)){
  
  cat(sprintf("pair %s/%s: %s", i, length(chrpairs_x),chrpairs_x[i]), "\n")
  
  # all data #
  # select for only your pairwise interaction
  pwdata <- zscore_df_Y %>% filter(chrs == chrpairs_x[i])
  # get mean zscore for graph
  mzscore_df_sum <- pwdata %>% group_by(st1,st2) %>% dplyr::summarize(mzscore = mean(zscore, na.rm = TRUE))
  
  # sig data #
  # select for only your pairwise interaction
  pwdata_sig <- zscore_df_sig_Y %>% filter(chrs == chrpairs_x[i])
  
  
  # heat maps #
  x=unique(gsub("chr","",pwdata$chrA))
  y=unique(gsub("chr","",pwdata$chrB))
  if (x == "X") { x = 23 }
  if (y == "X") { y = 23 }
  if (y == "Y") { y = 24 }
  if (x == "Y") { x = 24 }
  x=as.numeric(x)
  y=as.numeric(y)
  
  if (x < y){
    
    # all
    hm1 <- (ggplot(mzscore_df_sum, aes(x=st1, y=st2, fill = mzscore))
            + geom_tile(aes(fill = mzscore))
            + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score", na.value = "white")
            + labs(x = paste0("Chromosome ", x),y = paste0("Chromosome ", y))
            + scale_y_continuous(expand = c(0, 0))
            + scale_x_continuous(expand = c(0, 0))
    )
    assign(paste0("hm_",chrpairs_x[i],"_all_Y"), hm1) ; rm(hm1)
    
    # sig  
    hm2 <- (ggplot(pwdata_sig, aes(y=st1, x=st2, fill = percent))
            + geom_tile(aes(fill = percent))
            + scale_fill_gradientn(colors = sigcolour, name = "Mean z-score", na.value = "white")
            + labs(x = paste0("Chromosome ", y),y = paste0("Chromosome ", x))
            + scale_x_continuous(expand = c(0, 0))
            + scale_y_reverse(expand = c(0, 0))
    )
    assign(paste0("hm_",chrpairs_x[i],"_sig_Y"), hm2) ; rm(hm2)
    
    
  } else {
    
    
    # all
    hm1 <- (ggplot(mzscore_df_sum, aes(y=st1, x=st2, fill = mzscore))
            + geom_tile(aes(fill = mzscore))
            + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score", na.value = "white")
            + labs(x = paste0("Chromosome ", y),y = paste0("Chromosome ", x))
            + scale_y_continuous(expand = c(0, 0))
            + scale_x_continuous(expand = c(0, 0))
    )
    assign(paste0("hm_",chrpairs_x[i],"_all_Y"), hm1) ; rm(hm1)
    
    # sig  
    hm2 <- (ggplot(pwdata_sig, aes(x=st1, y=st2, fill = percent))
            + geom_tile(aes(fill = percent))
            + scale_fill_gradientn(colors = sigcolour, name = "Mean z-score", na.value = "white")
            + labs(x = paste0("Chromosome ", x),y = paste0("Chromosome ", y))
            + scale_x_continuous(expand = c(0, 0))
            + scale_y_reverse(expand = c(0, 0))
    )
    assign(paste0("hm_",chrpairs_x[i],"_sig_Y"), hm2) ; rm(hm2)  
    
  }
}




#______generate fake cis plot___________________________________________________
mzscore_df_sum <- data.frame(st1=1:50,st2=1:50,mzscore=c("NA"))
hm_cis <- (ggplot(mzscore_df_sum, aes(y=st1, x=st2, fill = mzscore))
           + geom_blank()
           + scale_y_continuous(expand = c(0, 0))
           + scale_x_continuous(expand = c(0, 0)))



#______facet all plots together into grid_______________________________________
print("Facet all heatmaps together into genomewide matrix")

final_hm <- ggarrange(
  #chr1
  hm_cis,hm_chr1chr2_sig,hm_chr1chr3_sig,hm_chr1chr4_sig,hm_chr1chr5_sig,hm_chr1chr6_sig,hm_chr1chr7_sig,hm_chr1chr8_sig,hm_chr1chr9_sig,hm_chr1chr10_sig,hm_chr1chr11_sig,hm_chr1chr12_sig,hm_chr1chr13_sig,hm_chr1chr14_sig,hm_chr1chr15_sig,hm_chr1chr16_sig,hm_chr1chr17_sig,hm_chr1chr18_sig,hm_chr1chr19_sig,hm_chr1chr20_sig,hm_chr1chr21_sig,hm_chr1chr22_sig,hm_chr1chrX_sig_X,hm_chr1chrX_sig_Y,hm_chr1chrY_sig,
  #chr2
  hm_chr1chr2_all,hm_cis,hm_chr2chr3_sig,hm_chr2chr4_sig,hm_chr2chr5_sig,hm_chr2chr6_sig,hm_chr2chr7_sig,hm_chr2chr8_sig,hm_chr2chr9_sig,hm_chr2chr10_sig,hm_chr2chr11_sig,hm_chr2chr12_sig,hm_chr2chr13_sig,hm_chr2chr14_sig,hm_chr2chr15_sig,hm_chr2chr16_sig,hm_chr2chr17_sig,hm_chr2chr18_sig,hm_chr2chr19_sig,hm_chr2chr20_sig,hm_chr2chr21_sig,hm_chr2chr22_sig,hm_chr2chrX_sig_X,hm_chr2chrX_sig_Y,hm_chr2chrY_sig,
  #chr3
  hm_chr1chr3_all,hm_chr2chr3_all,hm_cis,hm_chr3chr4_sig,hm_chr3chr5_sig,hm_chr3chr6_sig,hm_chr3chr7_sig,hm_chr3chr8_sig,hm_chr3chr9_sig,hm_chr3chr10_sig,hm_chr3chr11_sig,hm_chr3chr12_sig,hm_chr3chr13_sig,hm_chr3chr14_sig,hm_chr3chr15_sig,hm_chr3chr16_sig,hm_chr3chr17_sig,hm_chr3chr18_sig,hm_chr3chr19_sig,hm_chr3chr20_sig,hm_chr3chr21_sig,hm_chr3chr22_sig,hm_chr3chrX_sig_X,hm_chr3chrX_sig_Y,hm_chr3chrY_sig,
  #chr4
  hm_chr1chr4_all,hm_chr2chr4_all,hm_chr3chr4_all,hm_cis,hm_chr4chr5_sig,hm_chr4chr6_sig,hm_chr4chr7_sig,hm_chr4chr8_sig,hm_chr4chr9_sig,hm_chr4chr10_sig,hm_chr4chr11_sig,hm_chr4chr12_sig,hm_chr4chr13_sig,hm_chr4chr14_sig,hm_chr4chr15_sig,hm_chr4chr16_sig,hm_chr4chr17_sig,hm_chr4chr18_sig,hm_chr4chr19_sig,hm_chr4chr20_sig,hm_chr4chr21_sig,hm_chr4chr22_sig,hm_chr4chrX_sig_X,hm_chr4chrX_sig_Y,hm_chr4chrY_sig,
  #chr5
  hm_chr1chr5_all,hm_chr2chr5_all,hm_chr3chr5_all,hm_chr4chr5_all,hm_cis,hm_chr5chr6_sig,hm_chr5chr7_sig,hm_chr5chr8_sig,hm_chr5chr9_sig,hm_chr5chr10_sig,hm_chr5chr11_sig,hm_chr5chr12_sig,hm_chr5chr13_sig,hm_chr5chr14_sig,hm_chr5chr15_sig,hm_chr5chr16_sig,hm_chr5chr17_sig,hm_chr5chr18_sig,hm_chr5chr19_sig,hm_chr5chr20_sig,hm_chr5chr21_sig,hm_chr5chr22_sig,hm_chr5chrX_sig_X,hm_chr5chrX_sig_Y,hm_chr5chrY_sig,
  #chr6
  hm_chr1chr6_all,hm_chr2chr6_all,hm_chr3chr6_all,hm_chr4chr6_all,hm_chr5chr6_all,hm_cis,hm_chr6chr7_sig,hm_chr6chr8_sig,hm_chr6chr9_sig,hm_chr6chr10_sig,hm_chr6chr11_sig,hm_chr6chr12_sig,hm_chr6chr13_sig,hm_chr6chr14_sig,hm_chr6chr15_sig,hm_chr6chr16_sig,hm_chr6chr17_sig,hm_chr6chr18_sig,hm_chr6chr19_sig,hm_chr6chr20_sig,hm_chr6chr21_sig,hm_chr6chr22_sig,hm_chr6chrX_sig_X,hm_chr6chrX_sig_Y,hm_chr6chrY_sig,
  #chr7
  hm_chr1chr7_all,hm_chr2chr7_all,hm_chr3chr7_all,hm_chr4chr7_all,hm_chr5chr7_all,hm_chr6chr7_all,hm_cis,hm_chr7chr8_sig,hm_chr7chr9_sig,hm_chr7chr10_sig,hm_chr7chr11_sig,hm_chr7chr12_sig,hm_chr7chr13_sig,hm_chr7chr14_sig,hm_chr7chr15_sig,hm_chr7chr16_sig,hm_chr7chr17_sig,hm_chr7chr18_sig,hm_chr7chr19_sig,hm_chr7chr20_sig,hm_chr7chr21_sig,hm_chr7chr22_sig,hm_chr7chrX_sig_X,hm_chr7chrX_sig_Y,hm_chr7chrY_sig,
  #chr8
  hm_chr1chr8_all,hm_chr2chr8_all,hm_chr3chr8_all,hm_chr4chr8_all,hm_chr5chr8_all,hm_chr6chr8_all,hm_chr7chr8_all,hm_cis,hm_chr8chr9_sig,hm_chr8chr10_sig,hm_chr8chr11_sig,hm_chr8chr12_sig,hm_chr8chr13_sig,hm_chr8chr14_sig,hm_chr8chr15_sig,hm_chr8chr16_sig,hm_chr8chr17_sig,hm_chr8chr18_sig,hm_chr8chr19_sig,hm_chr8chr20_sig,hm_chr8chr21_sig,hm_chr8chr22_sig,hm_chrXchr8_sig_X,hm_chrXchr8_sig_Y,hm_chr8chrY_sig,
  #chr9
  hm_chr1chr9_all,hm_chr2chr9_all,hm_chr3chr9_all,hm_chr4chr9_all,hm_chr5chr9_all,hm_chr6chr9_all,hm_chr7chr9_all,hm_chr8chr9_all,hm_cis,hm_chr9chr10_sig,hm_chr9chr11_sig,hm_chr9chr12_sig,hm_chr9chr13_sig,hm_chr9chr14_sig,hm_chr9chr15_sig,hm_chr9chr16_sig,hm_chr9chr17_sig,hm_chr9chr18_sig,hm_chr9chr19_sig,hm_chr9chr20_sig,hm_chr9chr21_sig,hm_chr9chr22_sig,hm_chrXchr9_sig_X,hm_chrXchr9_sig_Y,hm_chr9chrY_sig,
  #chr10
  hm_chr1chr10_all,hm_chr2chr10_all,hm_chr3chr10_all,hm_chr4chr10_all,hm_chr5chr10_all,hm_chr6chr10_all,hm_chr7chr10_all,hm_chr8chr10_all,hm_chr9chr10_all,hm_cis,hm_chr11chr10_sig,hm_chr10chr12_sig,hm_chr10chr13_sig,hm_chr10chr14_sig,hm_chr10chr15_sig,hm_chr10chr16_sig,hm_chr10chr17_sig,hm_chr10chr18_sig,hm_chr10chr19_sig,hm_chr10chr20_sig,hm_chr10chr21_sig,hm_chr10chr22_sig,hm_chrXchr10_sig_X,hm_chrXchr10_sig_Y,hm_chr10chrY_sig,
  #chr11
  hm_chr1chr11_all,hm_chr2chr11_all,hm_chr3chr11_all,hm_chr4chr11_all,hm_chr5chr11_all,hm_chr6chr11_all,hm_chr7chr11_all,hm_chr8chr11_all,hm_chr9chr11_all,hm_chr11chr10_all,hm_cis,hm_chr11chr12_sig,hm_chr11chr13_sig,hm_chr11chr14_sig,hm_chr11chr15_sig,hm_chr11chr16_sig,hm_chr11chr17_sig,hm_chr11chr18_sig,hm_chr11chr19_sig,hm_chr11chr20_sig,hm_chr11chr21_sig,hm_chr11chr22_sig,hm_chrXchr11_sig_X,hm_chrXchr11_sig_Y,hm_chr11chrY_sig,
  #chr12
  hm_chr1chr12_all,hm_chr2chr12_all,hm_chr3chr12_all,hm_chr4chr12_all,hm_chr5chr12_all,hm_chr6chr12_all,hm_chr7chr12_all,hm_chr8chr12_all,hm_chr9chr12_all,hm_chr10chr12_all,hm_chr11chr12_all,hm_cis,hm_chr12chr13_sig,hm_chr12chr14_sig,hm_chr12chr15_sig,hm_chr12chr16_sig,hm_chr12chr17_sig,hm_chr12chr18_sig,hm_chr12chr19_sig,hm_chr12chr20_sig,hm_chr12chr21_sig,hm_chr12chr22_sig,hm_chrXchr12_sig_X,hm_chrXchr12_sig_Y,hm_chr12chrY_sig,
  #chr13
  hm_chr1chr13_all,hm_chr2chr13_all,hm_chr3chr13_all,hm_chr4chr13_all,hm_chr5chr13_all,hm_chr6chr13_all,hm_chr7chr13_all,hm_chr8chr13_all,hm_chr9chr13_all,hm_chr10chr13_all,hm_chr11chr13_all,hm_chr12chr13_all,hm_cis,hm_chr13chr14_sig,hm_chr13chr15_sig,hm_chr13chr16_sig,hm_chr13chr17_sig,hm_chr13chr18_sig,hm_chr13chr19_sig,hm_chr13chr20_sig,hm_chr13chr21_sig,hm_chr13chr22_sig,hm_chrXchr13_sig_X,hm_chrXchr13_sig_Y,hm_chr13chrY_sig,
  #chr14
  hm_chr1chr14_all,hm_chr2chr14_all,hm_chr3chr14_all,hm_chr4chr14_all,hm_chr5chr14_all,hm_chr6chr14_all,hm_chr7chr14_all,hm_chr8chr14_all,hm_chr9chr14_all,hm_chr10chr14_all,hm_chr11chr14_all,hm_chr12chr14_all,hm_chr13chr14_all,hm_cis,hm_chr14chr15_sig,hm_chr14chr16_sig,hm_chr14chr17_sig,hm_chr14chr18_sig,hm_chr14chr19_sig,hm_chr14chr20_sig,hm_chr14chr21_sig,hm_chr14chr22_sig,hm_chrXchr14_sig_X,hm_chrXchr14_sig_Y,hm_chr14chrY_sig,
  #chr15
  hm_chr1chr15_all,hm_chr2chr15_all,hm_chr3chr15_all,hm_chr4chr15_all,hm_chr5chr15_all,hm_chr6chr15_all,hm_chr7chr15_all,hm_chr8chr15_all,hm_chr9chr15_all,hm_chr10chr15_all,hm_chr11chr15_all,hm_chr12chr15_all,hm_chr13chr15_all,hm_chr14chr15_all,hm_cis,hm_chr15chr16_sig,hm_chr15chr17_sig,hm_chr15chr18_sig,hm_chr15chr19_sig,hm_chr15chr20_sig,hm_chr15chr21_sig,hm_chr15chr22_sig,hm_chrXchr15_sig_X,hm_chrXchr15_sig_Y,hm_chr15chrY_sig,
  #chr16
  hm_chr1chr16_all,hm_chr2chr16_all,hm_chr3chr16_all,hm_chr4chr16_all,hm_chr5chr16_all,hm_chr6chr16_all,hm_chr7chr16_all,hm_chr8chr16_all,hm_chr9chr16_all,hm_chr10chr16_all,hm_chr11chr16_all,hm_chr12chr16_all,hm_chr13chr16_all,hm_chr14chr16_all,hm_chr15chr16_all,hm_cis,hm_chr16chr17_sig,hm_chr16chr18_sig,hm_chr16chr19_sig,hm_chr16chr20_sig,hm_chr16chr21_sig,hm_chr16chr22_sig,hm_chrXchr16_sig_X,hm_chrXchr16_sig_Y,hm_chr16chrY_sig,
  #chr17
  hm_chr1chr17_all,hm_chr2chr17_all,hm_chr3chr17_all,hm_chr4chr17_all,hm_chr5chr17_all,hm_chr6chr17_all,hm_chr7chr17_all,hm_chr8chr17_all,hm_chr9chr17_all,hm_chr10chr17_all,hm_chr11chr17_all,hm_chr12chr17_all,hm_chr13chr17_all,hm_chr14chr17_all,hm_chr15chr17_all,hm_chr16chr17_all,hm_cis,hm_chr17chr18_sig,hm_chr17chr19_sig,hm_chr17chr20_sig,hm_chr17chr21_sig,hm_chr17chr22_sig,hm_chrXchr17_sig_X,hm_chrXchr17_sig_Y,hm_chr17chrY_sig,
  #chr18
  hm_chr1chr18_all,hm_chr2chr18_all,hm_chr3chr18_all,hm_chr4chr18_all,hm_chr5chr18_all,hm_chr6chr18_all,hm_chr7chr18_all,hm_chr8chr18_all,hm_chr9chr18_all,hm_chr10chr18_all,hm_chr11chr18_all,hm_chr12chr18_all,hm_chr13chr18_all,hm_chr14chr18_all,hm_chr15chr18_all,hm_chr16chr18_all,hm_chr17chr18_all,hm_cis,hm_chr18chr19_sig,hm_chr18chr20_sig,hm_chr18chr21_sig,hm_chr18chr22_sig,hm_chrXchr18_sig_X,hm_chrXchr18_sig_Y,hm_chr18chrY_sig,
  #chr19
  hm_chr1chr19_all,hm_chr2chr19_all,hm_chr3chr19_all,hm_chr4chr19_all,hm_chr5chr19_all,hm_chr6chr19_all,hm_chr7chr19_all,hm_chr8chr19_all,hm_chr9chr19_all,hm_chr10chr19_all,hm_chr11chr19_all,hm_chr12chr19_all,hm_chr13chr19_all,hm_chr14chr19_all,hm_chr15chr19_all,hm_chr16chr19_all,hm_chr17chr19_all,hm_chr18chr19_all,hm_cis,hm_chr20chr19_sig,hm_chr19chr21_sig,hm_chr19chr22_sig,hm_chrXchr19_sig_X,hm_chrXchr19_sig_Y,hm_chr19chrY_sig,
  #chr20
  hm_chr1chr20_all,hm_chr2chr20_all,hm_chr3chr20_all,hm_chr4chr20_all,hm_chr5chr20_all,hm_chr6chr20_all,hm_chr7chr20_all,hm_chr8chr20_all,hm_chr9chr20_all,hm_chr10chr20_all,hm_chr11chr20_all,hm_chr12chr20_all,hm_chr13chr20_all,hm_chr14chr20_all,hm_chr15chr20_all,hm_chr16chr20_all,hm_chr17chr20_all,hm_chr18chr20_all,hm_chr20chr19_all,hm_cis,hm_chr20chr21_sig,hm_chr20chr22_sig,hm_chrXchr20_sig_X,hm_chrXchr20_sig_Y,hm_chr20chrY_sig,
  #chr21
  hm_chr1chr21_all,hm_chr2chr21_all,hm_chr3chr21_all,hm_chr4chr21_all,hm_chr5chr21_all,hm_chr6chr21_all,hm_chr7chr21_all,hm_chr8chr21_all,hm_chr9chr21_all,hm_chr10chr21_all,hm_chr11chr21_all,hm_chr12chr21_all,hm_chr13chr21_all,hm_chr14chr21_all,hm_chr15chr21_all,hm_chr16chr21_all,hm_chr17chr21_all,hm_chr18chr21_all,hm_chr19chr21_all,hm_chr20chr21_all,hm_cis,hm_chr22chr21_sig,hm_chrXchr21_sig_X,hm_chrXchr21_sig_Y,hm_chrYchr21_sig,
  #chr22
  hm_chr1chr22_all,hm_chr2chr22_all,hm_chr3chr22_all,hm_chr4chr22_all,hm_chr5chr22_all,hm_chr6chr22_all,hm_chr7chr22_all,hm_chr8chr22_all,hm_chr9chr22_all,hm_chr10chr22_all,hm_chr11chr22_all,hm_chr12chr22_all,hm_chr13chr22_all,hm_chr14chr22_all,hm_chr15chr22_all,hm_chr16chr22_all,hm_chr17chr22_all,hm_chr18chr22_all,hm_chr19chr22_all,hm_chr20chr22_all,hm_chr22chr21_all,hm_cis,hm_chrXchr22_sig_X,hm_chrXchr22_sig_Y,hm_chrYchr22_sig,
  
  #chrX (XX)
  hm_chr1chrX_all_X,hm_chr2chrX_all_X,hm_chr3chrX_all_X,hm_chr4chrX_all_X,hm_chr5chrX_all_X,hm_chr6chrX_all_X,hm_chr7chrX_all_X,hm_chrXchr8_all_X,hm_chrXchr9_all_X,hm_chrXchr10_all_X,hm_chrXchr11_all_X,hm_chrXchr12_all_X,hm_chrXchr13_all_X,hm_chrXchr14_all_X,hm_chrXchr15_all_X,hm_chrXchr16_all_X,hm_chrXchr17_all_X,hm_chrXchr18_all_X,hm_chrXchr19_all_X,hm_chrXchr20_all_X,hm_chrXchr21_all_X,hm_chrXchr22_all_X,hm_cis,hm_cis,hm_chrXchrY_sig_X,
  #chrX (XY)
  hm_chr1chrX_all_Y,hm_chr2chrX_all_Y,hm_chr3chrX_all_Y,hm_chr4chrX_all_Y,hm_chr5chrX_all_Y,hm_chr6chrX_all_Y,hm_chr7chrX_all_Y,hm_chrXchr8_all_Y,hm_chrXchr9_all_Y,hm_chrXchr10_all_Y,hm_chrXchr11_all_Y,hm_chrXchr12_all_Y,hm_chrXchr13_all_Y,hm_chrXchr14_all_Y,hm_chrXchr15_all_Y,hm_chrXchr16_all_Y,hm_chrXchr17_all_Y,hm_chrXchr18_all_Y,hm_chrXchr19_all_Y,hm_chrXchr20_all_Y,hm_chrXchr21_all_Y,hm_chrXchr22_all_Y,hm_cis,hm_cis,hm_chrXchrY_sig_Y,
  
  #chrY
  hm_chr1chrY_all,hm_chr2chrY_all,hm_chr3chrY_all,hm_chr4chrY_all,hm_chr5chrY_all,hm_chr6chrY_all,hm_chr7chrY_all,hm_chr8chrY_all,hm_chr9chrY_all,hm_chr10chrY_all,hm_chr11chrY_all,hm_chr12chrY_all,hm_chr13chrY_all,hm_chr14chrY_all,hm_chr15chrY_all,hm_chr16chrY_all,hm_chr17chrY_all,hm_chr18chrY_all,hm_chr19chrY_all,hm_chr20chrY_all,hm_chrYchr21_all,hm_chrYchr22_all,hm_chrXchrY_all_X,hm_chrXchrY_all_Y,hm_cis,
  
  #number of columns and rows
  ncol=25, nrow=25
  
)


print("Exporting ...")

filename <- paste0(out,"/genomewide_matrix_graph_Xsep.pdf")
pdf(filename, width = 18, height = 18)
print(final_hm)
dev.off()
print("done")



