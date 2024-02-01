#______Read in arguments________________________________________________________
options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

qvalue_all <- read.table(args[1], header = T, sep = "\t")
qvalue_pos <- read.table(args[2], header = T, sep = "\t")
genomic_scale <- read.table(args[3], header = T, sep = "\t")

#_____Load required packages_____________________________________________________
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(gridExtra)

#____Data preparation_____________________________________________________________

#qvalue neg of IMR90_Dixon dataset
qvalue_all_IMR90_Dixon <- qvalue_all[,c("ID","IMR90_Dixon")]
qvalue_pos_IMR90_Dixon <- qvalue_pos[,c("ID","IMR90_Dixon")]
qvalue_pos_IMR90_Dixon_not_na <- qvalue_pos_IMR90_Dixon[!is.na(qvalue_pos_IMR90_Dixon$IMR90_Dixon),]
qvalue_all_IMR90_Dixon_not_na <- qvalue_all_IMR90_Dixon[!is.na(qvalue_all_IMR90_Dixon$IMR90_Dixon),]

rows_not_pos_IMR90_Dixon <- !(qvalue_all_IMR90_Dixon_not_na$ID %in% qvalue_pos_IMR90_Dixon_not_na$ID)
qvalue_neg_IMR90_Dixon<- qvalue_all_IMR90_Dixon_not_na[ rows_not_pos_IMR90_Dixon,]

#significant qvalue pos and neg of IMR90_Dixon dataset
qvalue_neg_IMR90_Dixon_sig <- qvalue_neg_IMR90_Dixon[qvalue_neg_IMR90_Dixon$IMR90_Dixon < 0.05,]
qvalue_pos_IMR90_Dixon_sig <- qvalue_pos_IMR90_Dixon_not_na[qvalue_pos_IMR90_Dixon_not_na$IMR90_Dixon < 0.05,]

#qvalue neg of IMR90_Rao dataset
qvalue_all_IMR90_Rao <- qvalue_all[,c("ID", "IMR90_Rao")]
qvalue_pos_IMR90_Rao <- qvalue_pos[,c("ID", "IMR90_Rao")]
qvalue_all_IMR90_Rao_not_na <- qvalue_all_IMR90_Rao[!is.na(qvalue_all_IMR90_Rao$IMR90_Rao),]
qvalue_pos_IMR90_Rao_not_na <- qvalue_pos_IMR90_Rao[!is.na(qvalue_pos_IMR90_Rao$IMR90_Rao),]
rows_not_posIMR90_Rao <- !(qvalue_all_IMR90_Rao_not_na$ID %in% qvalue_pos_IMR90_Rao_not_na$ID)
qvalue_neg_IMR90_Rao<- qvalue_all_IMR90_Rao_not_na[ rows_not_posIMR90_Rao,]

#significant qvalue pos and neg for IMR90_Dixon
qvalue_neg_IMR90_Rao_sig <- qvalue_neg_IMR90_Rao[qvalue_neg_IMR90_Rao$IMR90_Rao < 0.05,]
qvalue_pos_IMR90_Rao_sig <- qvalue_pos_IMR90_Rao_not_na[qvalue_pos_IMR90_Rao_not_na$IMR90_Rao < 0.05,]

#___________________________________________split columns_______________________________________________________

###################################################
#split columns for table "qvalue_pos_IMR90_Rao_sig" 
###################################################

colnm <- c("chrA", "st1", "end1", "chrB", "st2", "end2")
qvalue_pos_IMR90_Rao_sig$ID <- sub("B", "\\.B", as.character(qvalue_pos_IMR90_Rao_sig$ID))
qvalue_pos_IMR90_Rao_sig <- qvalue_pos_IMR90_Rao_sig %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
head(qvalue_pos_IMR90_Rao_sig)
qvalue_pos_IMR90_Rao_sig$chrA <-  gsub("A","",qvalue_pos_IMR90_Rao_sig$chrA)
qvalue_pos_IMR90_Rao_sig$chrB <-  gsub("B","",qvalue_pos_IMR90_Rao_sig$chrB)

####################################################
#split columns for table "qvalue_neg_IMR90_Rao_sig" 
####################################################

colnm <- c("chrA", "st1", "end1", "chrB", "st2", "end2")
qvalue_neg_IMR90_Rao_sig$ID <- sub("B", "\\.B", as.character(qvalue_neg_IMR90_Rao_sig$ID))
#qvalue_pos_IMR90_Dixon_not_na$ID <- sub("A", "\\.A", as.character(qvalue_pos_IMR90_Dixon_not_na$ID))
qvalue_neg_IMR90_Rao_sig <- qvalue_neg_IMR90_Rao_sig %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
head(qvalue_neg_IMR90_Rao_sig)
qvalue_neg_IMR90_Rao_sig$chrA <-  gsub("A","",qvalue_neg_IMR90_Rao_sig$chrA)
qvalue_neg_IMR90_Rao_sig$chrB <-  gsub("B","",qvalue_neg_IMR90_Rao_sig$chrB)

#####################################################
#split columns for table "qvalue_neg_IMR90_Dixon_sig" 
#####################################################

colnm <- c("chrA", "st1", "end1", "chrB", "st2", "end2")
qvalue_neg_IMR90_Dixon_sig$ID <- sub("B", "\\.B", as.character(qvalue_neg_IMR90_Dixon_sig$ID))
qvalue_neg_IMR90_Dixon_sig <- qvalue_neg_IMR90_Dixon_sig %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
head(qvalue_neg_IMR90_Dixon_sig)
qvalue_neg_IMR90_Dixon_sig$chrA <-  gsub("A","",qvalue_neg_IMR90_Dixon_sig$chrA)
qvalue_neg_IMR90_Dixon_sig$chrB <-  gsub("B","",qvalue_neg_IMR90_Dixon_sig$chrB)

######################################################
#split columns for table "qvalue_pos_IMR90_Dixon_sig" 
######################################################

colnm <- c("chrA", "st1", "end1", "chrB", "st2", "end2")
qvalue_pos_IMR90_Dixon_sig$ID <- sub("B", "\\.B", as.character(qvalue_pos_IMR90_Dixon_sig$ID))
qvalue_pos_IMR90_Dixon_sig <- qvalue_pos_IMR90_Dixon_sig %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
head(qvalue_pos_IMR90_Dixon_sig)
qvalue_pos_IMR90_Dixon_sig$chrA <-  gsub("A","",qvalue_pos_IMR90_Dixon_sig$chrA)
qvalue_pos_IMR90_Dixon_sig$chrB <-  gsub("B","",qvalue_pos_IMR90_Dixon_sig$chrB)

#___________ merge pos datasets______________________________________________________________

merged_pos <- merge(qvalue_pos_IMR90_Dixon_sig,qvalue_pos_IMR90_Rao_sig,by=c("ID","chrA","st1", "end1", "chrB", "st2", "end2"), all = TRUE)
#__________merge neg datasets__________________________________________________________________

merged_neg <- merge(qvalue_neg_IMR90_Dixon_sig,qvalue_neg_IMR90_Rao_sig,by=c("ID","chrA","st1", "end1", "chrB", "st2", "end2"), all = TRUE)
#___________read merfish data and merged_pos and neg____________________________

genomic_scale_split <-separate(genomic_scale,genomic.coordinate,into =c("chr","start","end"),sep = "[:-]")
genomic_scale_split <- as.data.frame(genomic_scale_split)

class(genomic_scale_split$start)
genomic_scale_split$start <- as.numeric(genomic_scale_split$start)
min(genomic_scale_split$start)

genomic_scale_split$index <- seq_len((nrow(genomic_scale_split)))
merged_pos$index <- seq_len((nrow(merged_pos)))
merfish <- genomic_scale_split[,c(4,5,6)]

##########distance function ##############
calculate_distance <- function(x1, y1, z1, x2, y2, z2) {
  distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
  return(distance)
}

#_____________________________________________find distance value for merged_pos_________________________________________
summary_sig <- list()
summary_genomic_scale_split <- list()

genomic_scale_split$end <- as.numeric(genomic_scale_split$end)
genomic_scale_split$start <- as.numeric(genomic_scale_split$start)
genomic_scale_split$z.nm. <- as.numeric(genomic_scale_split$z.nm.)
genomic_scale_split$x.nm. <- as.numeric(genomic_scale_split$x.nm)
genomic_scale_split$y.nm. <- as.numeric(genomic_scale_split$y.nm.)

merged_pos$MerDis <- NA 

for (i in 1:nrow(merged_pos)) {
  
  
  cat(sprintf("row_num = %s", i), sep="\n")
  
  ## (chrA, st1, end1)
  
  condition1_anch <- which(genomic_scale_split$chr == merged_pos$chrA[i] &
                 genomic_scale_split$start >= merged_pos$st1[i] &
                 genomic_scale_split$start <= merged_pos$end1[i])
  
  
  
  condition2_anch <- which(genomic_scale_split$chr == merged_pos$chrA[i] &
                 genomic_scale_split$end >= merged_pos$st1[i] &
                 genomic_scale_split$end <= merged_pos$end1[i])
  
  
  
  anch <- union(condition1_anch, condition2_anch)
  
  mean_x_anchor <- mean(genomic_scale_split$x.nm.[anch], na.rm = T)
  #non_na_x_value_anchor  <- x_value_anchor[!is.na(x_value_anchor)]
  #mean_X_anchor <- mean(non_na_x_value_anchor, na.rm = T)
  
  mean_y_anchor <- mean(genomic_scale_split$y.nm.[anch], na.rm = T)
  
  
  mean_z_anchor <- mean(genomic_scale_split$z.nm.[anch], na.rm = T)
  
  
  ## (chrB, st2, end2) 
  
  condition1_tar <- which(genomic_scale_split$chr == merged_pos$chrB[i] &
                 genomic_scale_split$start >= merged_pos$st2[i] &
                 genomic_scale_split$start <= merged_pos$end2[i]) 
  
  condition2_tar <- which(genomic_scale_split$chr == merged_pos$chrB[i] &
                 genomic_scale_split$end >= merged_pos$st2[i] &
                 genomic_scale_split$end <= merged_pos$end2[i])
  
  tar <- union(condition1_tar, condition2_tar)
  
  mean_x_target <- mean(genomic_scale_split$x.nm.[tar], na.rm = T)
  
  
  mean_y_target <- mean(genomic_scale_split$y.nm.[tar], na.rm = T)
  
  
  mean_z_target <- mean(genomic_scale_split$z.nm.[tar], na.rm = T)
  
  ###
  
  #distance
  
  
  if(!any(is.na(c(mean_x_anchor,mean_y_anchor,mean_z_anchor,mean_z_target,mean_x_target,mean_y_target,mean_z_target)))){
    merged_pos$MerDis[i] <- calculate_distance(mean_x_anchor, mean_y_anchor, mean_z_anchor, mean_x_target, mean_y_target, mean_z_target)
  }
}

#_________________________________________ find distance value for merged_neg_____________________________________________
summary_sig <- list()
summary_genomic_scale_split <- list()

genomic_scale_split$end <- as.numeric(genomic_scale_split$end)
genomic_scale_split$start <- as.numeric(genomic_scale_split$start)
genomic_scale_split$z.nm. <- as.numeric(genomic_scale_split$z.nm.)
genomic_scale_split$x.nm. <- as.numeric(genomic_scale_split$x.nm)
genomic_scale_split$y.nm. <- as.numeric(genomic_scale_split$y.nm.)

merged_neg$MerDis <- NA 

for (i in 1:nrow(merged_neg)) {
  
  
  cat(sprintf("row_num = %s", i), sep="\n")
  
  ## (chrA, st1, end1)
  
  condition_1_anchor_neg <- which(genomic_scale_split$chr == merged_neg$chrA[i] &
                                    genomic_scale_split$start >= merged_neg$st1[i] &
                                    genomic_scale_split$start <= merged_neg$end1[i])
  
  
  
  condition_2_anchor_neg <- which(genomic_scale_split$chr == merged_neg$chrA[i] &
                                    genomic_scale_split$end >= merged_neg$st1[i] &
                                    genomic_scale_split$end <= merged_neg$end1[i])
  
  
  
  union_anchor_neg<- union(condition_1_anchor_neg, condition_2_anchor_neg)
  
  mean_X_anchor_neg <- mean(genomic_scale_split$x.nm.[union_anchor_neg], na.rm = T)
  
  mean_y_anchor_neg <- mean(genomic_scale_split$y.nm.[union_anchor_neg], na.rm = T)
  
  mean_z_anchor_neg <- mean(genomic_scale_split$z.nm.[union_anchor_neg], na.rm = T)
  
  
  ## (chrB, st2, end2) 
  
  condition_1_target_neg <- which(genomic_scale_split$chr == merged_neg$chrB[i] &
                                    genomic_scale_split$start >= merged_neg$st2[i] &
                                    genomic_scale_split$start <= merged_neg$end2[i]) 
  
  condition_1_target_neg <- which(genomic_scale_split$chr == merged_neg$chrB[i] &
                                    genomic_scale_split$end >= merged_neg$st2[i] &
                                    genomic_scale_split$end <= merged_neg$end2[i])
  
  target_neg <- union(condition_1_target_neg, condition_1_target_neg)
  
  mean_X_target_neg <- mean(genomic_scale_split$x.nm.[target_neg], na.rm = T)
  
  mean_y_target_neg <- mean(genomic_scale_split$y.nm.[target_neg], na.rm = T)
  
  
  mean_z_target_neg <- mean(genomic_scale_split$z.nm.[target_neg], na.rm = T)
  
  ###
  
  #distance
  
  
  if(!any(is.na(c(mean_X_anchor_neg,mean_y_anchor_neg,mean_z_anchor_neg,mean_z_target_neg,mean_X_target_neg,mean_y_target_neg,mean_z_target_neg)))){
    merged_neg$MerDis[i] <- calculate_distance(mean_X_anchor_neg, mean_y_anchor_neg, mean_z_anchor_neg, mean_X_target_neg, mean_y_target_neg, mean_z_target_neg)
  }
}

#___________________________________________visualization________________________________________________
merged_pos_with_distance <- merged_pos
merged_neg_with_distance <- merged_neg
pos <- na.omit(merged_pos_with_distance$MerDis)
neg <- na.omit(merged_neg_with_distance$MerDis) 
pos <- as.numeric(pos)
neg <- as.numeric(neg)

length(pos)
pos_50 <-length( pos[pos<=50])
pos_50
pos_100 <- length(pos[pos>50 & pos<=100])
pos_100
pos_150 <- length(pos[pos>100 & pos<=150])
pos_150
pos_200 <- length(pos[pos>150 & pos<=200])
pos_200
pos_250 <- length(pos[pos>200 & pos<=250])
pos_250
pos_300 <- length(pos[pos>250 & pos<=300])
pos_300
pos_350 <- length(pos[pos>300 & pos<=350])
pos_350
pos_400 <- length(pos[ pos>350 &pos<=400])
pos_400


neg_50 <- length(neg[neg<=50])
neg_50
neg_100 <- length(neg[ neg>50 & neg<=100])
neg_100
neg_150 <- length(neg[neg>100 & neg<=150])
neg_150
neg_200 <- length(neg[neg>150 & neg<=200])
neg_200
neg_250 <- length(neg[neg>200 & neg<=250])
neg_250
neg_300 <- length(neg[neg>250 & neg<=300])
neg_300
neg_350 <- length(neg[neg>300 & neg<=350])
neg_350
neg_400 <- length(neg[neg>350 & neg<=400])
neg_400 

#_________________________________________number of interactions in each region (not cumulative)________________________________________________

pos_data <- c(pos_50,pos_100,pos_150,pos_200,pos_250,pos_300,pos_350,pos_400)
neg_data <- c(neg_50,neg_100,neg_150,neg_200,neg_250,neg_300,neg_350,neg_400) 
df  <- data.frame(distance = c(pos_data,neg_data),
                  group = c(50,100,150,200,250,300,350,400,50,100,150,200,250,300,350,400),type = rep(c(rep("pos",8),rep("neg", 8))))


number_plot_not_cumulative<- ggplot(df,aes(x=factor(group),y= distance,fill= type))+
  geom_bar(position = "dodge",stat= "identity")+
  labs(title = "number of positive and negative significant interactions across different distances",
       x = "Distance (Nm)", y= "Number of interactions")+
  scale_fill_manual(values = c("purple","yellow"))+ theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
                                                                theme(plot.title = element_text(hjust = 0.5, size = 16),
                                                                      panel.background = element_rect(fill = "white", colour = NA),
                                                                      panel.grid.major = element_blank(),
                                                                      panel.grid.minor = element_blank(),
                                                                      # panel.spacing = unit(0.25, "lines"),
                                                                      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                                                                      #change the colour of facet label background
                                                                      #strip.background = element_rect(fill = "#E6E1EA"),
                                                                      #remove space between facets
                                                                      # panel.spacing.x=unit(0, "lines"),
                                                                      #legend.key = element_blank(),
                                                                      legend.background=element_blank(),
                                                                      #legend background
                                                                      legend.key = element_rect(fill = NA),
                                                                      legend.title = element_text(vjust = 0.8),
                                                                      legend.position="top",
                                                                      axis.title.x = element_text(size = 14, vjust =-1),
                                                                      axis.text.x = element_text(size = 12),
                                                                      axis.text.y = element_text(size = 12),
                                                                      axis.title.y = element_text(size = 14, vjust = -1)))

filename <- paste0("pos_neg_number_not_cumulative.pdf")
pdf(filename, width = 14, height = 8)
print(number_plot_not_cumulative)
dev.off()
print("done")


#test significance
df_per  <- data.frame(distance = c(pos_data,neg_data),
                      group = c(50,100,150,200,250,300,350,400,50,100,150,200,250,300,350,400),type = rep(c(rep("pos",8),rep("neg", 8))))


pos <- subset(df_per, type == "pos")
neg <- subset(df_per, type == "neg")

wilcox.test(pos$distance,neg$distance, paired = TRUE,exact = FALSE,correct = TRUE)
