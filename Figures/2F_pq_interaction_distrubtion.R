suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
options(scipen = 999)


#______Read in arguments________________________________________________________
# running with a scheduler #

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

datafile <- args[1]
out <- args[2]


#_____get data ready____________________________________________________________

# read in input data (pos sig qval)
dat <- read.table(datafile, header = T)


# get cell count
dat$count <- NA 
for (i in 1:nrow(dat)){dat$count[i] <- sum(!is.na(dat[i,2:ncol(dat)]))}

# keep only ID and count
dat <- dat[,c(1,64)]

# split up ID
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
dat$ID <- sub("B", "\\.B", as.character(dat$ID))
dat <- dat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
dat$chrA <- gsub("A", "", dat$chrA)
dat$chrB <- gsub("B", "", dat$chrB)



#_____annotate with p and q arm_________________________________________________

# make vertical
dat1 <- dat[,c(2:4,8)]
colnames(dat1) <- c("chr","st","end","count")
dat2 <- dat[,5:8]
colnames(dat2) <- c("chr","st","end","count")
data <- rbind(dat1,dat2)

# set up files for conversion
data$arm <- NA
data$st_conv <- NA
data$end_conv <- NA



#_____linear conversion to uniform scale________________________________________

chrInf <- data.frame( chrom = c("chr1","chr2","chr3","chr4","chr5","chr6",
                                "chr7","chrX","chr8","chr9","chr11","chr10",
                                "chr12","chr13","chr14","chr15","chr16","chr17",
                                "chr18","chr20","chr19","chrY","chr22","chr21"),
                      centromere = c(123433028,93163954.5,92119691,50396643.5,
                                     48181167.5,59340548,59998520,60601375,
                                     44672548,44381224.5,52656433.5,40067375.5, 
                                     36021355,17234414.5,17131293,18609622.5,
                                     37197922.5,24826902.5,18148357.5,27788943.5,
                                     25888610.5,10440387.5,14475908.5,11941022),
                      size = c(248956422,242193529,198295559,190214555,
                               181538259,170805979,159345973,156040895,145138636,
                               138394717,135086622,133797422,133275309,
                               114364328,107043718,101991189,90338345,
                               83257441,80373285,64444167,58617616,57227415,
                               50818468,46709983))

data$st <- as.numeric(data$st)
data$end <- as.numeric(data$end)

for (i in 1:nrow(data)){
  
  print(paste0("Index ",i))
  # get centromere for that row's bin
  centr <- chrInf$centromere[which(chrInf$chrom == data$chr[i])]
  max <- chrInf$size[which(chrInf$chrom == data$chr[i])]
  
  # assign pq arm, variables, and calculate new values
  if (data$st[i] <= centr){
    
    # new_max limit calculated from 4.2% of unified length
    data$arm[i] <- "p_arm"
    new_min = 0
    new_max = 958
    old_min = 0
    old_max = centr
    
    # calculate new length based on p or q arm
    #start
    value_to_change = data$st[i]
    new_value = ( (value_to_change - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min
    data$st_conv[i] <- new_value
    #end
    value_to_change = data$end[i]
    new_value = ( (value_to_change - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min
    data$end_conv[i] <- new_value
    
    
  } else {
    
    # new_min limit calculated from 4.2% of unified length 
    data$arm[i] <- "q_arm"
    new_min = 1042
    new_max = 2000
    old_min = centr
    old_max = max
    
    # calculate new length based on p or q arm
    #start
    value_to_change = data$st[i]
    new_value = ( (value_to_change - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min
    data$st_conv[i] <- new_value
    #end
    value_to_change = data$end[i]
    new_value = ( (value_to_change - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min
    data$end_conv[i] <- new_value
    
    }
}
  
 
# re-organize columns and prepare for graph
pqdist <- data[,c("chr","st_conv","end_conv","count")]
pqdist[,2:3] <- sapply(pqdist[,2:3], round)

# fix chrY linearization to match small centromere
fix <- which(pqdist$end_conv < 1042 & pqdist$end_conv > 958)
pqdist[fix,2:3] <- round(sapply(pqdist[fix,2:3], function(x) {x/1.055}))



#_____plotting__________________________________________________________________

# null dataframe
total_count_df <- data.frame(unit=c(0:1999),total=0)

# fill in dataframe
for (i in 1:nrow(pqdist)){
  
  print(paste0("Index ",i))
  bins <- seq(from = pqdist$st_conv[i], to = pqdist$end_conv[i]-1, by = 1)
  
  for (b in 1:length(bins)){
    total_count_df$total[which(total_count_df$unit == bins[b])] <- total_count_df$total[which(total_count_df$unit == bins[b])] + pqdist$count[i]
  }
  
}


# set theme

theme_set(theme_bw() + theme(panel.background = element_rect(fill = "white", colour = NA),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             plot.margin = unit(c(0.25,0.75,0.5,0.75), "cm"),
                             legend.background=element_blank(),
                             legend.title = element_text(vjust = 0.8),
                             legend.position="right")
)


# plot

sp <- (ggplot(total_count_df, aes(x=unit+0.05, y=total))
       + geom_line(linewidth = 0.6)
       + labs(x = "Unified Chromosome", y = "Total Number of Significant Interactions", title = "P/Q Arm Interaction Distrubtion")
       + expand_limits(x = c(0,2000),y=c(0,max(total_count_df$total)+50))
       + scale_y_continuous(expand = c(0, 0))
       + scale_x_continuous(expand = c(0, 0))
       + theme(#title
         plot.title = element_text(size=14),
         #x-axis
         axis.title.x = element_text(size = 14, vjust = -1),
         axis.text.x = element_text(size = 12),
         #y-axis
         axis.title.y = element_text(size = 14, vjust = 4),
         axis.text.y = element_text(size = 12),
         
       )
)

file=paste0(out,"/PQ_Interaction_Distrubtion.pdf")
pdf(file, width = 16, height = 8)
print(sp)
dev.off()
