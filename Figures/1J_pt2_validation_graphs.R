################################################################################
# NOTE: THIS SCRIPT PRODUCES GRAPHS FROM RDS OBJECTS
################################################################################

#_____Load required packages____________________________________________________

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))
suppressMessages(library(pracma))
suppressMessages(library(ggpubr))
suppressMessages(library(harrypotter))
suppressMessages(library(nortest))
suppressMessages(library(hexbin))
suppressMessages(library(readxl))




#_____Read in arguments_________________________________________________________
print("Read in dataframe")

args = commandArgs(trailingOnly = TRUE)

input <- args[1]
int_name <- args[2]
pval_in <- args[3]
roi1_file <- args[4]
roi2_file <- args[5]
metadata_condensed <- args[6]
out <- args[7]


hm_zscore_df <- readRDS(file=input)



#______Set ggplot graph theme___________________________________________________

theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
            theme(plot.title = element_text(hjust = 0.5, size = 16),
                  panel.background = element_rect(fill = "white", colour = NA),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.spacing = unit(0.25, "lines"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                  #change the colour of facet label background
                  strip.background = element_rect(fill = "#E6E1EA"),
                  #remove space between facets
                  panel.spacing.x=unit(0, "lines"),
                  #legend.key = element_blank(),
                  legend.background=element_blank(),
                  #legend background
                  legend.key = element_rect(fill = NA),
                  legend.title = element_text(vjust = 0.8),
                  legend.position="top")
)




#___Heatmap along chromosome positions: ALL INTERACTIONS (incl. non-sig)________

#roi1
roi1 <- read.table(roi1_file)
colnames(roi1) <- c("chrom","start","end")
roi1$chrom <- as.factor(roi1$chrom)
print("ROI 1")
print(roi1)
#roi2
roi2 <- read.table(roi2_file)
colnames(roi2) <- c("chrom","start","end")
roi2$chrom <- as.factor(roi2$chrom)
print("ROI 2")
print(roi2)

#filter based on ROIs in both directions (i.e. chr12chr17 AND chr17chr12)
print("filter in both directions")
dirA <- paste0(unique(roi1$chrom),unique(roi2$chrom))
print(dirA)
dirB <- paste0(unique(roi2$chrom),unique(roi1$chrom))
print(dirB)
hm_df <- hm_zscore_df %>% filter(chrs == dirA | chrs == dirB)
head(hm_df)

#scale genomic positions by 1Mb
hm_df$st1 <- hm_df$st1/1000000
hm_df$end1 <- hm_df$end1/1000000
hm_df$st2 <- hm_df$st2/1000000
hm_df$end2 <- hm_df$end2/1000000
head(hm_df)



#___Generate dataframe containing chromosomal territories info__________________
#chrom info: centromere midpoint calculated from UCSC, approx.)

print("chromosome info")
chrInf <- data.frame( chrom = c("chr1","chr2",
                                "chr3","chr4",
                                "chr5","chr6",
                                "chr7","chrX",
                                "chr8","chr9",
                                "chr11","chr10",
                                "chr12","chr13",
                                "chr14","chr15",
                                "chr16","chr17",
                                "chr18","chr20",
                                "chr19","chrY",
                                "chr22","chr21"),
                      centromere = c(123252373.5,93787431.5,
                                     90856062,50074452.5,
                                     48585285.5,60557102.5,
                                     60058972.5,61016889,
                                     45249872,43893383.5,
                                     53454152,39800499.5,
                                     35764400,17692000.5,
                                     17117352,19037747.5,
                                     36878628.5,25067566.5,
                                     18464134,28099979.5,
                                     26161912,10470308,
                                     15520235.5,11917946),
                      chrClass = c("Metacentric","Metacentric",
                                   "Metacentric","Submetacentric",
                                   "Submetacentric","Submetacentric",
                                   "Submetacentric","Submetacentric",
                                   "Submetacentric","Submetacentric",
                                   "Submetacentric","Submetacentric",
                                   "Submetacentric","Acrocentric",
                                   "Acrocentric","Acrocentric",
                                   "Metacentric","Submetacentric",
                                   "Submetacentric","Metacentric",
                                   "Metacentric","Acrocentric",
                                   "Acrocentric","Acrocentric"),
                      size = c(248956422,242193529,
                               198295559,190214555,
                               181538259,170805979,
                               159345973,156040895,145138636,
                               138394717,135086622,
                               133797422,133275309,
                               114364328,107043718,
                               101991189,90338345,
                               83257441,80373285,
                               64444167,58617616,57227415,
                               50818468,46709983))

head(chrInf)



#____Preparing data for graphs (1)______________________________________________

# getting the valid inter positions to match with the x and y chroms in the graphs
print("get chromosomal positions for each ROI")

beds <- rbind(roi1,roi2)

CA <- unique(hm_df$chrA)
CA <- gsub("chr", "", CA)
cat("Chr A =",CA)
#Astart <- roi1$start /1000000
Astart <- beds$start[match(paste0("chr",CA),beds$chrom)] / 1000000
cat("ROI location:",Astart)

CB <- unique(hm_df$chrB)
CB <- gsub("chr", "", CB)
cat("Chr B =",CB)
#Bstart <- roi2$start /1000000
Bstart <- beds$start[match(paste0("chr",CB),beds$chrom)] / 1000000
cat("ROI location:",Bstart)

#reformat dataframe, count each inter twice
tp_A <- hm_df %>% select(chrA,st1,end1,cell,zscore,pvalue) 
colnames(tp_A) <- c("chr","st","end","cell","zscore","pvalue")
tp_B <- hm_df %>% select(chrB,st2,end2,cell,zscore,pvalue) 
colnames(tp_B) <- c("chr","st","end","cell","zscore","pvalue")
tp_dat <- rbind(tp_A,tp_B)

#get mean zscore
tp_dat_sum_all = tp_dat %>% group_by(chr,st,cell) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
tp_dat_sum_all$chr <- gsub("chr", "", tp_dat_sum_all$chr)

print("mean zscore per bin")
summary(tp_dat_sum_all)


# change cell name (dataset) with tissue ID
tissue_groups <- read.table(metadata_condensed, header=T)
for (t in 1:nrow(tissue_groups)){
tp_dat_sum_all$cell <- gsub(tissue_groups$Dataset[t], tissue_groups$Tissue_Identifier[t], tp_dat_sum_all$cell)
}



#____1: Tickplot/heatmap of mean zscore (all interactions)______________________________

print("#mean zscore tickplot for all interactions")
output=paste0(out,"/",int_name)
if (dir.exists(output) == FALSE) {dir.create(output) ; print("creating output directory")} else {print("directory exists")}

for(i in unique(tp_dat_sum_all$chr)) {
  interpos <- Astart
  if (CA == i) {
    interpos = Astart
  } else {
    interpos = Bstart
  }
  xmax <- chrInf$size[match(paste0("chr",i),chrInf$chrom)] / 1000000
  tpp <- tp_dat_sum_all %>% filter(chr == i)
  #plot
  print("plotting ...")
  cat("Chromosome",i,"\n")
  cat("Length of chromosome:",xmax,"\n")
  cat("Interaction location:",interpos,"\n")
  tp <- (ggplot(tpp, aes(x=st+0.5, y=1, fill = mzscore))
         + geom_tile(aes(fill = mzscore), width = 1, height = 1)
         + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score", na.value = "grey")
         + labs(x = paste0("Chromosome ", i, " position [Mb]"), title = "Trans-chromosomal interactions z-scores")
         + facet_grid(cell ~ .)
         + geom_vline(xintercept = interpos, colour = "red", size = 0.75)
         + expand_limits(x = c(0,xmax))
         + scale_y_continuous(expand = c(0, 0))
         + scale_x_continuous(expand = c(0, 0))
         + theme(strip.background = element_rect(fill = "white"),
                 strip.text.y.right = element_text(angle = 0,size=6),
                 panel.spacing = unit(0, "lines"),
                 #title
                 plot.title = element_text(size=14),
                 #legend
                 legend.title = element_text(size = 10, angle = 0),
                 legend.text=element_text(size = 12),
                 legend.position="top",
                 #x-axis
                 axis.title.x = element_text(size = 14, vjust = -1),
                 axis.text.x = element_text(size = 12),
                 #remove y-axis 
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank(),
                 
         )
  )
  filename <- paste0(output,"/StripePlot_mean_zscore_chrom",i,"_",int_name,".pdf")
  pdf(filename, width = 14, height = 8)
  print(tp)
  dev.off()
  print("done")
}



#____Preparing data for graphs (2)______________________________________________

#select only sig inters
tp_dat_sig <- tp_dat %>% filter(pvalue <= pval_in)
tp_dat_sig_pos <- tp_dat_sig %>% filter(zscore > 0)
tp_dat_sum_sig_pos = tp_dat_sig_pos %>% group_by(chr,st,cell) %>% dplyr::summarize(numSig=n())
tp_dat_sum_sig_pos$chr <- gsub("chr", "", tp_dat_sum_sig_pos$chr)
summary(tp_dat_sum_sig_pos)

# save sig df from bubble graph for quantification test
for(i in unique(tp_dat_sum_sig_pos$chr)) {
  interpos <- Astart
  if (CA == i) {
    interpos = Astart
  } else {
    interpos = Bstart
  }
  quant <- tp_dat_sum_sig_pos %>% filter(chr == i)
  write.table(quant, file = as.character(paste0(output,"/chrom",i,"_bubble_data.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

#df of mean zscore for bubble graph
bgz_dat <- tp_dat_sum_all %>% group_by(chr,st) %>% dplyr::summarize(mmzscore=mean(mzscore, na.rm = TRUE), .groups = "keep")

#df of mean number of interactions for bubble graph
bgi_dat <- tp_dat_sum_sig_pos %>% group_by(chr,st) %>% dplyr::summarize(mnumSig=mean(numSig, na.rm = TRUE), .groups = "keep")

#combine dfs
bg_dat <- merge(bgz_dat, bgi_dat, by = c('chr','st'))

print("mean zscore and number of significant interactions per bin")
head(bg_dat)

print("# Pearson correlation between mean zscore and mean number of sig inters per chrom per bin across cells")
cor.test(bg_dat$mmzscore, bg_dat$mnumSig, method = 'pearson')




#____2: Bubble plot of total number of significant interactions_________________

print("#total number of significant interactions bubble plot v1")

for(i in unique(bg_dat$chr)) {
  interpos <- Astart
  if (CA == i) {
    interpos = Astart
  } else {
    interpos = Bstart
  }
  xmax <- chrInf$size[match(paste0("chr",i),chrInf$chrom)] / 1000000
  bgg <- bg_dat %>% filter(chr == i)
  #plot
  print("plotting ...")
  cat("Chromosome",i,"\n")
  cat("Length of chromosome:",xmax,"\n")
  cat("Interaction location:",interpos,"\n")
  bgv1 <- (ggplot(bgg, aes(st, y= mmzscore, size=mnumSig, colour=mnumSig))
         + geom_point(alpha = 0.5)
         #adjusting size of bubbles to make them bigger
         + scale_size(range = c(0, 18))
         + scale_color_hp(discrete = FALSE, option = "ronweasley2", guide="legend")
         + labs(x = paste0("Chromosome ", i, " position [Mb]"),
                y = "Mean z-score per bin across cells",
                title = "Trans-chromosomal significant interactions (all cells)",
                size = "Mean number of \npositive significant \ninteractions per bin",
                color = "Mean number of \npositive significant \ninteractions per bin")
         + geom_vline(xintercept = interpos, colour = "red")
         + geom_smooth(aes(x = st, y = mmzscore), orientation = "x",span = 0.3,linewidth=0.75)
         + expand_limits(x = c(0,xmax))
         + scale_y_continuous(expand = c(0, 0))
         + scale_x_continuous(expand = c(0, 0))
         + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
                 strip.background = element_rect(fill = "white"),
                 panel.spacing = unit(0, "lines"),
                 legend.position="right",
                 legend.title = element_text(angle = 0, size = 8, vjust = 3))
  )
  filename <- paste0(output,"/BubblePlot_zscore_numSigInters_chrom",i,"_",int_name,".pdf")
  pdf(filename, width = 14, height = 8)
  print(bgv1)
  dev.off()
}




