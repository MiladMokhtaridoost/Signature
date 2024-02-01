#______Read in arguments________________________________________________________
options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

qvalue_pos <- read.table(args[1], header = T, sep = "\t")

#______Load required packages___________________________________________________
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)

#__________________________Set the theme__________________________________________________
theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
            theme(plot.title = element_text(hjust = 0.5, size = 13),
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
                  legend.position="right")
)

#________________________________________split column ID______________________________

colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
qvalue_pos$ID <- sub("B", "\\.B", as.character(qvalue_pos$ID))
qvalue_pos <- qvalue_pos %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
qvalue_pos$chrA <- sub("Achr", "\\chr", as.character(qvalue_pos$chrA))
qvalue_pos$chrB <- sub("Bchr", "\\chr", as.character(qvalue_pos$chrB))
head(qvalue_pos)
 
qvalue_pos <- qvalue_pos[,c(1,2,3,4,5,6,7)]
#_____________________________________chromosome mapped regions_______________________
#preparing chrInf table to find the mapping length of each chromosome 
#unmapped reads calculated from signature output
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

chrInf$st_centromere_unmapped <- c(120000000,92000000,
                          91000000,50000000,
                          47000000,59000000,
                          58000000,59000000,
                          44000000,43000000,
                          51000000,40000000,
                          35000000,16000000,
                         16000000,17000000,
                         36000000,23000000,
                         15000000,26000000,
                         25000000,10000000,
                         13000000,11000000)

chrInf$end_centromere_unmapped <- c(143000000,94000000,
                             94000000,52000000,
                             49000000,60000000,
                             61000000,62000000,
                             46000000,68000000,
                             54000000,41000000,
                             37000000,19000000,
                           19000000,25000000,
                          46000000,26000000,
                          21000000,28000000,
                         27000000,11000000,
                         16000000,13000000)
  
chrInf$st_unmapped_region1 <- c(144000000,87000000,
                                  NA,189000000,
                                  69000000,29000000,
                                  142000000,0,
                                  1000000,39000000,
                                  1000000,38000000,
                                  NA,0,
                                  0,0,
                                  15000000,0,
                                  NA,NA,
                                  54000000,0,
                                  0,0)
                                  
chrInf$end_unmapped_region1 <- c(145000000,88000000,
                                   NA,190000000,
                                   71000000,33000000,
                                   144000000,3000000,
                                   3000000,41000000,
                                   2000000,39000000,
                                   NA,19000000,
                                   19000000, 25000000,
                                  17000000,1000000,
                                 NA, NA,
                                   55000000,2000000,
                                    11000000,  9000000
                              
                                   )
  
chrInf$st_unmapped_region2 <- c(146000000,89000000,
                                  NA,NA,
                                  179000000,NA,
                                  NA,NA,
                                  7000000,NA,
                                  NA,NA,
                                  NA,NA,
                                  93000000, 28000000,
                                 21000000, 36000000,
                                 NA,NA,
                                  NA,23000000,
                                  18000000,NA
                                  )
  chrInf$end_unmapped_region2 <- c(147000000,91000000,
                                   NA,NA,
                                   180000000,NA,
                                   NA,NA,
                                   13000000,NA,
                                   NA, NA,
                                   NA,NA,
                                   94000000,33000000,
                                   23000000,39000000,
                                   NA,NA,
                                   NA, 26000000,
                                  19000000,NA
                                   
                                   )
  
chrInf$st_unmapped_region3 <- c(148000000,NA,
                                  NA,NA,
                                  NA,NA,
                                  NA,NA,
                                  39000000,NA,
                                  NA,NA,
                                  NA,NA, 
                                  106000000,NA,
                                  32000000,46000000,
                                  NA,NA,
                                  NA,27000000,
                                  NA, NA
                                 
                                  )
  chrInf$end_unmapped_region3 <- c(150000000,NA,
                                   NA,NA,
                                   NA,NA,
                                   NA, NA,
                                   40000000,NA,
                                   NA,NA,
                                   NA,NA,
                                   107000000,NA,
                                   34000000,47000000,
                                   NA,NA,
                                   NA,57000000,
                                   NA,NA
                                   )
  
  chrInf$st_unmapped_region4 <- c(243000000,NA,
                                  NA,NA,
                                  NA,NA, 
                                  NA,NA,
                                  NA,NA, 
                                  NA,NA,
                                  NA, NA,
                                  NA,NA,
                                  NA,NA,
                                  NA,NA,
                                  NA,NA,
                                  NA,NA
                                  
                                  )
  
  chrInf$end_unmapped_region4 <- c(244000000,NA,
                                   NA,NA,
                                   NA,NA,
                                   NA,NA,
                                   NA, NA,
                                   NA, NA,
                                   NA,NA,
                                   NA,NA,
                                   NA,NA,
                                   NA,NA,
                                   NA,NA,
                                   NA,NA)
  

 

  chrInf[is.na(chrInf)]<- 0
                            
   chrInf$unmapped_centromere_length <- (chrInf$end_centromere) - (chrInf$st_centromere) 
   
   chrInf$unmapped_region1_length <- (chrInf$end_unmapped_region1) - (chrInf$st_unmapped_region1)   
   chrInf$unmapped_region2_length <- (chrInf$end_unmapped_region2) - (chrInf$st_unmapped_region2)  
   chrInf$unmapped_region3_length <- (chrInf$end_unmapped_region3) - (chrInf$st_unmapped_region3)  
   chrInf$unmapped_region4_length <- (chrInf$end_unmapped_region4) - (chrInf$st_unmapped_region4)  
   chrInf$total_unmapped_length <- (chrInf$unmapped_centromere_length)+(chrInf$unmapped_region1_length)+(chrInf$unmapped_region2_length)+(chrInf$unmapped_region3_length)+(chrInf$unmapped_region4_length)
   chrInf$total_length_mapped <-  (chrInf$size)-(chrInf$total_unmapped_length )
   chrInf$five_percent  <- ((chrInf$total_length_mapped)*5)/100
   chrInf
   chrInf$two_point_five_percent <- chrInf$five_percent/2
   
   chrInf$q_arm <-  chrInf$size - chrInf$five_percent
   
#calculate q_arm 5% seperately for each chr
   #If any of the unmapped region length or their summation is gereater than 2.5 percent of mapped length of chr then put NA for that value 
   #chr1 
   chrInf$q_telo <- NA 
   chrInf$centromere_q_telo <- NA
   chrInf$p_telo <- NA
   chrInf$centromere_p_telo <- NA
     ifelse(chrInf$unmapped_region1_length[1] < chrInf$two_point_five_percent[1],TRUE,FALSE)
     ifelse(chrInf$unmapped_region2_length[1] < chrInf$two_point_five_percent[1],TRUE,FALSE)
     ifelse(chrInf$unmapped_region3_length[1] < chrInf$two_point_five_percent[1],TRUE,FALSE)
     ifelse(chrInf$unmapped_region4_length[1] < chrInf$two_point_five_percent[1],TRUE,FALSE)
     chrInf$summation_unmapped_regions_except_centromere<- NA
   chrInf$summation_unmapped_regions_except_centromere[1] <- chrInf$unmapped_region1_length[1] +chrInf$unmapped_region2_length[1] + chrInf$unmapped_region3_length[1] + chrInf$unmapped_region4_length[1] 
    ifelse(chrInf$summation_unmapped_regions_except_centromere[1] < chrInf$two_point_five_percent[1],TRUE,FALSE)
    
   (chrInf$st_unmapped_region1[1]-chrInf$end_centromere_unmapped[1] )+(chrInf$st_unmapped_region2[1]-chrInf$end_unmapped_region1[1])+(chrInf$st_unmapped_region3[1]-chrInf$end_unmapped_region2[1])+(chrInf$st_unmapped_region4[1]-chrInf$end_unmapped_region3[1])
    chrInf$five_percent[1]  -  ((chrInf$st_unmapped_region1[1]-chrInf$end_centromere_unmapped[1])+(chrInf$st_unmapped_region2[1]-chrInf$end_unmapped_region1[1])+(chrInf$st_unmapped_region3[1]-chrInf$end_unmapped_region2[1]))

   chrInf$centromere_q_telo[1]<- chrInf$end_unmapped_region3[1] + (chrInf$five_percent[1]  -  ((chrInf$st_unmapped_region1[1]-chrInf$end_centromere_unmapped[1])+(chrInf$st_unmapped_region2[1]-chrInf$end_unmapped_region1[1])+(chrInf$st_unmapped_region3[1]-chrInf$end_unmapped_region2[1])))
   chrInf$q_telo[1]<- ifelse(chrInf$unmapped_region4_length[1]< chrInf$two_point_five_percent[1],chrInf$st_unmapped_region4[1]  -(chrInf$five_percent[1]-(chrInf$size[1]-chrInf$end_unmapped_region4[1])),NA)
   
   chrInf$p_telo[1] <- chrInf$five_percent[1]
   chrInf$centromere_p_telo[1] <- chrInf$st_centromere_unmapped[1]-chrInf$five_percent[1]
   p_telo_chr1 <- c(0,chrInf$p_telo[1])
   centromere_p_telo_chr1 <- c( chrInf$centromere_p_telo[1],chrInf$st_centromere_unmapped[1])
   centromere_q_telo_chr1<- c(chrInf$end_centromere_unmapped[1],chrInf$centromere_q_telo[1])
   q_telo_chr1 <- c(chrInf$q_telo[1],chrInf$size[1])
   q_telo_chr1 
   
   #chr2 
   chrInf$p_telo[2] <- chrInf$five_percent[2]
   ifelse(chrInf$unmapped_region1_length[2] < chrInf$two_point_five_percent[2],TRUE,FALSE)
   ifelse(chrInf$unmapped_region2_length[2] < chrInf$two_point_five_percent[2],TRUE,FALSE)
   chrInf$summation_unmapped_regions_except_centromere[2]<- chrInf$unmapped_region1_length[2] +chrInf$unmapped_region2_length[2]
   chrInf$five_percent[2]-(chrInf$st_centromere_unmapped[2]-chrInf$end_unmapped_region2[2])
   chrInf$st_unmapped_region1[2]-(chrInf$five_percent[2]-(chrInf$st_centromere_unmapped[2]-chrInf$end_unmapped_region2[2]))
   chrInf$five_percent[2]-((chrInf$st_unmapped_region2[2] - chrInf$end_unmapped_region1[2])+(chrInf$st_centromere_unmapped[2]-chrInf$end_unmapped_region2[2]))
   chrInf$st_unmapped_region1[2]-(chrInf$five_percent[2]-((chrInf$st_unmapped_region2[2] - chrInf$end_unmapped_region1[2])+(chrInf$st_centromere_unmapped[2]-chrInf$end_unmapped_region2[2])))
   chrInf$centromere_p_telo[2]<- ifelse(chrInf$summation_unmapped_regions_except_centromere[2] < chrInf$two_point_five_percent[2],chrInf$st_unmapped_region1[2]-(chrInf$five_percent[2]-((chrInf$st_unmapped_region2[2] - chrInf$end_unmapped_region1[2])+(chrInf$st_centromere_unmapped[2]-chrInf$end_unmapped_region2[2]))) ,NA)
   chrInf$centromere_q_telo[2] <- chrInf$end_centromere_unmapped[2]+chrInf$five_percent[2]
   chrInf$q_telo[2]  <- chrInf$size[2] - chrInf$five_percent[2]
   
   p_telo_chr2 <- c(0,chrInf$p_telo[2])
   centromere_p_telo_chr2 <- c( chrInf$centromere_p_telo[2],chrInf$st_centromere_unmapped[2])
   centromere_q_telo_chr2<- c(chrInf$end_centromere_unmapped[2],chrInf$centromere_q_telo[2])
   q_telo_chr2 <- c(chrInf$q_telo[2],chrInf$size[2])
    
   
   #chr3
   chrInf$p_telo[3]<-  chrInf$five_percent[3]
   chrInf$centromere_p_telo[3]<- chrInf$st_centromere_unmapped[3]- chrInf$five_percent[3]
   chrInf$centromere_q_telo[3]<- chrInf$end_centromere_unmapped[3]+ chrInf$five_percent[3]
   chrInf$q_telo[3] <- chrInf$size[3]-chrInf$five_percent[3]
   p_telo_chr3 <- c(0,chrInf$p_telo[3])
   centromere_p_telo_chr3 <- c( chrInf$centromere_p_telo[3],chrInf$st_centromere_unmapped[3])
   centromere_q_telo_chr3<- c(chrInf$end_centromere_unmapped[3],chrInf$centromere_q_telo[3])
   q_telo_chr3 <- c(chrInf$q_telo[3],chrInf$size[3])
   
   
   
   #chr4 
   chrInf$summation_unmapped_regions_except_centromere[4] <- chrInf$end_unmapped_region1[4]- chrInf$st_unmapped_region1[4]
   chrInf$p_telo[4]<-  chrInf$five_percent[4]
   chrInf$centromere_p_telo[4]<- chrInf$st_centromere_unmapped[4]- chrInf$five_percent[4]
   chrInf$centromere_q_telo[4]<- chrInf$end_centromere_unmapped[4]+ chrInf$five_percent[4]
   chrInf$q_telo[4] <-   ifelse(chrInf$unmapped_region1_length[4]< chrInf$two_point_five_percent[4], chrInf$st_unmapped_region1[4]- chrInf$five_percent[4], NA)
   p_telo_chr4 <- c(0,chrInf$p_telo[4])
   centromere_p_telo_chr4 <- c( chrInf$centromere_p_telo[4],chrInf$st_centromere_unmapped[4])
   centromere_q_telo_chr4 <- c(chrInf$end_centromere_unmapped[4],chrInf$centromere_q_telo[4])
   q_telo_chr4 <- c(chrInf$q_telo[4],chrInf$st_unmapped_region1[4])    
   
                            
   #chr5 
   chrInf$summation_unmapped_regions_except_centromere[5] <- chrInf$unmapped_region1_length[5]+ chrInf$unmapped_region2_length[5]
   chrInf$p_telo[5]<-  chrInf$five_percent[5]
   chrInf$centromere_p_telo[5]<- chrInf$st_centromere_unmapped[5]- chrInf$five_percent[5]
   chrInf$centromere_q_telo[5]<- chrInf$end_centromere_unmapped[5]+ chrInf$five_percent[5]
   ifelse(chrInf$unmapped_region1_length[5] < chrInf$two_point_five_percent[5],TRUE,FALSE)
   ifelse(chrInf$unmapped_region2_length[5] < chrInf$two_point_five_percent[5],TRUE,FALSE)
   chrInf$q_telo[5]<-  chrInf$st_unmapped_region2[5]-(chrInf$five_percent[5]-(chrInf$size[5]-chrInf$end_unmapped_region2[5]))
   p_telo_chr5 <- c(0,chrInf$p_telo[5])
   centromere_p_telo_chr5 <- c( chrInf$centromere_p_telo[5],chrInf$st_centromere_unmapped[5])
   centromere_q_telo_chr5 <- c(chrInf$end_centromere_unmapped[5],chrInf$centromere_q_telo[5])
   q_telo_chr5 <- c(chrInf$q_telo[5],chrInf$size[5]) 
   
   
   
   #chr6
   chrInf$summation_unmapped_regions_except_centromere[6] <- chrInf$unmapped_region1_length[6]
   ifelse(chrInf$unmapped_region1_length[6] < chrInf$two_point_five_percent[6],TRUE,FALSE)
   chrInf$p_telo[6]<-  chrInf$five_percent[6]
   (chrInf$st_unmapped_region1[6]-(chrInf$five_percent[6]-(chrInf$st_centromere_unmapped[6]-chrInf$end_unmapped_region1[6]))) 
   chrInf$centromere_p_telo[6]<- ifelse(chrInf$st_centromere_unmapped[6]- chrInf$five_percent[6] > chrInf$end_unmapped_region1[6] ,chrInf$st_centromere_unmapped[6]- chrInf$five_percent[6] ,(chrInf$st_unmapped_region1[6]-(chrInf$five_percent[6]-(chrInf$st_centromere_unmapped[6]-chrInf$end_unmapped_region1[6]))) )
   chrInf$centromere_q_telo[6]<- chrInf$end_centromere_unmapped[6]+ chrInf$five_percent[6]
   chrInf$q_telo[6]<- chrInf$size[6]- chrInf$five_percent[6]
   p_telo_chr6 <- c(0,chrInf$p_telo[6])
   centromere_p_telo_chr6 <- c( chrInf$centromere_p_telo[6],chrInf$st_centromere_unmapped[6])
   centromere_q_telo_chr6 <- c(chrInf$end_centromere_unmapped[6],chrInf$centromere_q_telo[6])
   q_telo_chr6 <- c(chrInf$q_telo[6],chrInf$size[6]) 
   
    
#chr7   
   chrInf$summation_unmapped_regions_except_centromere[7] <- chrInf$unmapped_region1_length[7]
   chrInf$p_telo[7]<-  chrInf$five_percent[7]
   ifelse(chrInf$unmapped_region1_length[7] < chrInf$two_point_five_percent[7],TRUE,FALSE)
   chrInf$centromere_q_telo[7]<- chrInf$end_centromere_unmapped[7]+ chrInf$five_percent[7]
   chrInf$q_telo[7]<- chrInf$size[7]- chrInf$five_percent[7]
   chrInf$centromere_p_telo[7] <- chrInf$st_centromere_unmapped[7]- chrInf$five_percent[7]
   p_telo_chr7 <- c(0,chrInf$p_telo[7])
   centromere_p_telo_chr7 <- c( chrInf$centromere_p_telo[7],chrInf$st_centromere_unmapped[7])
   centromere_q_telo_chr7 <- c(chrInf$end_centromere_unmapped[7],chrInf$centromere_q_telo[7])
   q_telo_chr7 <- c(chrInf$q_telo[7],chrInf$size[7]) 
   
#chrX
   chrInf$summation_unmapped_regions_except_centromere[8] <- chrInf$unmapped_region1_length[8]
   ifelse(chrInf$unmapped_region1_length[8] < chrInf$two_point_five_percent[8],TRUE,FALSE)  
   chrInf$p_telo[8]<- ifelse(chrInf$unmapped_region1_length[8] < chrInf$two_point_five_percent[8],chrInf$end_unmapped_region1[8]+chrInf$five_percent[8],NA)  
   chrInf$centromere_p_telo[8]<- chrInf$st_centromere_unmapped[8] - chrInf$five_percent[8]
   chrInf$centromere_q_telo[8]<- chrInf$end_centromere_unmapped[8]+ chrInf$five_percent[8]
   chrInf$q_telo[8]<- chrInf$size[8]- chrInf$five_percent[8]
   p_telo_chrX <- c(chrInf$end_unmapped_region1[8],chrInf$p_telo[8])
   centromere_p_telo_chrX <- c( chrInf$centromere_p_telo[8],chrInf$st_centromere_unmapped[8])
   centromere_q_telo_chrX <- c(chrInf$end_centromere_unmapped[8],chrInf$centromere_q_telo[8])
   q_telo_chrX <- c(chrInf$q_telo[8],chrInf$size[8]) 

#chr8
   chrInf$summation_unmapped_regions_except_centromere[9] <- chrInf$unmapped_region1_length[9]+ chrInf$unmapped_region2_length[9]+ chrInf$unmapped_region3_length[9]
   ifelse(chrInf$unmapped_region1_length[9] < chrInf$two_point_five_percent[9],TRUE,FALSE)
   ifelse(chrInf$unmapped_region2_length[9] < chrInf$two_point_five_percent[9],TRUE,FALSE)
   ifelse(chrInf$unmapped_region3_length[9] < chrInf$two_point_five_percent[9],TRUE,FALSE)
   chrInf$st_unmapped_region1[9] - 0
   chrInf$st_unmapped_region2[9] -  chrInf$end_unmapped_region1[9]
   (chrInf$st_unmapped_region1[9] - 0)+ (chrInf$st_unmapped_region2[9] -  chrInf$end_unmapped_region1[9])
  chrInf$five_percent [9]- ((chrInf$st_unmapped_region1[9] - 0)+ (chrInf$st_unmapped_region2[9] -  chrInf$end_unmapped_region1[9]))
  # we do not have  chrInf$p_telo[9] because 5% will be after unmapped region 2 and the size of unmapped region2 is greater than 2.5% of mapped length of chromosome 8
  chrInf$p_telo[9]<- NA
  chrInf$five_percent[9]- (chrInf$st_centromere_unmapped[9]-chrInf$end_unmapped_region3[9])
  chrInf$st_unmapped_region3[9]-(chrInf$five_percent[9]- (chrInf$st_centromere_unmapped[9]-chrInf$end_unmapped_region3[9]))
  chrInf$centromere_p_telo[9]<- ifelse(chrInf$unmapped_region3_length[9]< chrInf$two_point_five_percent[9],chrInf$st_unmapped_region3[9]-(chrInf$five_percent[9]- (chrInf$st_centromere_unmapped[9]-chrInf$end_unmapped_region3[9])),NA)
  chrInf$centromere_q_telo[9]<- chrInf$end_centromere_unmapped[9]+ chrInf$five_percent[9]
  chrInf$q_telo[9]<- chrInf$size[9]- chrInf$five_percent[9]
  centromere_p_telo_chr8 <- c( chrInf$centromere_p_telo[9],chrInf$st_centromere_unmapped[9])
  centromere_q_telo_chr8 <- c(chrInf$end_centromere_unmapped[9],chrInf$centromere_q_telo[9])
  q_telo_chr8 <- c(chrInf$q_telo[9],chrInf$size[9]) 
  
  
  #chr9 
  ifelse(chrInf$unmapped_region1_length[10] < chrInf$two_point_five_percent[10],TRUE,FALSE)
  chrInf$summation_unmapped_regions_except_centromere[10] <- chrInf$unmapped_region1_length[10]
  chrInf$p_telo[10]<- chrInf$five_percent[10]
  chrInf$q_telo[10]<- chrInf$size[10] - chrInf$five_percent[10]
  chrInf$size[10]
  chrInf$centromere_q_telo[10]<- chrInf$end_centromere_unmapped[10]+ chrInf$five_percent[10]
  chrInf$five_percent[10] -(chrInf$st_centromere_unmapped[10] - chrInf$end_unmapped_region1[10])
  chrInf$st_unmapped_region1[10] - (chrInf$five_percent[10] -(chrInf$st_centromere_unmapped[10] - chrInf$end_unmapped_region1[10]))
  chrInf$centromere_p_telo[10]<- ifelse(chrInf$unmapped_region1_length[10] < chrInf$two_point_five_percent[10],chrInf$st_unmapped_region1[10] - (chrInf$five_percent[10] -(chrInf$st_centromere_unmapped[10] - chrInf$end_unmapped_region1[10])),NA)
  p_telo_chr9 <- c(0,chrInf$p_telo[10])
  centromere_p_telo_chr9 <- c( chrInf$centromere_p_telo[10],chrInf$st_centromere_unmapped[10])
  centromere_q_telo_chr9 <- c(chrInf$end_centromere_unmapped[10],chrInf$centromere_q_telo[10])
  q_telo_chr9 <- c(chrInf$q_telo[10],chrInf$size[10]) 
 
  #chr11
  chrInf$summation_unmapped_regions_except_centromere[11] <- chrInf$unmapped_region1_length[11]
  ifelse(chrInf$unmapped_region1_length[11] < chrInf$two_point_five_percent[11],TRUE,FALSE)
  chrInf$q_telo[11]<- chrInf$size[11] - chrInf$five_percent[11]
  chrInf$centromere_q_telo[11]<- chrInf$end_centromere_unmapped[11]+ chrInf$five_percent[11]
  chrInf$centromere_p_telo[11]<-chrInf$st_centromere_unmapped[11] - chrInf$five_percent[11]
  chrInf$five_percent[11]-(chrInf$st_unmapped_region1[11] - 0) 
  chrInf$end_unmapped_region1[11] +(chrInf$five_percent[11]-(chrInf$st_unmapped_region1[11] - 0)) 
  chrInf$p_telo[11]<- ifelse(chrInf$unmapped_region1_length[11] < chrInf$two_point_five_percent[11],chrInf$end_unmapped_region1[11] +(chrInf$five_percent[11]-(chrInf$st_unmapped_region1[11] - 0)),NA)
  p_telo_chr11 <- c(0,chrInf$p_telo[11])
  centromere_p_telo_chr11 <- c( chrInf$centromere_p_telo[11],chrInf$st_centromere_unmapped[11])
  centromere_q_telo_chr11 <- c(chrInf$end_centromere_unmapped[11],chrInf$centromere_q_telo[11])
  q_telo_chr11 <- c(chrInf$q_telo[11],chrInf$size[11]) 
  
  #chr10
  chrInf$summation_unmapped_regions_except_centromere[12] <- chrInf$unmapped_region1_length[12]
  ifelse(chrInf$unmapped_region1_length[12] < chrInf$two_point_five_percent[12],TRUE,FALSE)
  chrInf$q_telo[12]<- chrInf$size[12] - chrInf$five_percent[12]
  chrInf$centromere_q_telo[12]<- chrInf$end_centromere_unmapped[12]+ chrInf$five_percent[12]
  chrInf$p_telo[12]<- chrInf$five_percent[12]
  chrInf$five_percent[12] -( chrInf$st_centromere_unmapped[12]- chrInf$end_unmapped_region1[12])
  chrInf$st_unmapped_region1[12]- (chrInf$five_percent[12] -( chrInf$st_centromere_unmapped[12]- chrInf$end_unmapped_region1[12]))
  chrInf$centromere_p_telo[12] <-  ifelse(chrInf$unmapped_region1_length[12] < chrInf$two_point_five_percent[12],chrInf$st_unmapped_region1[12]- (chrInf$five_percent[12] -( chrInf$st_centromere_unmapped[12]- chrInf$end_unmapped_region1[12])),NA)
  p_telo_chr10 <- c(0,chrInf$p_telo[12])
  centromere_p_telo_chr10 <- c( chrInf$centromere_p_telo[12],chrInf$st_centromere_unmapped[12])
  centromere_q_telo_chr10 <- c(chrInf$end_centromere_unmapped[12],chrInf$centromere_q_telo[12])
  q_telo_chr10 <- c(chrInf$q_telo[12],chrInf$size[12]) 
  
  #chr12
  chrInf$summation_unmapped_regions_except_centromere[13] <- NA
  chrInf$q_telo[13]<- chrInf$size[13] - chrInf$five_percent[13]
  chrInf$centromere_q_telo[13]<- chrInf$end_centromere_unmapped[13]+ chrInf$five_percent[13]
  chrInf$p_telo[13]<- chrInf$five_percent[13]
  chrInf$centromere_p_telo[13]<- chrInf$st_centromere_unmapped[13]- chrInf$five_percent[13]
  p_telo_chr12 <- c(0,chrInf$p_telo[13])
  centromere_p_telo_chr12 <- c( chrInf$centromere_p_telo[13],chrInf$st_centromere_unmapped[13])
  centromere_q_telo_chr12 <- c(chrInf$end_centromere_unmapped[13],chrInf$centromere_q_telo[13])
  q_telo_chr12 <- c(chrInf$q_telo[13],chrInf$size[13]) 
  
  #chr13
  #acrocentric chromosome and does not have any p_arm length (all p_arm is unmapped)
  chrInf$summation_unmapped_regions_except_centromere[14] <- chrInf$unmapped_region1_length[14]
  ifelse(chrInf$unmapped_region1_length[14] < chrInf$two_point_five_percent[14],TRUE,FALSE)
  chrInf$p_telo[14]<- NA
  chrInf$centromere_p_telo[14]<- NA
  chrInf$centromere_q_telo[14]<- chrInf$end_centromere_unmapped[14]+ chrInf$five_percent[14]
  chrInf$q_telo[14]<- chrInf$size[14] - chrInf$five_percent[14]
  centromere_q_telo_chr13 <- c(chrInf$end_centromere_unmapped[14],chrInf$centromere_q_telo[14])
  q_telo_chr13 <- c(chrInf$q_telo[14],chrInf$size[14]) 
  
  #chr14
  #acrocentric chromosome and does not have any p_arm length (all p_arm is unmapped)
  chrInf$summation_unmapped_regions_except_centromere[15] <- chrInf$unmapped_region1_length[15]+ chrInf$unmapped_region2_length[15]+ chrInf$unmapped_region3_length[15]
  ifelse(chrInf$unmapped_region1_length[15] < chrInf$two_point_five_percent[15],TRUE,FALSE)
  ifelse(chrInf$unmapped_region2_length[15] < chrInf$two_point_five_percent[15],TRUE,FALSE)
  ifelse(chrInf$unmapped_region3_length[15] < chrInf$two_point_five_percent[15],TRUE,FALSE)
  chrInf$p_telo[15]<- NA
  chrInf$centromere_p_telo[15]<- NA
  chrInf$centromere_q_telo[15]<- chrInf$end_centromere_unmapped[15]+ chrInf$five_percent[15]
  chrInf$st_unmapped_region3[15] - chrInf$five_percent[15] 
  chrInf$q_telo[15]<-  chrInf$st_unmapped_region3[15] - chrInf$five_percent[15] 
  centromere_q_telo_chr14 <- c(chrInf$end_centromere_unmapped[15],chrInf$centromere_q_telo[15])
  q_telo_chr14 <- c(chrInf$q_telo[15], chrInf$st_unmapped_region3[15]) 
  
  
  #chr15
  #acrocentric chromosome and does not have any p_arm length (all p_arm is unmapped) 
  chrInf$summation_unmapped_regions_except_centromere[16] <- chrInf$unmapped_region1_length[16]+ chrInf$unmapped_region2_length[16]
  ifelse(chrInf$unmapped_region1_length[16] < chrInf$two_point_five_percent[16],TRUE,FALSE)
  ifelse(chrInf$unmapped_region2_length[16] < chrInf$two_point_five_percent[16],TRUE,FALSE)
# we do not have p_arm and also after centromere before 5% unmapped region 2 is greater than 2.5 % of the mapped length of the chromnosomes 
  chrInf$q_telo[16]<-  chrInf$size[16] -chrInf$five_percent[16]
  chrInf$p_telo[16]<- NA
  chrInf$centromere_p_telo[16]<- NA
  chrInf$centromere_q_telo[16]<- NA
  q_telo_chr15 <- c(chrInf$q_telo[16], chrInf$size[16])
  
  #chr16
  chrInf$summation_unmapped_regions_except_centromere[17] <- chrInf$unmapped_region1_length[17]+ chrInf$unmapped_region2_length[17]+ chrInf$unmapped_region3_length[17]
  ifelse(chrInf$unmapped_region1_length[17] < chrInf$two_point_five_percent[17],TRUE,FALSE)
  ifelse(chrInf$unmapped_region2_length[17] < chrInf$two_point_five_percent[17],TRUE,FALSE)
  ifelse(chrInf$unmapped_region3_length[17] < chrInf$two_point_five_percent[17],TRUE,FALSE)
  chrInf$q_telo[17]<- chrInf$size[17]- chrInf$five_percent[17]
  chrInf$centromere_q_telo[17]<- chrInf$end_centromere_unmapped[17]+chrInf$five_percent[17]
  chrInf$p_telo[17]<- chrInf$five_percent[17]
  chrInf$st_centromere_unmapped[17]  - chrInf$end_unmapped_region3[17]
  chrInf$five_percent[17]-(chrInf$st_centromere_unmapped[17]  - chrInf$end_unmapped_region3[17])
  chrInf$st_unmapped_region3[17] -(chrInf$five_percent[17]-(chrInf$st_centromere_unmapped[17]  - chrInf$end_unmapped_region3[17]))
  chrInf$centromere_p_telo[17]<- ifelse(chrInf$unmapped_region3_length[17] < chrInf$two_point_five_percent[17],chrInf$st_unmapped_region3[17] -(chrInf$five_percent[17]-(chrInf$st_centromere_unmapped[17]  - chrInf$end_unmapped_region3[17])),NA)
  p_telo_chr16 <- c(0,chrInf$p_telo[17])
  centromere_q_telo_chr16 <- c(chrInf$end_centromere_unmapped[17],chrInf$centromere_q_telo[17])
  q_telo_chr16 <- c(chrInf$q_telo[17],chrInf$size[17]) 
  
#chr17
chrInf$summation_unmapped_regions_except_centromere[18] <- chrInf$unmapped_region1_length[18]+ chrInf$unmapped_region2_length[18]+ chrInf$unmapped_region3_length[18]
ifelse(chrInf$unmapped_region1_length[18] < chrInf$two_point_five_percent[18],TRUE,FALSE)
ifelse(chrInf$unmapped_region2_length[18] < chrInf$two_point_five_percent[18],TRUE,FALSE)
ifelse(chrInf$unmapped_region3_length[18] < chrInf$two_point_five_percent[18],TRUE,FALSE)
chrInf$end_unmapped_region1[18]+chrInf$five_percent[18]
chrInf$p_telo[18]<- ifelse(chrInf$unmapped_region1_length[18] < chrInf$two_point_five_percent[18],chrInf$end_unmapped_region1[18]+chrInf$five_percent[18],NA)
chrInf$centromere_p_telo[18]<- chrInf$st_centromere_unmapped[18] - chrInf$five_percent[18]
chrInf$centromere_q_telo[18]<- chrInf$end_centromere_unmapped[18] + chrInf$five_percent[18]
chrInf$q_telo[18]<- chrInf$size[18]- chrInf$five_percent[18]
p_telo_chr17 <- c(chrInf$end_unmapped_region1[18],chrInf$p_telo[18])
centromere_p_telo_chr17 <- c( chrInf$centromere_p_telo[18],chrInf$st_centromere_unmapped[18])
centromere_q_telo_chr17 <- c(chrInf$end_centromere_unmapped[18],chrInf$centromere_q_telo[18])
q_telo_chr17 <- c(chrInf$q_telo[18],chrInf$size[18])

#chr18
chrInf$summation_unmapped_regions_except_centromere[19] <- NA
chrInf$p_telo[19]<-  chrInf$five_percent[19]
chrInf$centromere_q_telo[19]<- chrInf$end_centromere_unmapped[19]+ chrInf$five_percent[19]
chrInf$q_telo[19]<- chrInf$size[19]- chrInf$five_percent[19]
chrInf$centromere_p_telo[19] <- chrInf$st_centromere_unmapped[19]- chrInf$five_percent[19]
p_telo_chr18 <- c(0,chrInf$p_telo[19])
centromere_p_telo_chr18 <- c( chrInf$centromere_p_telo[19],chrInf$st_centromere_unmapped[19])
centromere_q_telo_chr18 <- c(chrInf$end_centromere_unmapped[19],chrInf$centromere_q_telo[19])
q_telo_chr18 <- c(chrInf$q_telo[19],chrInf$size[19]) 

#Chr20
chrInf$summation_unmapped_regions_except_centromere[20] <- NA
chrInf$p_telo[20]<-  chrInf$five_percent[20]
chrInf$centromere_q_telo[20]<- chrInf$end_centromere_unmapped[20]+ chrInf$five_percent[20]
chrInf$q_telo[20]<- chrInf$size[20]- chrInf$five_percent[20]
chrInf$centromere_p_telo[20] <- chrInf$st_centromere_unmapped[20]- chrInf$five_percent[20]
p_telo_chr20 <- c(0,chrInf$p_telo[20])
centromere_p_telo_chr20 <- c( chrInf$centromere_p_telo[20],chrInf$st_centromere_unmapped[20])
centromere_q_telo_chr20 <- c(chrInf$end_centromere_unmapped[20],chrInf$centromere_q_telo[20])
q_telo_chr20 <- c(chrInf$q_telo[20],chrInf$size[20]) 

#Chr19
chrInf$summation_unmapped_regions_except_centromere[21] <- chrInf$unmapped_region1_length[21]
ifelse(chrInf$unmapped_region1_length[21] < chrInf$two_point_five_percent[21],TRUE,FALSE)
chrInf$size[21]-chrInf$end_unmapped_region1[21]
ifelse(chrInf$size[21]- chrInf$five_percent[21] > chrInf$end_unmapped_region1[21],TRUE,FALSE)
chrInf$p_telo[21]<- chrInf$five_percent[21]
chrInf$centromere_p_telo[21]<- chrInf$st_centromere_unmapped[21] - chrInf$five_percent[21]
chrInf$centromere_q_telo[21]<- chrInf$end_centromere_unmapped[21] + chrInf$five_percent[21]
chrInf$q_telo[21]<- chrInf$size[21]- chrInf$five_percent[21]
p_telo_chr19 <- c(0,chrInf$p_telo[21])
centromere_p_telo_chr19 <- c( chrInf$centromere_p_telo[21],chrInf$st_centromere_unmapped[21])
centromere_q_telo_chr19 <- c(chrInf$end_centromere_unmapped[21],chrInf$centromere_q_telo[21])
q_telo_chr19 <- c(chrInf$q_telo[21],chrInf$size[21]) 

#chrY
chrInf$summation_unmapped_regions_except_centromere[22] <- chrInf$unmapped_region1_length[22]+chrInf$unmapped_region2_length[22]+chrInf$unmapped_region3_length[22]
ifelse(chrInf$unmapped_region1_length[22] < chrInf$two_point_five_percent[22],TRUE,FALSE)
ifelse(chrInf$unmapped_region2_length[22] < chrInf$two_point_five_percent[22],TRUE,FALSE)
ifelse(chrInf$unmapped_region3_length[22] < chrInf$two_point_five_percent[22],TRUE,FALSE)
chrInf$p_telo[22]<- NA
chrInf$centromere_p_telo[22] <- chrInf$st_centromere_unmapped[22]- chrInf$five_percent[22]
chrInf$centromere_q_telo[22] <- chrInf$end_centromere_unmapped[22]+ chrInf$five_percent[22]
chrInf$q_telo[22]<- NA
centromere_p_telo_chrY <- c( chrInf$centromere_p_telo[22],chrInf$st_centromere_unmapped[22])
centromere_q_telo_chrY <- c(chrInf$end_centromere_unmapped[22],chrInf$centromere_q_telo[22])


#chr22
chrInf$summation_unmapped_regions_except_centromere[23] <- chrInf$unmapped_region1_length[23]+chrInf$unmapped_region2_length[23]
ifelse(chrInf$unmapped_region1_length[23] < chrInf$two_point_five_percent[23],TRUE,FALSE)
ifelse(chrInf$unmapped_region2_length[23] < chrInf$two_point_five_percent[23],TRUE,FALSE)
chrInf$centromere_p_telo[23] <-ifelse(chrInf$st_centromere_unmapped[23]-chrInf$five_percent[23] > chrInf$end_unmapped_region1[23],chrInf$st_centromere_unmapped[23]-chrInf$five_percent[23],NA)
chrInf$centromere_q_telo[23] <-ifelse(chrInf$end_centromere_unmapped[23]+chrInf$five_percent[23] < chrInf$st_unmapped_region2[23], chrInf$end_centromere_unmapped[23]+chrInf$five_percent[23],NA)
chrInf$p_telo[23]<- NA
chrInf$q_telo[23]<- chrInf$size[23]- chrInf$five_percent[23]
centromere_p_telo_chr22 <- c( chrInf$centromere_p_telo[23],chrInf$st_centromere_unmapped[23])
centromere_q_telo_chr22 <- c(chrInf$end_centromere_unmapped[23],chrInf$centromere_q_telo[23])
q_telo_chr22 <- c(chrInf$q_telo[23],chrInf$size[23]) 


#chr21
chrInf$summation_unmapped_regions_except_centromere[24] <- chrInf$unmapped_region1_length[24]
ifelse(chrInf$unmapped_region1_length[24] < chrInf$two_point_five_percent[24],TRUE,FALSE)
chrInf$p_telo[24]<- NA
chrInf$centromere_p_telo[24] <-ifelse(chrInf$st_centromere_unmapped[24]-chrInf$five_percent[24] > chrInf$end_unmapped_region1[24],chrInf$st_centromere_unmapped[24]-chrInf$five_percent[24],NA)
chrInf$centromere_q_telo[24] <- chrInf$end_centromere_unmapped[24]+chrInf$five_percent[24]
chrInf$q_telo[24]<- chrInf$end_centromere_unmapped[24]+ chrInf$five_percent[24]
centromere_p_telo_chr21 <- c( chrInf$centromere_p_telo[24],chrInf$st_centromere_unmapped[24])
centromere_q_telo_chr21 <- c(chrInf$end_centromere_unmapped[24],chrInf$centromere_q_telo[24])
q_telo_chr21 <- c(chrInf$q_telo[24],chrInf$size[24]) 

#________________________________build a table that has regions for all chromosomes___________

chr_mapped_regions <- data.frame(
  chrom = rep(NA,24),
  p_telo_st = rep(NA,24),
  p_telo_end = rep(NA,24),
  centromere_p_arm_st = rep(NA,24),
  centromere_q_arm_end = rep(NA,24),
 q_telo_st = rep(NA,24),
  q_telo_end = rep(NA,24)
)
chr_mapped_regions$chrom = c("chr1","chr2",
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
          "chr22","chr21")

chr_mapped_regions$st_centromere_unmapped <- c(120000000,92000000,
                                   91000000,50000000,
                                   47000000,59000000,
                                   58000000,59000000,
                                   44000000,43000000,
                                   51000000,40000000,
                                   35000000,16000000,
                                   16000000,17000000,
                                   36000000,23000000,
                                   15000000,26000000,
                                   25000000,10000000,
                                   13000000,11000000)

chr_mapped_regions$end_centromere_unmapped <- c(143000000,94000000,
                                    94000000,52000000,
                                    49000000,60000000,
                                    61000000,62000000,
                                    46000000,68000000,
                                    54000000,41000000,
                                    37000000,19000000,
                                    19000000,25000000,
                                    46000000,26000000,
                                    21000000,28000000,
                                    27000000,11000000,
                                    16000000,13000000)
#p_telo_chr21
chr_mapped_regions$p_telo_st = c(0,0,
                                     0,0,
                                     0,0,
                                     0,3000000,
                                     NA,0, 
                                     0, 0,
                                     0, NA,
                                     NA, NA,
                                     0,1000000,
                                   0,0,
                                   0,NA,
                                   NA,NA)

chr_mapped_regions$p_telo_end = c(11047821,11859676,
                                      9764778,9360728,
                                      8826913,8290299,
                                      7717299,10502045,
                                      NA,5569736,
                                      7554331,6589871,
                                      6563765, NA,
                                     NA, NA,
                                     3716917,4762872,
                                     3718664,3122208,
                                     2780881,NA,
                                     NA,NA)


#centromere_p_telo_chr21

chr_mapped_regions$centromere_p_arm_st = c(108952179,77140324, 
                      81235222, 40639272,
                      38173087, 50709701,
                      50282701,51497955,
                      36293068,35430264,
                      44445669,32410129,
                      28436235,NA,
                      NA,NA,
                      NA,19237128,
                      11281336,22877792,
                      22219119,8938629,
                      11209077,9214501 )
                    
  
#centromere_q_telo_chr3
chr_mapped_regions$centromere_q_arm_end = c(158047821,105859676,
                                         103764778,61360728,
                                         57826913,68290299,
                                         68717299,69502045,
                                         52706932,73569736,
                                         60554331,47589871,
                                         43563765,23618216,
                                         23152186,NA,
                                         49716917,29762872,
                                         24718664,31122208,
                                         29780881,12061371,
                                         17790923,14785499)

 
# q_telo_chr3
chr_mapped_regions$q_telo_st =c(236908601,230333853,
                                188530781,179639272,
                                171711346,162515680,
                                151628674,
                                148538850,138431704,
                                132824981,128532291,
                                127207551,126711544,
                                109746112,101847814,
                                98791630,86621428,
                                79494569, 76654621,
                                61321959,55836735,
                                NA,49027545,
                                14785499 )

chr_mapped_regions$q_telo_end =c(248956422,242193529,
                                 198295559,189000000,
                                 181538259,170805979,
                                 159345973,
                                 156040895,145138636,
                                 138394717,135086622,
                                 133797422,133275309, 
                                 114364328,106000000,
                                 101991189, 90338345,
                                 83257441,80373285,
                                 64444167,58617616,
                                 NA,50818468,
                                 46709983 )
chr_mapped_regions$total_mapped_bins<- c(220,237,
                       195,187,
                       176,165,
                       154,150,
                       134,111,
                       131,131,
                       131,95,
                       86,71,
                       74,75,
                       74,62,
                       55,22,
                       35,35 )
#round up the p_telo_end and round down the q_telo_st for centromere_p_telo down round and for centromere_q_telo round up to expand the region 
chr_mapped_regions$p_telo_end<-round(chr_mapped_regions$p_telo_end ,-6)
chr_mapped_regions$centromere_q_arm_end <- round( chr_mapped_regions$centromere_q_arm_end,-6)
chr_mapped_regions$centromere_p_arm_st <- floor(chr_mapped_regions$centromere_p_arm_st/1000000)*1000000
chr_mapped_regions$q_telo_st <- floor(chr_mapped_regions$q_telo_st/1000000)*1000000

qvalue_pos[,c(3,4,6,7)] <- lapply(qvalue_pos[,c(3,4,6,7)],as.numeric)

 qvalue_pos$interaction_type_chrA <- "None" 
 for (i in 1:nrow(qvalue_pos)) {   
   chrA <- qvalue_pos[i, "chrA"] 
   
   if(chrA %in% c("chr8","chr13","chr14","chr15","chr16","chrY","chr22","chr21")){
     next
   }
   
   region_row <- chr_mapped_regions[chr_mapped_regions$chrom == chrA, ]   
   if ( region_row$p_telo_end >= qvalue_pos[i, "end1"] ||      
        region_row$q_telo_st <= qvalue_pos[i, "st1"]) {   
     qvalue_pos[i, "interaction_type_chrA"] <- "T"   } 
   else if
   (          
    qvalue_pos[i, "st1"] >= region_row$centromere_p_arm_st &&           
    region_row$centromere_q_arm_end >= qvalue_pos[i, "end1"]) {   
     qvalue_pos[i, "interaction_type_chrA"] <- "C"   } }
 


#chr8,chr21,chr22
 for (i in 1:nrow(qvalue_pos)) {   
   chrA <- qvalue_pos[i, "chrA"] 
   
   if(!(chrA %in% c("chr8","chr21","chr22"))){
     next
   }
   
   region_row <- chr_mapped_regions[chr_mapped_regions$chrom == chrA, ]   
   if ( region_row$q_telo_st <= qvalue_pos[i, "st1"]) {   
     qvalue_pos[i, "interaction_type_chrA"] <- "T"   } 
   else if
   (          
     qvalue_pos[i, "st1"] >= region_row$centromere_p_arm_st &&           
     region_row$centromere_q_arm_end >= qvalue_pos[i, "end1"]) {   
     qvalue_pos[i, "interaction_type_chrA"] <- "C"   } }
 
 
#chr13,chr14
 for (i in 1:nrow(qvalue_pos)) {   
   chrA <- qvalue_pos[i, "chrA"] 
   
   if(!(chrA %in% c("chr13","chr14"))){
     next
   }
   
   region_row <- chr_mapped_regions[chr_mapped_regions$chrom == chrA, ]   
   if ( region_row$q_telo_st <= qvalue_pos[i, "st1"]) {   
     qvalue_pos[i, "interaction_type_chrA"] <- "T"   } 
   else if
   (  qvalue_pos[i, "st1"] >= region_row$st_centromere_unmapped &&           
     region_row$centromere_q_arm_end >= qvalue_pos[i, "end1"]) {   
     qvalue_pos[i, "interaction_type_chrA"] <- "C"   } }
 
#chr15 
 for (i in 1:nrow(qvalue_pos)) {   
   chrA <- qvalue_pos[i, "chrA"] 
   
   if(!(chrA %in% c("chr15"))){
     next
   }
   
   region_row <- chr_mapped_regions[chr_mapped_regions$chrom == chrA, ]   
   if ( region_row$q_telo_st <= qvalue_pos[i, "st1"]) {   
     qvalue_pos[i, "interaction_type_chrA"] <- "T"   } 
}

 #chr16
 for (i in 1:nrow(qvalue_pos)) {   
   chrA <- qvalue_pos[i, "chrA"] 
   
   if(!(chrA %in% c("chr16"))){
     next
   }
   
   region_row <- chr_mapped_regions[chr_mapped_regions$chrom == chrA, ]   
   if ( region_row$q_telo_st <= qvalue_pos[i, "st1"]) {   
     qvalue_pos[i, "interaction_type_chrA"] <- "T"   } 
   else if
   (          
     qvalue_pos[i, "st1"] >= region_row$st_centromere_unmapped  &&           
     region_row$centromere_q_arm_end >= qvalue_pos[i, "end1"]) {   
     qvalue_pos[i, "interaction_type_chrA"] <- "C"   } }
 
 
 #chrY

 for (i in 1:nrow(qvalue_pos)) {   
   chrA <- qvalue_pos[i, "chrA"] 
   
   if(!(chrA %in% c("chrY"))){
     next
   }
   
   region_row <- chr_mapped_regions[chr_mapped_regions$chrom == chrA, ]   
   if (  qvalue_pos[i, "st1"] >= region_row$st_centromere_unmapped  &&           
         region_row$centromere_q_arm_end >= qvalue_pos[i, "end1"])  {   
     qvalue_pos[i, "interaction_type_chrA"] <- "C"   } 
 }
 
 
#column chrB
 qvalue_pos$interaction_type_chrB <- "None" 
 for (i in 1:nrow(qvalue_pos)) {   
   chrB <- qvalue_pos[i, "chrB"] 
   
   if(chrB %in% c("chr8","chr13","chr14","chr15","chr16","chrY","chr22","chr21")){
     next
   }
   
   region_row <- chr_mapped_regions[chr_mapped_regions$chrom == chrB, ]   
   if ( region_row$p_telo_end >= qvalue_pos[i, "end2"] ||      
        region_row$q_telo_st <= qvalue_pos[i, "st2"]) {   
     qvalue_pos[i, "interaction_type_chrB"] <- "T"   } 
   else if
   (          
     qvalue_pos[i, "st2"] >= region_row$centromere_p_arm_st &&           
     region_row$centromere_q_arm_end >= qvalue_pos[i, "end2"]) {   
     qvalue_pos[i, "interaction_type_chrB"] <- "C"   } }
 
 
 
 #columnB#chr8,chr21,chr22
 for (i in 1:nrow(qvalue_pos)) {   
   chrB <- qvalue_pos[i, "chrB"] 
   
   if(!(chrB %in% c("chr8","chr21","chr22"))){
     next
   }
   
   region_row <- chr_mapped_regions[chr_mapped_regions$chrom == chrB, ]   
   if ( region_row$q_telo_st <= qvalue_pos[i, "st2"]) {   
     qvalue_pos[i, "interaction_type_chrA"] <- "T"   } 
   else if
   (          
     qvalue_pos[i, "st2"] >= region_row$centromere_p_arm_st &&           
     region_row$centromere_q_arm_end >= qvalue_pos[i, "end2"]) {   
     qvalue_pos[i, "interaction_type_chrA"] <- "C"   } }
 
 
#columnB #chr13,chr14
 for (i in 1:nrow(qvalue_pos)) {   
   chrB <- qvalue_pos[i, "chrB"] 
   
   if(!(chrB %in% c("chr13","chr14"))){
     next
   }
   
   region_row <- chr_mapped_regions[chr_mapped_regions$chrom == chrB, ]   
   if ( region_row$q_telo_st <= qvalue_pos[i, "st2"]) {   
     qvalue_pos[i, "interaction_type_chrB"] <- "T"   } 
   else if
   (  qvalue_pos[i, "st2"] >= region_row$st_centromere_unmapped &&           
      region_row$centromere_q_arm_end >= qvalue_pos[i, "end2"]) {   
     qvalue_pos[i, "interaction_type_chrB"] <- "C"   } }
 
 #columnB#chr15 
 for (i in 1:nrow(qvalue_pos)) {   
   chrB <- qvalue_pos[i, "chrB"] 
   
   if(!(chrB %in% c("chr15"))){
     next
   }
   
   region_row <- chr_mapped_regions[chr_mapped_regions$chrom == chrB, ]   
   if ( region_row$q_telo_st <= qvalue_pos[i, "st2"]) {   
     qvalue_pos[i, "interaction_type_chrB"] <- "T"   } 
 }
 
 #ColumnB#chr16
 for (i in 1:nrow(qvalue_pos)) {   
   chrB <- qvalue_pos[i, "chrB"] 
   
   if(!(chrB %in% c("chr16"))){
     next
   }
   
   region_row <- chr_mapped_regions[chr_mapped_regions$chrom == chrB, ]   
   if ( region_row$q_telo_st <= qvalue_pos[i, "st2"]) {   
     qvalue_pos[i, "interaction_type_chrB"] <- "T"   } 
   else if
   (          
     qvalue_pos[i, "st2"] >= region_row$st_centromere_unmapped  &&           
     region_row$centromere_q_arm_end >= qvalue_pos[i, "end2"]) {   
     qvalue_pos[i, "interaction_type_chrB"] <- "C"   } }
 
 
#ColumnB #chrY
 
 for (i in 1:nrow(qvalue_pos)) {   
   chrB <- qvalue_pos[i, "chrB"] 
   
   if(!(chrB %in% c("chrY"))){
     next
   }
   
   region_row <- chr_mapped_regions[chr_mapped_regions$chrom == chrB, ]   
   if (  qvalue_pos[i, "st2"] >= region_row$st_centromere_unmapped  &&           
         region_row$centromere_q_arm_end >= qvalue_pos[i, "end2"])  {   
     qvalue_pos[i, "interaction_type_chrB"] <- "C"   } 
 }
 
 qvalue_pos$total_interaction_type <- NA
 for (i in 1:nrow(qvalue_pos)) {

     if (qvalue_pos$interaction_type_chrA[i] == "T" && qvalue_pos$interaction_type_chrB[i] == "T") {
       qvalue_pos$total_interaction_type[i] <- "TT"
     } else if (qvalue_pos$interaction_type_chrA[i] == "C" && qvalue_pos$interaction_type_chrB[i] == "C") {
       qvalue_pos$total_interaction_type[i] <- "CC"
     } else if (qvalue_pos$interaction_type_chrA[i] == "C" && qvalue_pos$interaction_type_chrB[i] == "T") {
       qvalue_pos$total_interaction_type[i] <- "CT"
     } else if (qvalue_pos$interaction_type_chrA[i] == "T" && qvalue_pos$interaction_type_chrB[i] == "C") {
       qvalue_pos$total_interaction_type[i] <- "CT"
     } else if (qvalue_pos$interaction_type_chrA[i] == "C" && qvalue_pos$interaction_type_chrB[i] == "None") {
       qvalue_pos$total_interaction_type[i] <- "CN"
     } else if (qvalue_pos$interaction_type_chrA[i] == "T" && qvalue_pos$interaction_type_chrB[i] == "None") {
       qvalue_pos$total_interaction_type[i] <- "TN"
     } else if (qvalue_pos$interaction_type_chrA[i] == "None" && qvalue_pos$interaction_type_chrB[i] == "T") {
       qvalue_pos$total_interaction_type[i] <- "TN"
     } else if (qvalue_pos$interaction_type_chrA[i] == "None" && qvalue_pos$interaction_type_chrB[i] == "C") {
       qvalue_pos$total_interaction_type[i] <- "CN"
     } else {
       qvalue_pos$total_interaction_type[i] <- "NN"
     }
 }
 
 TT_count <- sum(qvalue_pos$total_interaction_type == "TT")
 CT_count <- sum(qvalue_pos$total_interaction_type == "CT")
 CC_count <- sum(qvalue_pos$total_interaction_type == "CC")
 NN_count <- sum(qvalue_pos$total_interaction_type == "NN")
 TN_count <- sum(qvalue_pos$total_interaction_type == "TN")
 CN_count <- sum(qvalue_pos$total_interaction_type == "CN")
 
 
 #expected number of each interaction type 
 #for chromosomes chr1,chr2,chr3,chr4,chr5,chr6,chr7,chrX,chr9,chr11,chr10,chr12,chr16,chr17,chr18,chr20,chr19 because they have both telomeric region
 
 T1 <- sum(chr_mapped_regions$total_mapped_bins[c(1:8,10:13,17:21)] * 0.1)
 
 #chr8,chr13,chr14,chr15,chr22,chr21
 
 T2 <- sum(chr_mapped_regions$total_mapped_bins[c(9,14:16,23:24)] * 0.05)
 #chrY does not have telomeric region 
 T = T1+T2
 
# chr1 to chr12 and chr17 to chr21 (all except chr13,chr14,chr15,chr16)
C1 <- sum(chr_mapped_regions$total_mapped_bins[c(1:13,18:24)] * 0.1)
C2 <- sum(chr_mapped_regions$total_mapped_bins[c(14:15,17)] * 0.05) 
#chr15 does not have centromeric region 
C = C1+C2

N1 <- sum(chr_mapped_regions$total_mapped_bins[c(1:8,10:13,18:21)] * 0.8)
N2 <- sum(chr_mapped_regions$total_mapped_bins[c(9,23:24)] * 0.85)
N3 <- sum(chr_mapped_regions$total_mapped_bins[14:15]*0.9)
N4 <- sum(chr_mapped_regions$total_mapped_bins[16]*0.95)
N4 <- sum(chr_mapped_regions$total_mapped_bins[17]*0.85)
N5 <- sum(chr_mapped_regions$total_mapped_bins[22]*0.9)
N=N1+N2+N3+N4+N5

#Expected_values
probability_T<- (T/(T+C+N))
probability_N<- (N/(T+C+N))
probability_C<- (C/(T+C+N))

expected_C_C <- (probability_C)*(probability_C)*40282
expected_T_T <- (probability_T)*(probability_T)*40282
expected_N_N <- (probability_N)*(probability_N)*40282
expected_C_N <- (probability_C)*(probability_N)*40282*2
expected_T_N <- (probability_T)*(probability_N)*40282*2
expected_C_T <- (probability_C)*(probability_T)*40282*2

expected_value <- c(expected_C_C,expected_T_T,expected_N_N,expected_C_N,expected_T_N,expected_C_T)
observed_value <- c(CC_count,TT_count,NN_count,CN_count,TN_count,CT_count)
df <- data.frame (expexted = expected_value, observed = observed_value, row.names =  c("CC","TT","NN","CN","TN","CT") )

data<- data.frame(total_interaction = c(df$expexted,df$observed),
                  type = rep(c("expected","observed"), each=6),
                  interaction_type = rep(c("CC","TT","NN","CN","TN","CT"),times=2)
                  ) 
expected_value <- c(expected_C_C,expected_T_T,expected_C_T)
observed_value <- c(CC_count,TT_count,CT_count)
df <- data.frame (expexted = expected_value, observed = observed_value, row.names =  c("CC","TT","CT") )

data<- data.frame(total_interaction = c(df$expexted,df$observed),
                  type = rep(c("expected","observed"), each=3),
                  interaction_type = rep(c("CC","TT","CT"),times=1)
)
 
 ggplot(data, aes(x = interaction_type, y = total_interaction, fill =type)) +
   geom_bar(width=0.7, position=position_dodge(width=0.75), stat="identity") +
   labs(title = " Expected and observed number of interactions in mapped region",  x = "interaction type", y = "Value") +   
   scale_fill_manual(values = c("red", "blue"),labels=c("expected","observed")) +
   scale_y_continuous(labels = function(x) paste0(round((x/40282)*100),"%"))+theme(plot.title = element_text(size=8))
 
ggsave("plot_without_N.pdf")

#testing 
#probability
probability_C_T <- (probability_C)*(probability_T)*2
probability_C_C <- (probability_C)*(probability_C)
probability_T_T <- (probability_T)*(probability_T)
probability_N_N <- (probability_N)*(probability_N)
probability_T_N <- (probability_N)*(probability_T)*2
probability_C_N <- (probability_C)*(probability_N)*2

Result_C_C<-binom.test(1595,40282,probability_C_C)
Result_C_C

Result_T_T<-binom.test(4118,40282,probability_T_T)
Result_T_T

Result_N_N<-binom.test(16213,40282,probability_N_N)
Result_N_N

Result_T_N<-binom.test(12669,40282,probability_T_N)
Result_T_N

Result_C_N<-binom.test(4907,40282,probability_C_N)
Result_C_N

Result_C_T<-binom.test(780,40282,probability_C_T)
Result_C_T


 df <- data.frame(   Group = c("CT", "CC-TT","CC-TT"),   Subgroup = c("CT", "CC", "TT"),   Value = c(780, 1595, 4118) ) 
 bar_plot <- ggplot(df, aes(x = Group, y = Value, fill = Subgroup)) +   geom_bar(stat = "identity") +   labs(title = "Bar Plot", x = "Group", y = "Value") 
 stacked_bar_plot <- ggplot(df, aes(x = Group, y = Value, fill = Subgroup)) +   geom_bar(stat = "identity") +   labs(title = "Stacked Bar Plot", x = "Group", y = "Value") +   scale_fill_manual(values = c("CT" = "blue", "CC" = "green", "TT" = "red")) 
 print(bar_plot)
 ggsave("plot_observed_two_bars.pdf")
