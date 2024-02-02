################################################################################
# NOTE: THIS SCRIPT IS THE PRE-PROCESSING STEP AND ONLY SAVES AN RData OBJECT
################################################################################
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
options(scipen = 999)


#______Read in arguments________________________________________________________
# running with a scheduler #
options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

qdat_file <- args[1]
qalldat_file <- args[2]
f_in <- args[3]
m_in <- args[4]
output <- args[5]
pqdat <- args[6]


qdat <- read.table(qdat_file, header = T)
qalldat <- read.table(qalldat_file, header = T)

f_list <- read.table(f_in)
f_list <- c("ID",f_list$V1)
m_list <- read.table(m_in)
m_list <- c("ID",m_list$V1)

qdat_f <- qdat %>% select(contains(f_list))
qdat_m <- qdat %>% select(contains(m_list))

qalldat_f <- qalldat %>% select(contains(f_list))
qalldat_m <- qalldat %>% select(contains(m_list))



#_____Get percent of cells sig F________________________________________________

# prepare null dataframe (how many cells is the bin significant in (q<0.05), in percentage)
list2_f <- c()
list2_f$ID <- qdat_f$ID
list2_f <- as.data.frame(list2_f)
list2_f$count <- NA
list2_f$NAs <- NA
list2_f$available <- NA
list2_f$percent <- NA


total <- ncol(qdat_f)-1

for (i in 1:nrow(qdat_f)){
  row <- which(qalldat_f$ID == qdat_f$ID[i])

  list2_f$count[i] <- sum(!is.na(qdat_f[i,2:ncol(qdat_f)]))
  list2_f$NAs[i] <- sum(is.na(qalldat_f[row,]))
  list2_f$available[i] <- total - list2_f$NAs[i]
  list2_f$percent[i] <- ( list2_f$count[i] / list2_f$available[i] ) * 100
  cat(sprintf("index=%s",i), sep = "\n")

}

list2_f <- list2_f %>% select(c("ID","count","percent"))

head(list2_f)

saveRDS(list2_f, paste0(output,"/list2_f.rds"))



#_____Get percent of cells sig M________________________________________________

# prepare null dataframe (how many cells is the bin significant in (q<0.05), in percentage)
list2_m <- c()
list2_m$ID <- qdat_m$ID
list2_m <- as.data.frame(list2_m)
list2_m$count <- NA
list2_m$NAs <- NA
list2_m$available <- NA
list2_m$percent <- NA


total <- ncol(qdat_m)-1

for (i in 1:nrow(qdat_m)){
  row <- which(qalldat_m$ID == qdat_m$ID[i])

  list2_m$count[i] <- sum(!is.na(qdat_m[i,2:ncol(qdat_m)]))
  list2_m$NAs[i] <- sum(is.na(qalldat_m[row,]))
  list2_m$available[i] <- total - list2_m$NAs[i]
  list2_m$percent[i] <- ( list2_m$count[i] / list2_m$available[i] ) * 100
  cat(sprintf("index=%s",i), sep = "\n")

}

list2_m <- list2_m %>% select(c("ID","count","percent"))

head(list2_m)

saveRDS(list2_m, paste0(output,"/list2_m.rds"))



#_____split up data frames______________________________________________________

# remove zeros 
list2_f <- list2_f[c(which(list2_f$count >= 1)), ]
list2_m <- list2_m[c(which(list2_m$count >= 1)), ]

print("# subsetting data into 4 dataframes")

# data for plot A ( chrX vs ...) #

# filter for chromosome x interactions
grep <- grep("chrX", list2_f$ID)
data_A_female <- list2_f[c(grep),]
grep <- grep("chrX", list2_m$ID)
data_A_male <- list2_m[c(grep),]
grep <- grep("chrY", data_A_male$ID)
data_A_male <- data_A_male[-c(grep),]


# data for plot B/C ( chr Y vs ...) #

# filter for chromosome Y interactions
grep <- grep("chrY", list2_m$ID)
data_B <- list2_m[c(grep),]
# remove any chrX interaction and save it as data_C
grep <- grep("chrX", data_B$ID)
data_C <- data_B[c(grep),]
data_B <- data_B[-c(grep),]

print("A) chrX (f) vs autosomes ...")
head(data_A_female$ID)
print("A) chrX (m) vs autosomes ...")
head(data_A_male$ID)
print("B) chrY vs autosomes ...")
head(data_B$ID)
print("C) chrX vs chrY ...")
head(data_C$ID)



#_____Separate ID and sort into p and q arms____________________________________
# split up ID column #
print("# separate ID")
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")

# A - female
data_A_female$ID <- sub("B", "\\.B", as.character(data_A_female$ID))
data_A_female <- data_A_female %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
data_A_female$chrA <- gsub("A", "", data_A_female$chrA)
data_A_female$chrB <- gsub("B", "", data_A_female$chrB)
# A - male
data_A_male$ID <- sub("B", "\\.B", as.character(data_A_male$ID))
data_A_male <- data_A_male %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
data_A_male$chrA <- gsub("A", "", data_A_male$chrA)
data_A_male$chrB <- gsub("B", "", data_A_male$chrB)

# B
data_B$ID <- sub("B", "\\.B", as.character(data_B$ID))
data_B <- data_B %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
data_B$chrA <- gsub("A", "", data_B$chrA)
data_B$chrB <- gsub("B", "", data_B$chrB)

# C
data_C$ID <- sub("B", "\\.B", as.character(data_C$ID))
data_C <- data_C %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
data_C$chrA <- gsub("A", "", data_C$chrA)
data_C$chrB <- gsub("B", "", data_C$chrB)



# annotate with p and q arm #
print("# annotate with p and q arm")
pqarm <- fread(pqdat, select = c("id","chr1","st1","end1","chr2","st2","end2","target","interaction_type_target","anchor","interaction_type_anchor"))
colnames(pqarm)[2] <- "chrA"
colnames(pqarm)[5] <- "chrB"

# A - female
data_A_female$target <- NA
data_A_female$interaction_type_target <- NA
data_A_female$anchor <- NA
data_A_female$interaction_type_anchor <- NA

for (i in 1:nrow(data_A_female)){
index <- which(pqarm$chrA == data_A_female$chrA[i] &
               pqarm$st1 == data_A_female$st1[i] &
               pqarm$end1 == data_A_female$end1[i] &
               pqarm$chrB == data_A_female$chrB[i] &
               pqarm$st2 == data_A_female$st2[i] &
               pqarm$end2 == data_A_female$end2[i])

data_A_female$target[i] = pqarm$target[index[1]] 
data_A_female$interaction_type_target[i] = pqarm$interaction_type_target[index[1]] 
data_A_female$anchor[i] = pqarm$anchor[index[1]] 
data_A_female$interaction_type_anchor[i] = pqarm$interaction_type_anchor[index[1]] 
cat(paste0("Index ",i),"\n")
}
# A - male
data_A_male$target <- NA
data_A_male$interaction_type_target <- NA
data_A_male$anchor <- NA
data_A_male$interaction_type_anchor <- NA

for (i in 1:nrow(data_A_male)){
  index <- which(pqarm$chrA == data_A_male$chrA[i] &
                   pqarm$st1 == data_A_male$st1[i] &
                   pqarm$end1 == data_A_male$end1[i] &
                   pqarm$chrB == data_A_male$chrB[i] &
                   pqarm$st2 == data_A_male$st2[i] &
                   pqarm$end2 == data_A_male$end2[i])
  
  data_A_male$target[i] = pqarm$target[index[1]] 
  data_A_male$interaction_type_target[i] = pqarm$interaction_type_target[index[1]] 
  data_A_male$anchor[i] = pqarm$anchor[index[1]] 
  data_A_male$interaction_type_anchor[i] = pqarm$interaction_type_anchor[index[1]] 
  cat(paste0("Index ",i),"\n")
}

# B
data_B$target <- NA
data_B$interaction_type_target <- NA
data_B$anchor <- NA
data_B$interaction_type_anchor <- NA

for (i in 1:nrow(data_B)){
  index <- which(pqarm$chrA == data_B$chrA[i] &
                   pqarm$st1 == data_B$st1[i] &
                   pqarm$end1 == data_B$end1[i] &
                   pqarm$chrB == data_B$chrB[i] &
                   pqarm$st2 == data_B$st2[i] &
                   pqarm$end2 == data_B$end2[i])
  
  data_B$target[i] = pqarm$target[index[1]] 
  data_B$interaction_type_target[i] = pqarm$interaction_type_target[index[1]] 
  data_B$anchor[i] = pqarm$anchor[index[1]] 
  data_B$interaction_type_anchor[i] = pqarm$interaction_type_anchor[index[1]] 
  cat(paste0("Index ",i),"\n")
}


# C
data_C$target <- NA
data_C$interaction_type_target <- NA
data_C$anchor <- NA
data_C$interaction_type_anchor <- NA

for (i in 1:nrow(data_C)){
  index <- which(pqarm$chrA == data_C$chrA[i] &
                   pqarm$st1 == data_C$st1[i] &
                   pqarm$end1 == data_C$end1[i] &
                   pqarm$chrB == data_C$chrB[i] &
                   pqarm$st2 == data_C$st2[i] &
                   pqarm$end2 == data_C$end2[i])
  
  data_C$target[i] = pqarm$target[index[1]] 
  data_C$interaction_type_target[i] = pqarm$interaction_type_target[index[1]] 
  data_C$anchor[i] = pqarm$anchor[index[1]] 
  data_C$interaction_type_anchor[i] = pqarm$interaction_type_anchor[index[1]] 
  cat(paste0("Index ",i),"\n")
}


print("A) chrX (f) vs autosomes ...")
head(data_A_female[,c(1:7,10:13)])
print("A) chrX (m) vs autosomes ...")
head(data_A_male[,c(1:7,10:13)])
print("B) chrY vs autosomes ...")
head(data_B[,c(1:7,10:13)])
print("C) chrX vs chrY ...")
head(data_C[,c(1:7,10:13)])


#_____linear conversion to uniform scale________________________________________

chrInf <- data.frame( chrom = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chrX",
                                "chr8","chr9","chr11","chr10","chr12","chr13","chr14",
                                "chr15","chr16","chr17","chr18","chr20","chr19","chrY",
                                "chr22","chr21"),
                      centromere = c(123252373.5,93787431.5,90856062,50074452.5,
                                     48585285.5,60557102.5,60058972.5,61016889,
                                     45249872,43893383.5,53454152,39800499.5,
                                     35764400,17692000.5,17117352,19037747.5,
                                     36878628.5,25067566.5,18464134,28099979.5,
                                     26161912,10470308,15520235.5,11917946),
                      size = c(248956422,242193529,198295559,190214555,181538259,170805979,
                               159345973,156040895,145138636,138394717,135086622,133797422,
                               133275309,114364328,107043718,101991189,90338345,83257441,
                               80373285,64444167,58617616,57227415,50818468,46709983))

# set up bed files
bed_A1_f <- data_A_female[,2:4]
bed_A1_f$st1_conv <- NA
bed_A1_f$end1_conv <- NA
bed_A2_f <- data_A_female[,5:7]
bed_A2_f$st2_conv <- NA
bed_A2_f$end2_conv <- NA

bed_A1_m <- data_A_male[,2:4]
bed_A1_m$st1_conv <- NA
bed_A1_m$end1_conv <- NA
bed_A2_m <- data_A_male[,5:7]
bed_A2_m$st2_conv <- NA
bed_A2_m$end2_conv <- NA

bed_B1 <- data_B[,2:4]
bed_B1$st1_conv <- NA
bed_B1$end1_conv <- NA
bed_B2 <- data_B[,5:7]
bed_B2$st2_conv <- NA
bed_B2$end2_conv <- NA

bed_C1 <- data_C[,2:4]
bed_C1$st1_conv <- NA
bed_C1$end1_conv <- NA
bed_C2 <- data_C[,5:7]
bed_C2$st2_conv <- NA
bed_C2$end2_conv <- NA



for(data in c('data_A_female', 'data_A_male', 'data_B', 'data_C')) {
  
  # set the variables 
  df <- get(data)
  if (data == "data_A_female"){
    df_bed1 <- get("bed_A1_f")
    df_bed2 <- get("bed_A2_f")
  }
  if (data == "data_A_male"){
    df_bed1 <- get("bed_A1_m")
    df_bed2 <- get("bed_A2_m")
  }
  if (data == "data_B"){
    df_bed1 <- get("bed_B1")
    df_bed2 <- get("bed_B2")
  }
  if (data == "data_C"){
    df_bed1 <- get("bed_C1")
    df_bed2 <- get("bed_C2")
  }
  
  ##################################
  # linear conversion bed 1 (chrA) #
  ##################################
  for (i in 1:nrow(df)){
  
    # determine if bin is p or q arm
    df$chrA[i]
    pq <- which(df[i,10:13] == df$chrA[i])
    pq <- pq + 10
    arm = df[i,pq]
    
    # set new length based on p or q arm
    if (arm == "p_arm"){
      
      new_min = 0
      new_max = 1000
      old_min = 0
      old_max = chrInf$centromere[which(chrInf$chrom == df$chrA[i])]
      
      } else {
  
      new_min = 1000
      new_max = 2000
      old_min = chrInf$centromere[which(chrInf$chrom == df$chrA[i])]
      old_max = chrInf$size[which(chrInf$chrom == df$chrA[i])]
    
    }
    
    # calculate new length based on p or q arm
    #start
    value_to_change = as.numeric(df$st1[i])
    new_value = ( (value_to_change - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min
    df_bed1$st1_conv[i] <- new_value
    #end
    value_to_change = as.numeric(df$end1[i])
    new_value = ( (value_to_change - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min
    df_bed1$end1_conv[i] <- new_value
 
  }
  
  # save "bed" file as own dataframe
  assign(paste0(data,"_bed1"),df_bed1)
  
  
  ##################################
  # linear conversion bed 2 (chrB) #
  ##################################
  for (i in 1:nrow(df)){
    
    # determine if bin is p or q arm
    df$chrB[i]
    pq <- which(df[i,10:13] == df$chrB[i])
    pq <- pq + 10
    arm = df[i,pq]
    
    # set new length based on p or q arm
    if (arm == "p_arm"){
      
      new_min = 0
      new_max = 1000
      old_min = 0
      old_max = chrInf$centromere[which(chrInf$chrom == df$chrB[i])]
      
    } else {
      
      new_min = 1000
      new_max = 2000
      old_min = chrInf$centromere[which(chrInf$chrom == df$chrB[i])]
      old_max = chrInf$size[which(chrInf$chrom == df$chrB[i])]
      
    }
    
    # calculate new length based on p or q arm
    #start
    value_to_change = as.numeric(df$st2[i])
    new_value = ( (value_to_change - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min
    df_bed2$st2_conv[i] <- new_value
    #end
    value_to_change = as.numeric(df$end2[i])
    new_value = ( (value_to_change - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min
    df_bed2$end2_conv[i] <- new_value
    
  }
  
  # save "bed" file as own dataframe
  assign(paste0(data,"_bed2"),df_bed2)
  
}


# save entire workspace image with all BED files
filename <- paste0(output,"/pqarm_bedfiles_all.Rdata")
save.image(filename)

