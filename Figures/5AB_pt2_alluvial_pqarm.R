################################################################################
#  Circos ("alluvial") plot of gonosome interactions for p and q arms
################################################################################
library(circlize)
library(tidyr)
library(dplyr)


#______Read in arguments________________________________________________________
# running with a scheduler #
options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

output <- args[1]

load(paste0(output,"/pqarm_bedfiles_all.Rdata"))




#_____data processing___________________________________________________________

# remove old chromosomal length of bin, and keep only the unified length
data_A_female_bed1 <- data_A_female_bed1 %>% select(-c("st1","end1"))
data_A_female_bed2 <- data_A_female_bed2 %>% select(-c("st2","end2"))

data_A_male_bed1 <- data_A_male_bed1 %>% select(-c("st1","end1"))
data_A_male_bed2 <- data_A_male_bed2 %>% select(-c("st2","end2"))

data_B_bed1 <- data_B_bed1 %>% select(-c("st1","end1"))
data_B_bed2 <- data_B_bed2 %>% select(-c("st2","end2"))

data_C_bed1 <- data_C_bed1 %>% select(-c("st1","end1"))
data_C_bed2 <- data_C_bed2 %>% select(-c("st2","end2"))


#######################################
#  Plot A: chromosome X vs autosomes  #
#######################################
# rename link to either "autosome" or "X"

autosomes = c(paste0("chr",seq(1:22)))

# bed A1 - female
for (i in 1:nrow(data_A_female_bed1)){

  if (length(grep(data_A_female_bed1[i,1], autosomes) > 0)) { data_A_female_bed1[i,1] = "autosome" }
  if (length(grep(data_A_female_bed1[i,1], "chrX") > 0)) { data_A_female_bed1[i,1] = "X" }
  
}

# bed A2 - female
for (i in 1:nrow(data_A_female_bed2)){
  
  if (length(grep(data_A_female_bed2[i,1], autosomes) > 0)) { data_A_female_bed2[i,1] = "autosome" }
  if (length(grep(data_A_female_bed2[i,1], "chrX") > 0)) { data_A_female_bed2[i,1] = "X" }
  
}
  
# bed A1 - male
for (i in 1:nrow(data_A_male_bed1)){
  
  if (length(grep(data_A_male_bed1[i,1], autosomes) > 0)) { data_A_male_bed1[i,1] = "autosome" }
  if (length(grep(data_A_male_bed1[i,1], "chrX") > 0)) { data_A_male_bed1[i,1] = "X" }
  
}

# bed A2 - male #
for (i in 1:nrow(data_A_male_bed2)){
  
  if (length(grep(data_A_male_bed2[i,1], autosomes) > 0)) { data_A_male_bed2[i,1] = "autosome" }
  if (length(grep(data_A_male_bed2[i,1], "chrX") > 0)) { data_A_male_bed2[i,1] = "X" }
  
}
 

#######################################
#  Plot B: chromosome Y vs autosomes  #
#######################################
# rename link to either "autosome" or "Y"

autosomes = c(paste0("chr",seq(1:22)))

# bed B1
for (i in 1:nrow(data_B_bed1)){
  
  if (length(grep(data_B_bed1[i,1], autosomes) > 0)) { data_B_bed1[i,1] = "autosome" }
  if (length(grep(data_B_bed1[i,1], "chrY") > 0)) { data_B_bed1[i,1] = "Y" }
  
}

# bed B2 #
for (i in 1:nrow(data_B_bed2)){
  
  if (length(grep(data_B_bed2[i,1], autosomes) > 0)) { data_B_bed2[i,1] = "autosome" }
  if (length(grep(data_B_bed2[i,1], "chrY") > 0)) { data_B_bed2[i,1] = "Y" }
  
}  


##########################################
#  Plot C: chromosome Y vs chromosome X  #
##########################################
# rename link to either "X" or "Y"

# bed C1 #
for (i in 1:nrow(data_C_bed1)){
  
  if (length(grep(data_C_bed1[i,1], "chrX") > 0)) { data_C_bed1[i,1] = "X" }
  if (length(grep(data_C_bed1[i,1], "chrY") > 0)) { data_C_bed1[i,1] = "Y" }
  
}

# bed C2 #
for (i in 1:nrow(data_C_bed2)){
  
  if (length(grep(data_C_bed2[i,1], "chrX") > 0)) { data_C_bed2[i,1] = "X" }
  if (length(grep(data_C_bed2[i,1], "chrY") > 0)) { data_C_bed2[i,1] = "Y" }
  
}  




#_____generating and saving plots_______________________________________________

rev_x = function(x, sector.index = NULL) {
  if(!is.null(sector.index)) set.current.cell(sector.index, 1)
  xrange = CELL_META$xlim
  xrange[2] - x + xrange[1]
}


################################################
#  Plot A - female: chromosome X vs autosomes  #
################################################
filename <- "pqarm_Chromosome_X_vs_Autosomes_female.pdf"

# clear plots and initiate pdf printer
circos.clear()
pdf(paste0(output,filename), width = 10, height = 10)


      chromosomes <- data.frame(chr= c("autosome","X"), start= c(0,0), end= c(2000,2000))  
      
      circos.par("start.degree" = 165,"gap.degree" = c(30,30))
      circos.genomicInitialize(chromosomes, plotType = NULL)
  
      circos.track(ylim = c(0, 1),
                   panel.fun = function(x, y) {
                     major.by = seq(0, 2000, by = 250)
                     if(CELL_META$sector.index == "autosome") {
                       circos.axis(major.at = major.by, labels = paste0(major.by, "bp"),labels.cex = 0.4 * par("cex"))
                     } else {
                       circos.axis(major.at = rev_x(major.by), labels = paste0(major.by, "bp"),labels.cex = 0.4 * par("cex"))
                     }
                   }
      )
      
      
      {
      circos.text(500, 0.5, "P", sector.index = chromosomes$chr[1])
      circos.text(1000, 0.5, "|", sector.index = chromosomes$chr[1])
      circos.text(1500, 0.5, "Q", sector.index = chromosomes$chr[1])
      
      circos.text(rev_x(500), 0.5, "P", sector.index = chromosomes$chr[2], facing = "outside")
      circos.text(rev_x(1000), 0.5, "|", sector.index = chromosomes$chr[2], facing = "outside")
      circos.text(rev_x(1500), 0.5, "Q", sector.index = chromosomes$chr[2], facing = "outside")
      }

      
      # links
      for (i in 1:nrow(data_A_female_bed1)){
        
        circos.link(data_A_female_bed1[i,1],
                    if(data_A_female_bed1[i,1] == chromosomes$chr[1]) {c(data_A_female_bed1[i,2],data_A_female_bed1[i,3])} else {rev_x(c(data_A_female_bed1[i,2],data_A_female_bed1[i,3]))},
                    data_A_female_bed2[i,1],
                    if(data_A_female_bed2[i,1] == chromosomes$chr[1]) {c(data_A_female_bed2[i,2],data_A_female_bed2[i,3])} else {rev_x(c(data_A_female_bed2[i,2],data_A_female_bed2[i,3]))},
                    col = "#cc0000BF", border = NA
                    )
      }
      
      
# close pdf printer
dev.off()


##############################################
#  Plot A - male: chromosome X vs autosomes  #
##############################################
filename <- "pqarm_Chromosome_X_vs_Autosomes_male.pdf"

# clear plots and initiate pdf printer
circos.clear()
pdf(paste0(output,filename), width = 10, height = 10)


      chromosomes <- data.frame(chr= c("autosome","X"), start= c(0,0), end= c(2000,2000))  
      
      circos.par("start.degree" = 165,"gap.degree" = c(30,30))
      circos.genomicInitialize(chromosomes, plotType = NULL)
      
      circos.track(ylim = c(0, 1),
                   panel.fun = function(x, y) {
                     major.by = seq(0, 2000, by = 250)
                     if(CELL_META$sector.index == "autosome") {
                       circos.axis(major.at = major.by, labels = paste0(major.by, "bp"),labels.cex = 0.4 * par("cex"))
                     } else {
                       circos.axis(major.at = rev_x(major.by), labels = paste0(major.by, "bp"),labels.cex = 0.4 * par("cex"))
                     }
                   }
      )
      
      
      {
        circos.text(500, 0.5, "P", sector.index = chromosomes$chr[1])
        circos.text(1000, 0.5, "|", sector.index = chromosomes$chr[1])
        circos.text(1500, 0.5, "Q", sector.index = chromosomes$chr[1])
        
        circos.text(rev_x(500), 0.5, "P", sector.index = chromosomes$chr[2], facing = "outside")
        circos.text(rev_x(1000), 0.5, "|", sector.index = chromosomes$chr[2], facing = "outside")
        circos.text(rev_x(1500), 0.5, "Q", sector.index = chromosomes$chr[2], facing = "outside")
      }
      
      
      # links
      for (i in 1:nrow(data_A_male_bed1)){
        
        circos.link(data_A_male_bed1[i,1],
                    if(data_A_male_bed1[i,1] == chromosomes$chr[1]) {c(data_A_male_bed1[i,2],data_A_male_bed1[i,3])} else {rev_x(c(data_A_male_bed1[i,2],data_A_male_bed1[i,3]))},
                    data_A_male_bed2[i,1],
                    if(data_A_male_bed2[i,1] == chromosomes$chr[1]) {c(data_A_male_bed2[i,2],data_A_male_bed2[i,3])} else {rev_x(c(data_A_male_bed2[i,2],data_A_male_bed2[i,3]))},
                    col = "#cc0000BF", border = NA
        )
      }


# close pdf printer
dev.off()


#######################################
#  Plot B: chromosome Y vs autosomes  #
#######################################
filename <- "pqarm_Chromosome_Y_vs_Autosomes.pdf"

# clear plots and initiate pdf printer
circos.clear()
pdf(paste0(output,filename), width = 10, height = 10)


      chromosomes <- data.frame(chr= c("autosome","Y"), start= c(0,0), end= c(2000,2000))  
      
      circos.par("start.degree" = 165,"gap.degree" = c(30,30))
      circos.genomicInitialize(chromosomes, plotType = NULL)
      
      circos.track(ylim = c(0, 1),
                   panel.fun = function(x, y) {
                     major.by = seq(0, 2000, by = 250)
                     if(CELL_META$sector.index == "autosome") {
                       circos.axis(major.at = major.by, labels = paste0(major.by, "bp"),labels.cex = 0.4 * par("cex"))
                     } else {
                       circos.axis(major.at = rev_x(major.by), labels = paste0(major.by, "bp"),labels.cex = 0.4 * par("cex"))
                     }
                   }
      )
      
      
      {
        circos.text(500, 0.5, "P", sector.index = chromosomes$chr[1])
        circos.text(1000, 0.5, "|", sector.index = chromosomes$chr[1])
        circos.text(1500, 0.5, "Q", sector.index = chromosomes$chr[1])
        
        circos.text(rev_x(500), 0.5, "P", sector.index = chromosomes$chr[2], facing = "outside")
        circos.text(rev_x(1000), 0.5, "|", sector.index = chromosomes$chr[2], facing = "outside")
        circos.text(rev_x(1500), 0.5, "Q", sector.index = chromosomes$chr[2], facing = "outside")
      }
      
      
      # links
      for (i in 1:nrow(data_B_bed1)){
        
        circos.link(data_B_bed1[i,1],
                    if(data_B_bed1[i,1] == chromosomes$chr[1]) {c(data_B_bed1[i,2],data_B_bed1[i,3])} else {rev_x(c(data_B_bed1[i,2],data_B_bed1[i,3]))},
                    data_B_bed2[i,1],
                    if(data_B_bed2[i,1] == chromosomes$chr[1]) {c(data_B_bed2[i,2],data_B_bed2[i,3])} else {rev_x(c(data_B_bed2[i,2],data_B_bed2[i,3]))},
                    col = "#cc0000BF", border = NA
        )
      }


# close pdf printer
dev.off()


##########################################
#  Plot C: chromosome X vs chromosome Y  #
##########################################
filename <- "pqarm_Chromosome_X_vs_Chromosome_Y.pdf"

# clear plots and initiate pdf printer
circos.clear()
pdf(paste0(output,filename), width = 10, height = 10)

      
      chromosomes <- data.frame(chr= c("X","Y"), start= c(0,0), end= c(2000,2000))  
      
      circos.par("start.degree" = 165,"gap.degree" = c(30,30))
      circos.genomicInitialize(chromosomes, plotType = NULL)
      
      circos.track(ylim = c(0, 1),
                   panel.fun = function(x, y) {
                     major.by = seq(0, 2000, by = 250)
                     if(CELL_META$sector.index == "X") {
                       circos.axis(major.at = major.by, labels = paste0(major.by, "bp"),labels.cex = 0.4 * par("cex"))
                     } else {
                       circos.axis(major.at = rev_x(major.by), labels = paste0(major.by, "bp"),labels.cex = 0.4 * par("cex"))
                     }
                   }
      )
      
      
      {
        circos.text(500, 0.5, "P", sector.index = chromosomes$chr[1])
        circos.text(1000, 0.5, "|", sector.index = chromosomes$chr[1])
        circos.text(1500, 0.5, "Q", sector.index = chromosomes$chr[1])
        
        circos.text(rev_x(500), 0.5, "P", sector.index = chromosomes$chr[2], facing = "outside")
        circos.text(rev_x(1000), 0.5, "|", sector.index = chromosomes$chr[2], facing = "outside")
        circos.text(rev_x(1500), 0.5, "Q", sector.index = chromosomes$chr[2], facing = "outside")
      }
      
      
      # links
      for (i in 1:nrow(data_C_bed1)){
        
        circos.link(data_C_bed1[i,1],
                    if(data_C_bed1[i,1] == chromosomes$chr[1]) {c(data_C_bed1[i,2],data_C_bed1[i,3])} else {rev_x(c(data_C_bed1[i,2],data_C_bed1[i,3]))},
                    data_C_bed2[i,1],
                    if(data_C_bed2[i,1] == chromosomes$chr[1]) {c(data_C_bed2[i,2],data_C_bed2[i,3])} else {rev_x(c(data_C_bed2[i,2],data_C_bed2[i,3]))},
                    col = "#cc0000BF", border = NA
        )
      }


# close pdf printer
dev.off()





